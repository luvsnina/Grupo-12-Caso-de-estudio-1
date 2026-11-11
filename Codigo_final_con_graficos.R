# ==============================================================================
# USACH - Estadística Aplicada | Caso 1: regresión
# Predicción de TEY (energía por período horario, MWh) - Turbina de gas UCI
# Integrantes: 
#Maria Ignacia Cordova 
#Martín Peñailillo Salazar 
#Millaray Rodríguez Oses
#Vicente Urra Cisternas


# PASO 1. CONFIGURACIÓN Y PAQUETES ----------------------------------------------
SEMILLA <- 2026L
ARBOLES_RF <- 150L
PROFUNDIDAD_MAX_RF <- 14L
nucleos <- parallel::detectCores(logical = TRUE)
HILOS <- if (length(nucleos) != 1L || is.na(nucleos)) 1L else as.integer(min(2L, max(1L, nucleos)))
set.seed(SEMILLA)

directorio_script <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  ruta_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  if (length(ruta_arg)) return(dirname(normalizePath(ruta_arg[1], mustWork = TRUE)))
  rutas_source <- lapply(sys.frames(), function(f) f$ofile)
  rutas_source <- Filter(function(x) !is.null(x) && nzchar(x), rutas_source)
  if (length(rutas_source)) {
    return(dirname(normalizePath(tail(rutas_source, 1)[[1]], mustWork = TRUE)))
  }
  getwd()
}
raiz_proyecto <- directorio_script()
carpeta_datos <- file.path(raiz_proyecto, "dataset")
carpeta_salida <- file.path(raiz_proyecto, "resultados_finales")
if (!dir.exists(carpeta_salida)) dir.create(carpeta_salida, recursive = TRUE)

paquetes <- c("readr", "dplyr", "ggplot2", "glmnet", "mgcv", "rpart",
              "rpart.plot", "ranger", "xgboost", "car")
faltantes <- paquetes[!vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltantes)) {
  stop(paste0("Faltan paquetes. Ejecuta una sola vez: install.packages(c(",
              paste(shQuote(faltantes), collapse = ", "), "))"))
}
suppressPackageStartupMessages(library(mgcv))

salvar_csv <- function(x, archivo) {
  readr::write_csv(x, file.path(carpeta_salida, archivo), na = "")
}
salvar_grafico <- function(grafico, archivo, ancho = 9, alto = 5.5) {
  ggplot2::ggsave(file.path(carpeta_salida, archivo), plot = grafico,
                  width = ancho, height = alto, dpi = 200)
}
metricas <- function(y, yhat) {
  y <- as.numeric(y)
  yhat <- as.numeric(yhat)
  if (!length(y) || length(y) != length(yhat) || any(!is.finite(yhat))) {
    stop("Predicciones vacías, no finitas o con longitud incorrecta.")
  }
  sst <- sum((y - mean(y))^2)
  data.frame(RMSE = sqrt(mean((y - yhat)^2)),
             MAE = mean(abs(y - yhat)),
             R2 = if (sst > 0) 1 - sum((y - yhat)^2) / sst else NA_real_)
}

# PASO 2. IMPORTACIÓN ESTRICTA DE LOS CINCO CSV -------------------------------
anios_esperados <- 2011:2015
archivos_esperados <- paste0("gt_", anios_esperados, ".csv")
columnas_originales <- c("AT", "AP", "AH", "AFDP", "GTEP", "TIT", "TAT",
                         "TEY", "CDP", "CO", "NOX")
predictores <- c("AT", "AP", "AH", "AFDP", "GTEP", "TIT", "TAT", "CDP")
if (!dir.exists(carpeta_datos)) stop("Falta la carpeta 'dataset' junto al script.")
archivos_encontrados <- list.files(carpeta_datos, pattern = "\\.csv$",
                                   full.names = FALSE, ignore.case = TRUE)
if (length(archivos_encontrados) != 5L ||
    !setequal(archivos_encontrados, archivos_esperados)) {
  stop("La carpeta dataset debe contener EXACTAMENTE gt_2011.csv a gt_2015.csv.")
}

lecturas <- lapply(anios_esperados, function(anio) {
  ruta <- file.path(carpeta_datos, paste0("gt_", anio, ".csv"))
  df <- readr::read_csv(ruta, show_col_types = FALSE,
                        col_types = readr::cols(.default = readr::col_double()))
  if (nrow(readr::problems(df))) stop("Problemas de lectura en: ", ruta)
  if (!setequal(names(df), columnas_originales) || ncol(df) != 11L) {
    stop("Columnas no esperadas en: ", ruta)
  }
  df <- df[, columnas_originales]
  if (anyNA(df) || any(!is.finite(as.matrix(df)))) {
    stop("Valores NA/Inf en: ", ruta, ". Se requiere auditoría antes de modelar.")
  }
  df$Year <- as.integer(anio)
  df$orden_en_anio <- seq_len(nrow(df))
  df
})
datos_completos <- dplyr::bind_rows(lecturas)
rutas_originales <- file.path(carpeta_datos, archivos_esperados)
md5_original <- unname(tools::md5sum(rutas_originales))
salvar_csv(data.frame(archivo = archivos_esperados, md5_original = md5_original),
           "00_huella_archivos_originales.csv")
stopifnot(identical(sort(unique(datos_completos$Year)), anios_esperados))

# PASO 3. AUDITORÍA DE REGISTROS Y DUPLICADOS --------------------------------
auditoria <- dplyr::bind_rows(lapply(lecturas, function(df) {
  data.frame(Year = unique(df$Year), registros = nrow(df),
             duplicados_exactos = if (unique(df$Year) == 2015L) NA_integer_
             else sum(duplicated(df[, columnas_originales])),
             faltantes = sum(is.na(df)),
             TEY_min_MWh = if (unique(df$Year) == 2015L) NA_real_ else min(df$TEY),
             TEY_max_MWh = if (unique(df$Year) == 2015L) NA_real_ else max(df$TEY),
             AT_min_C = if (unique(df$Year) == 2015L) NA_real_ else min(df$AT),
             AT_max_C = if (unique(df$Year) == 2015L) NA_real_ else max(df$AT))
}))
salvar_csv(auditoria, "01_auditoria_por_anio.csv")
cat("\nAUDITORÍA DE ARCHIVOS\n")
print(auditoria)

# PASO 4. EXPLORACIÓN: SOLO AÑOS ANTERIORES A TEST ---------------------------
datos_2011_2014 <- dplyr::filter(datos_completos, Year <= 2014)
resumen_eda <- dplyr::bind_rows(lapply(c(predictores, "TEY"), function(v) {
  z <- datos_2011_2014[[v]]
  data.frame(Variable = v, minimo = min(z), p25 = unname(quantile(z, .25)),
             media = mean(z), mediana = median(z), p75 = unname(quantile(z, .75)),
             maximo = max(z), sd = sd(z))
}))
salvar_csv(resumen_eda, "02_resumen_descriptivo_train.csv")

extremos_train <- dplyr::bind_rows(lapply(c(predictores, "TEY"), function(v) {
  z <- datos_2011_2014[[v]]
  limites <- as.numeric(stats::quantile(z, c(.25, .75), names = FALSE))
  amplitud <- limites[2] - limites[1]
  li <- limites[1] - 3 * amplitud
  ls <- limites[2] + 3 * amplitud
  data.frame(Variable = v, limite_inferior = li, limite_superior = ls,
             n_extremos_3IQR = sum(z < li | z > ls))
}))
salvar_csv(extremos_train, "02b_revision_extremos_train.csv")

correlaciones <- data.frame(
  Variable = predictores,
  Pearson_con_TEY = vapply(predictores, function(v) {
    cor(datos_2011_2014[[v]], datos_2011_2014$TEY)
  }, numeric(1))
)
salvar_csv(correlaciones, "03_correlaciones_train.csv")

p_eda <- ggplot2::ggplot(datos_2011_2014, ggplot2::aes(x = TEY)) +
  ggplot2::geom_histogram(bins = 40, fill = "#355C7D", color = "white") +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "Distribución de TEY - entrenamiento 2011-2014",
                x = "Energía TEY (MWh por período horario)", y = "Frecuencia")
salvar_grafico(p_eda, "01_histograma_tey_train.png")

# PASO 5. VARIABLES Y DEFINICIÓN DE LA TAREA ---------------------------------
datos_modelo <- datos_completos[, c("Year", "orden_en_anio", "TEY", predictores)]
datos_modelo$id_registro <- seq_len(nrow(datos_modelo))

# PASO 6. PARTICIONES TEMPORALES, BASELINE Y CALIBRACIÓN -----------------------
train_total <- dplyr::filter(datos_modelo, Year <= 2014)
stopifnot(all(predictores %in% names(datos_modelo)), length(predictores) == 8L)
test_2015 <- dplyr::filter(datos_modelo, Year == 2015)
corte_2014 <- floor(sum(train_total$Year == 2014) * .75)
desarrollo <- dplyr::filter(train_total,
                            Year < 2014 | (Year == 2014 & orden_en_anio <= corte_2014))
calibracion <- dplyr::filter(train_total,
                             Year == 2014 & orden_en_anio > corte_2014)
stopifnot(nrow(calibracion) > 500L,
          length(intersect(desarrollo$id_registro, calibracion$id_registro)) == 0L,
          length(intersect(train_total$id_registro, test_2015$id_registro)) == 0L)

folds <- list(
  list(nombre = "F1_2012", n = 1L,
       train = dplyr::filter(desarrollo, Year == 2011),
       valid = dplyr::filter(desarrollo, Year == 2012)),
  list(nombre = "F2_2013", n = 2L,
       train = dplyr::filter(desarrollo, Year <= 2012),
       valid = dplyr::filter(desarrollo, Year == 2013)),
  list(nombre = "F3_2014_inicial", n = 3L,
       train = dplyr::filter(desarrollo, Year <= 2013),
       valid = dplyr::filter(desarrollo, Year == 2014))
)
for (f in folds) {
  stopifnot(max(f$train$id_registro) < min(f$valid$id_registro),
            max(f$train$Year) < min(f$valid$Year),
            !any(f$train$Year == 2015), !any(f$valid$Year == 2015),
            length(intersect(f$train$id_registro, f$valid$id_registro)) == 0L)
}
resumen_particiones <- dplyr::bind_rows(lapply(folds, function(f) {
  data.frame(Fold = f$nombre, hasta_anio_train = max(f$train$Year),
             anio_validacion = unique(f$valid$Year),
             n_train = nrow(f$train), n_validacion = nrow(f$valid))
}))
salvar_csv(resumen_particiones, "04_folds_temporales.csv")
cat("\nFOLDS TEMPORALES (sin datos futuros)\n")
print(resumen_particiones)

# PASO 6B. MATRIZ X-X Y DIAGNÓSTICO DE REDUNDANCIA ---------------------------
UMBRAL_COR <- 0.90
matriz_xx <- stats::cor(desarrollo[, predictores, drop = FALSE])
salvar_csv(data.frame(Variable = rownames(matriz_xx), matriz_xx,
                      check.names = FALSE), "03b_matriz_correlacion_XX_desarrollo.csv")
pares_xx <- as.data.frame(as.table(matriz_xx))
names(pares_xx) <- c("Variable_1", "Variable_2", "Pearson_r")
pares_xx <- pares_xx[as.integer(factor(pares_xx$Variable_1, levels = predictores)) <
                       as.integer(factor(pares_xx$Variable_2, levels = predictores)), ]
pares_xx$Correlacion_absoluta <- abs(pares_xx$Pearson_r)
pares_xx$Alerta <- pares_xx$Correlacion_absoluta >= UMBRAL_COR
pares_xx <- pares_xx[order(-pares_xx$Correlacion_absoluta), ]
salvar_csv(pares_xx, "03c_pares_correlacion_XX.csv")

p_cor <- ggplot2::ggplot(
  transform(as.data.frame(as.table(matriz_xx)),
            Var1 = factor(Var1, levels = predictores),
            Var2 = factor(Var2, levels = rev(predictores))),
  ggplot2::aes(x = Var1, y = Var2, fill = Freq)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Freq)), size = 3) +
  ggplot2::scale_fill_gradient2(low = "#286B81", mid = "white",
                                high = "#B4534B", midpoint = 0,
                                limits = c(-1, 1), name = "r Pearson") +
  ggplot2::coord_equal() + ggplot2::theme_minimal() +
  ggplot2::labs(title = "Correlaciones entre las ocho predictoras (X-X)",
                subtitle = "Solo desarrollo 2011-2014 inicial; sin calibración ni test",
                x = NULL, y = NULL)
salvar_grafico(p_cor, "02_heatmap_correlacion_XX.png", ancho = 9, alto = 7)

conservar_predictores <- predictores
stopifnot(length(conservar_predictores) == 8L)
salvar_csv(data.frame(Variable = predictores,
                      Conservada_como_candidata = TRUE),
           "03e_variables_conservadas.csv")
cat("\nSe conservan las OCHO variables predictoras como candidatas.\n")

lambdas <- 10^seq(2, -4, length.out = 11L)
form_lineal <- stats::as.formula(paste("TEY ~", paste(predictores, collapse = " + ")))

ajustar_modelo <- function(nombre, config, entreno, semilla = SEMILLA) {
  set.seed(as.integer(semilla))
  variables_ajuste <- predictores
  d <- as.data.frame(entreno[, c("TEY", variables_ajuste), drop = FALSE])
  form_ajuste <- stats::reformulate(variables_ajuste, response = "TEY")
  if (nombre == "Baseline") {
    return(list(tipo = nombre, modelo = mean(d$TEY), config = config))
  }
  if (nombre == "OLS") {
    return(list(tipo = nombre, modelo = stats::lm(form_ajuste, data = d),
                config = config, variables = variables_ajuste))
  }
  if (nombre %in% c("Ridge", "Lasso", "Elastic Net")) {
    alfa <- if (nombre == "Lasso") 1 else
      if (nombre == "Ridge") 0 else as.numeric(config$alpha)
    m <- glmnet::glmnet(x = as.matrix(d[, variables_ajuste, drop = FALSE]), y = d$TEY,
                        family = "gaussian", alpha = alfa,
                        lambda = lambdas, standardize = TRUE)
    return(list(tipo = nombre, modelo = m, config = config,
                variables = variables_ajuste))
  }
  if (nombre == "GAM") {
    f <- stats::as.formula(paste0("TEY ~ ", paste(
      sprintf("s(%s, k=%d)", predictores, as.integer(config$k)), collapse = " + ")))
    m <- mgcv::gam(f, data = d, method = "REML", select = TRUE)
    return(list(tipo = nombre, modelo = m, config = config))
  }
  if (nombre == "Arbol") {
    ctrl <- rpart::rpart.control(cp = config$cp,
                                 maxdepth = as.integer(config$maxdepth),
                                 minsplit = 20L, xval = 0L)
    m <- rpart::rpart(form_lineal, data = d, method = "anova", control = ctrl)
    return(list(tipo = nombre, modelo = m, config = config))
  }
  if (nombre == "Random Forest") {
    m <- ranger::ranger(formula = form_lineal, data = d,
                        num.trees = ARBOLES_RF, mtry = as.integer(config$mtry),
                        min.node.size = as.integer(config$min.node.size),
                        max.depth = PROFUNDIDAD_MAX_RF,
                        importance = "none", splitrule = "variance",
                        seed = as.integer(semilla), num.threads = HILOS,
                        verbose = FALSE)
    return(list(tipo = nombre, modelo = m, config = config))
  }
  if (nombre == "XGBoost") {
    params <- list(booster = "gbtree", objective = "reg:squarederror",
                   eta = as.numeric(config$eta),
                   max_depth = as.integer(config$max_depth),
                   subsample = .85, colsample_bytree = .85, nthread = HILOS)
    dm <- xgboost::xgb.DMatrix(as.matrix(d[, predictores]), label = d$TEY)
    m <- xgboost::xgb.train(params = params, data = dm,
                            nrounds = as.integer(config$nrounds), verbose = 0L)
    return(list(tipo = nombre, modelo = m, config = config))
  }
  if (nombre == "Averaging") {
    modelos <- list(
      rf = ajustar_modelo("Random Forest", config$rf, entreno, semilla),
      xgb = ajustar_modelo("XGBoost", config$xgb, entreno, semilla),
      gam = ajustar_modelo("GAM", config$gam, entreno, semilla))
    return(list(tipo = nombre, modelo = modelos, config = config))
  }
  stop("Modelo desconocido: ", nombre)
}

predecir_modelo <- function(obj, nuevos) {
  pred_en_modelo <- if (is.null(obj$variables)) predictores else obj$variables
  d <- as.data.frame(nuevos[, pred_en_modelo, drop = FALSE])
  tipo <- obj$tipo
  if (tipo == "Baseline") return(rep(obj$modelo, nrow(d)))
  if (tipo %in% c("OLS", "GAM", "Arbol")) {
    return(as.numeric(stats::predict(obj$modelo, newdata = d)))
  }
  if (tipo %in% c("Ridge", "Lasso", "Elastic Net")) {
    return(as.numeric(stats::predict(obj$modelo, newx = as.matrix(d),
                                     s = as.numeric(obj$config$lambda))))
  }
  if (tipo == "Random Forest") {
    return(as.numeric(predict(obj$modelo, data = d)$predictions))
  }
  if (tipo == "XGBoost") {
    return(as.numeric(predict(obj$modelo, newdata = as.matrix(d))))
  }
  if (tipo == "Averaging") {
    p <- obj$config$pesos
    return(as.numeric(p[1] * predecir_modelo(obj$modelo$rf, nuevos) +
                        p[2] * predecir_modelo(obj$modelo$xgb, nuevos) +
                        p[3] * predecir_modelo(obj$modelo$gam, nuevos)))
  }
  stop("Modelo no programado en predecir_modelo().")
}

# PASO 7. OLS: COEFICIENTES, INTERVALOS Y MULTICOLINEALIDAD --------------------
ols_exploratorio <- ajustar_modelo("OLS", list(), desarrollo)$modelo
ols_coef <- as.data.frame(summary(ols_exploratorio)$coefficients)
ols_coef$Variable <- rownames(ols_coef)
ols_ic <- as.data.frame(stats::confint(ols_exploratorio, level = .95))
ols_ic$Variable <- rownames(ols_ic)
salvar_csv(dplyr::left_join(ols_coef, ols_ic, by = "Variable"),
           "05_ols_coeficientes_ic95.csv")
vif_ols <- car::vif(ols_exploratorio)
salvar_csv(data.frame(Variable = names(vif_ols), VIF = as.numeric(vif_ols)),
           "06_vif_ols.csv")

cook <- stats::cooks.distance(ols_exploratorio)
indices_cook <- head(order(cook, decreasing = TRUE), 15L)
salvar_csv(data.frame(indice_desarrollo = indices_cook,
                      distancia_cook = cook[indices_cook],
                      umbral_orientativo = 4 / nrow(desarrollo)),
           "07_observaciones_influyentes_ols.csv")

df_vif <- data.frame(Variable = names(vif_ols), VIF = as.numeric(vif_ols))

# Filtrar o seleccionar las dos variables con colinealidad extrema (VIF > 100)
df_vif_altos <- data.frame(
  Variable = names(vif_ols),
  VIF = as.numeric(vif_ols)
) |> 
  dplyr::filter(VIF > 100 | Variable %in% c("CDP", "GTEP"))

# Asegurar orden visual: CDP a la izquierda, GTEP a la derecha
df_vif_altos$Variable <- factor(df_vif_altos$Variable, levels = c("CDP", "GTEP"))

p_vif_extremo <- ggplot2::ggplot(df_vif_altos, ggplot2::aes(x = Variable, y = VIF, fill = Variable)) +
  ggplot2::geom_col(width = 0.55) +
  ggplot2::geom_hline(yintercept = 10, linetype = "dashed", color = "#D35400", linewidth = 0.8) +
  ggplot2::geom_text(
    ggplot2::aes(label = gsub("\\.", ",", sprintf("%.1f", VIF))),
    vjust = -0.5,
    size = 4,
    fontface = "bold"
  ) +
  ggplot2::annotate(
    "text", x = 2.1, y = 14, 
    label = "Referencia orientativa: VIF = 10", 
    color = "#D35400", size = 3.3, hjust = 1
  ) +
  ggplot2::scale_fill_manual(values = c("CDP" = "#2C4D63", "GTEP" = "#009688")) +
  ggplot2::scale_y_continuous(limits = c(0, 260), breaks = seq(0, 250, by = 50)) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "none",
    plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
    plot.caption = ggplot2::element_text(size = 8.5, hjust = 0, color = "#555555"),
    panel.grid.major.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(face = "bold", size = 10)
  ) +
  ggplot2::labs(
    title = "Multicolinealidad en OLS — variables con VIF elevado",
    x = NULL,
    y = "Factor de inflación de la varianza (VIF)",
    caption = "VIF alto = redundancia entre predictores; no implica falla del sensor."
  )

salvar_grafico(p_vif_extremo, "02b_diagnostico_vif_extremo.png", ancho = 7, alto = 5)

# PASO 8. REGRESIONES REGULARIZADAS --------------------------------------------
grid_ridge <- data.frame(lambda = lambdas)
grid_lasso <- data.frame(lambda = lambdas)
grid_enet <- expand.grid(alpha = c(.25, .50, .75), lambda = lambdas)

# PASO 9. GAM ------------------------------------------------------------------
grid_gam <- data.frame(k = c(5L, 9L))

# PASO 10. ÁRROL Y RANDOM FOREST -----------------------------------------------
grid_arbol <- expand.grid(cp = c(.001, .01), maxdepth = c(5L, 10L))
grid_rf <- expand.grid(mtry = c(3L, 6L),
                       min.node.size = c(10L, 20L))

# PASO 11. XGBOOST -------------------------------------------------------------
grid_xgb <- data.frame(eta = c(.05, .10, .05, .10),
                       max_depth = c(3L, 3L, 5L, 5L),
                       nrounds = c(120L, 180L, 180L, 120L))

# PASO 12. VALIDACIÓN TEMPORAL Y SELECCIÓN DE MODELOS -------------------------
grillas <- list(
  "Baseline" = data.frame(config_id = 1L),
  "OLS" = data.frame(config_id = 1L),
  "Ridge" = grid_ridge,
  "Lasso" = grid_lasso,
  "Elastic Net" = grid_enet,
  "GAM" = grid_gam,
  "Arbol" = grid_arbol,
  "Random Forest" = grid_rf,
  "XGBoost" = grid_xgb
)
grillas <- lapply(grillas, function(g) {
  g$config_id <- seq_len(nrow(g))
  g
})
configuracion <- function(nombre, id) {
  g <- grillas[[nombre]]
  fila <- g[g$config_id == id, setdiff(names(g), "config_id"), drop = FALSE]
  if (nrow(fila) != 1L) stop("Configuración no encontrada.")
  as.list(fila)
}

resultados_grilla <- list()
contador <- 0L
for (nombre in names(grillas)) {
  cat("\nValidando: ", nombre, "\n", sep = "")
  for (id in grillas[[nombre]]$config_id) {
    cfg <- configuracion(nombre, id)
    for (f in folds) {
      if (nombre %in% c("Random Forest", "XGBoost")) {
        cat(sprintf("  Configuración %d/%d, bloque %s\n",
                    id, nrow(grillas[[nombre]]), f$nombre))
        flush.console()
      }
      objeto <- ajustar_modelo(nombre, cfg, f$train, semilla = SEMILLA + f$n)
      pred <- predecir_modelo(objeto, f$valid)
      m <- metricas(f$valid$TEY, pred)
      contador <- contador + 1L
      resultados_grilla[[contador]] <- data.frame(
        Modelo = nombre, config_id = id, Fold = f$nombre,
        RMSE = m$RMSE, MAE = m$MAE, R2 = m$R2)
    }
  }
}
resultados_grilla <- dplyr::bind_rows(resultados_grilla)
salvar_csv(resultados_grilla, "08_metricas_todas_las_configuraciones_cv.csv")
ranking_grilla <- resultados_grilla |>
  dplyr::group_by(Modelo, config_id) |>
  dplyr::summarise(RMSE_medio = mean(RMSE), sd_RMSE = stats::sd(RMSE),
                   MAE_medio = mean(MAE), .groups = "drop") |>
  dplyr::arrange(Modelo, RMSE_medio, sd_RMSE)
mejor_config <- ranking_grilla |>
  dplyr::group_by(Modelo) |>
  dplyr::slice_head(n = 1L) |>
  dplyr::ungroup()
salvar_csv(mejor_config, "09_hiperparametros_seleccionados.csv")
config_elegida <- function(nombre) {
  id <- mejor_config$config_id[mejor_config$Modelo == nombre]
  if (length(id) != 1L) stop("No hay selección única para: ", nombre)
  configuracion(nombre, id)
}

variables_por_fold <- dplyr::bind_rows(lapply(folds, function(f) {
  data.frame(Fold = f$nombre, Variable = predictores, Conservada = TRUE)
}))
salvar_csv(variables_por_fold, "09c_variables_conservadas_por_fold.csv")

oof_list <- list()
contador <- 0L
for (nombre in names(grillas)) {
  for (f in folds) {
    objeto <- ajustar_modelo(nombre, config_elegida(nombre), f$train,
                             semilla = SEMILLA + f$n)
    pred <- predecir_modelo(objeto, f$valid)
    contador <- contador + 1L
    oof_list[[contador]] <- data.frame(
      Modelo = nombre, Fold = f$nombre, Year = f$valid$Year,
      id_registro = f$valid$id_registro,
      observado = f$valid$TEY, predicho = pred)
  }
}
oof <- dplyr::bind_rows(oof_list)

pesos_averaging <- c(.50, .30, .20)
stopifnot(abs(sum(pesos_averaging) - 1) < 1e-10)
rf_oof <- oof |>
  dplyr::filter(Modelo == "Random Forest") |>
  dplyr::select(Fold, Year, id_registro, observado, pred_rf = predicho)
xgb_oof <- oof |>
  dplyr::filter(Modelo == "XGBoost") |>
  dplyr::select(Fold, id_registro, pred_xgb = predicho)
gam_oof <- oof |>
  dplyr::filter(Modelo == "GAM") |>
  dplyr::select(Fold, id_registro, pred_gam = predicho)
mezcla <- rf_oof |>
  dplyr::inner_join(xgb_oof, by = c("Fold", "id_registro")) |>
  dplyr::inner_join(gam_oof, by = c("Fold", "id_registro"))
stopifnot(nrow(mezcla) == nrow(rf_oof))
oof_avg <- data.frame(
  Modelo = "Averaging", Fold = mezcla$Fold, Year = mezcla$Year,
  id_registro = mezcla$id_registro, observado = mezcla$observado,
  predicho = pesos_averaging[1] * mezcla$pred_rf +
    pesos_averaging[2] * mezcla$pred_xgb +
    pesos_averaging[3] * mezcla$pred_gam)
oof <- dplyr::bind_rows(oof, oof_avg)
oof$residuo <- oof$observado - oof$predicho
salvar_csv(oof, "10_predicciones_fuera_de_muestra_cv.csv")

por_fold <- oof |>
  dplyr::group_by(Modelo, Fold) |>
  dplyr::summarise(n = dplyr::n(),
                   RMSE = sqrt(mean((observado - predicho)^2)),
                   MAE = mean(abs(observado - predicho)),
                   R2 = if (sum((observado - mean(observado))^2) > 0)
                     1 - sum((observado - predicho)^2) /
                     sum((observado - mean(observado))^2) else NA_real_,
                   .groups = "drop")
base_por_fold <- por_fold |>
  dplyr::filter(Modelo == "Baseline") |>
  dplyr::select(Fold, RMSE_baseline_fold = RMSE)
por_fold <- por_fold |>
  dplyr::left_join(base_por_fold, by = "Fold") |>
  dplyr::mutate(Mejora_vs_baseline_pct =
                  100 * (RMSE_baseline_fold - RMSE) / RMSE_baseline_fold)
salvar_csv(por_fold, "11_resultados_por_fold.csv")

tabla_validacion <- por_fold |>
  dplyr::group_by(Modelo) |>
  dplyr::summarise(RMSE_medio = mean(RMSE), sd_RMSE = stats::sd(RMSE),
                   MAE_medio = mean(MAE), R2_medio = mean(R2),
                   mejora_media_pct = mean(Mejora_vs_baseline_pct),
                   n_folds = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(RMSE_medio, sd_RMSE, MAE_medio)
salvar_csv(tabla_validacion, "12_TABLA_OFICIAL_VALIDACION_TEMPORAL.csv")
cat("\nTABLA DE VALIDACIÓN TEMPORAL\n")
print(tabla_validacion, n = Inf)

ganador <- tabla_validacion$Modelo[1]
cfg_averaging <- list(rf = config_elegida("Random Forest"),
                      xgb = config_elegida("XGBoost"),
                      gam = config_elegida("GAM"), pesos = pesos_averaging)
cfg_ganador <- if (ganador == "Averaging") cfg_averaging else config_elegida(ganador)
cat("\nSISTEMA SELECCIONADO POR CV: ", ganador, "\n", sep = "")

lasso_seleccionado <- ajustar_modelo("Lasso", config_elegida("Lasso"),
                                     desarrollo, semilla = SEMILLA + 55L)
coef_lasso <- as.matrix(stats::coef(lasso_seleccionado$modelo,
                                    s = config_elegida("Lasso")$lambda))
salvar_csv(data.frame(Variable = rownames(coef_lasso),
                      Coeficiente = as.numeric(coef_lasso[, 1])),
           "09b_coeficientes_lasso.csv")
df_lasso_plot <- data.frame(
  Variable = rownames(coef_lasso),
  Coeficiente = as.numeric(coef_lasso[, 1])
) |> 
  dplyr::filter(Variable != "(Intercept)")

df_lasso_plot$Variable <- factor(df_lasso_plot$Variable, 
                                 levels = c("AT", "AP", "AH", "AFDP", "GTEP", "TIT", "TAT", "CDP"))
df_lasso_plot$Color_Grupo <- ifelse(df_lasso_plot$Coeficiente >= 0, "Positivo", "Negativo")

n_activos <- sum(df_lasso_plot$Coeficiente != 0)
n_ceros <- sum(df_lasso_plot$Coeficiente == 0)

p_lasso <- ggplot2::ggplot(df_lasso_plot, ggplot2::aes(x = Variable, y = Coeficiente, fill = Color_Grupo)) +
  ggplot2::geom_col(width = 0.55) +
  ggplot2::geom_hline(yintercept = 0, color = "#7F8C8D", linewidth = 0.8) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = ifelse(Coeficiente == 0, "0", sprintf("%+.2f", Coeficiente)),
      vjust = ifelse(Coeficiente >= 0, -0.6, 1.4)
    ),
    size = 3.8
  ) +
  ggplot2::scale_fill_manual(values = c("Positivo" = "#009688", "Negativo" = "#D35400")) +
  ggplot2::scale_y_continuous(limits = c(-1.5, 8.2), breaks = seq(0, 8, by = 2)) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "none",
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.caption = ggplot2::element_text(size = 8.5, color = "#555555", hjust = 0),
    panel.grid.major.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(face = "bold", size = 10)
  ) +
  ggplot2::labs(
    title = sprintf("Coeficientes Lasso — %d activos y %d en cero", n_activos, n_ceros),
    x = NULL,
    y = "Coeficiente estimado reportado",
    caption = "No comparar magnitudes sin comprobar escalas/unidades; asociación predictiva, no efecto causal."
  )

salvar_grafico(p_lasso, "02c_coeficientes_lasso.png", ancho = 9, alto = 5.5)

arbol_seleccionado <- ajustar_modelo("Arbol", config_elegida("Arbol"),
                                     desarrollo, semilla = SEMILLA + 56L)
grDevices::png(file.path(carpeta_salida, "03_arbol_regresion.png"),
               width = 1400, height = 1100, res = 150)
rpart.plot::rpart.plot(arbol_seleccionado$modelo, main = "Árbol de regresión seleccionado por CV")
grDevices::dev.off()

gam_seleccionado <- ajustar_modelo("GAM", config_elegida("GAM"),
                                   desarrollo, semilla = SEMILLA + 57L)
grDevices::png(file.path(carpeta_salida, "04_efectos_gam.png"),
               width = 1600, height = 1200, res = 150)
plot(gam_seleccionado$modelo, pages = 1, shade = TRUE, rug = FALSE)
grDevices::dev.off()

p_comp <- ggplot2::ggplot(tabla_validacion,
                          ggplot2::aes(x = reorder(Modelo, RMSE_medio),
                                       y = RMSE_medio)) +
  ggplot2::geom_col(fill = "#286B81") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(0, RMSE_medio - sd_RMSE),
                                      ymax = RMSE_medio + sd_RMSE),
                         width = .2) +
  ggplot2::coord_flip() + ggplot2::theme_minimal() +
  ggplot2::labs(title = "RMSE en validación temporal",
                subtitle = "Barras de error: desviación estándar entre 3 folds",
                x = "Modelo", y = "RMSE (MWh)")
salvar_grafico(p_comp, "05_comparacion_cv_rmse.png")

# PASO 13. CALIBRACIÓN Y EVALUACIÓN TEST 2015 ---------------------------------
modelo_para_calibrar <- ajustar_modelo(ganador, cfg_ganador, desarrollo,
                                       semilla = SEMILLA + 100L)
pred_cal <- predecir_modelo(modelo_para_calibrar, calibracion)
residuos_cal <- calibracion$TEY - pred_cal
nivel_cobertura <- .90
n_cal <- length(residuos_cal)
probabilidad <- min(1, ceiling((n_cal + 1) * nivel_cobertura) / n_cal)
q_error <- as.numeric(stats::quantile(abs(residuos_cal), probs = probabilidad,
                                      type = 1, names = FALSE))
calibracion_detalle <- data.frame(
  id_registro = calibracion$id_registro, Year = calibracion$Year,
  observado = calibracion$TEY, predicho = pred_cal,
  residuo = residuos_cal, abs_residuo = abs(residuos_cal))
salvar_csv(calibracion_detalle, "13_residuos_calibracion_2014.csv")

modelo_final <- modelo_para_calibrar
pred_final <- predecir_modelo(modelo_final, test_2015)
media_entrenamiento <- mean(desarrollo$TEY)
pred_baseline_test <- rep(media_entrenamiento, nrow(test_2015))
met_ganador <- metricas(test_2015$TEY, pred_final)
met_baseline <- metricas(test_2015$TEY, pred_baseline_test)

tabla_test <- dplyr::bind_rows(
  data.frame(Modelo = "Baseline", met_baseline),
  data.frame(Modelo = ganador, met_ganador)
)
tabla_test$Reduccion_RMSE_vs_baseline_pct <-
  100 * (met_baseline$RMSE - tabla_test$RMSE) / met_baseline$RMSE
salvar_csv(tabla_test, "14_TEST_FINAL_2015.csv")
df_test_grafico <- data.frame(
  Modelo = factor(rep(c("Baseline (media)", ganador), each = 2), 
                  levels = c("Baseline (media)", ganador)),
  Metrica = factor(rep(c("RMSE", "MAE"), times = 2), 
                   levels = c("RMSE", "MAE")),
  Valor = c(met_baseline$RMSE, met_baseline$MAE, met_ganador$RMSE, met_ganador$MAE)
)

# Cálculo dinámico de métricas para la leyenda inferior
pct_mejora <- 100 * (met_baseline$RMSE - met_ganador$RMSE) / met_baseline$RMSE
r2_val <- met_ganador$R2
n_test <- nrow(test_2015)

texto_pie <- sprintf("Reducción del RMSE frente al baseline: %.2f%%   |   R² de %s: %.4f",
                     pct_mejora, ganador, r2_val)
texto_pie <- gsub("\\.", ",", texto_pie)  # Usar formato con coma decimal

p_test_errores <- ggplot2::ggplot(df_test_grafico, 
                                  ggplot2::aes(x = Modelo, y = Valor, fill = Metrica)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.65) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = gsub("\\.", ",", sprintf("%.4f", Valor)),
      vjust = -0.5
    ),
    position = ggplot2::position_dodge(width = 0.75),
    size = 3.6,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(values = c("RMSE" = "#2C4D63", "MAE" = "#009688")) +
  ggplot2::scale_y_continuous(limits = c(0, max(df_test_grafico$Valor) * 1.15), 
                              breaks = seq(0, 17.5, by = 2.5)) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
    plot.caption = ggplot2::element_text(size = 9, face = "bold", hjust = 0.5, color = "#333333"),
    legend.position = c(0.88, 0.88),
    legend.title = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    title = sprintf("Prueba final 2015 (n = %s) — errores fuera de muestra", 
                    format(n_test, big.mark = ".")),
    x = NULL,
    y = "Error (MWh)",
    caption = texto_pie
  )

salvar_grafico(p_test_errores, "08b_test_2015_comparacion_errores.png", ancho = 9, alto = 5.5)

p75_at <- as.numeric(stats::quantile(desarrollo$AT, .75, names = FALSE))
en_alta_temp <- test_2015$AT > p75_at
if (!any(en_alta_temp)) stop("El subgrupo de alta temperatura quedó vacío.")
met_subgan <- metricas(test_2015$TEY[en_alta_temp], pred_final[en_alta_temp])
met_subbase <- metricas(test_2015$TEY[en_alta_temp], pred_baseline_test[en_alta_temp])
tabla_subgrupo <- data.frame(
  umbral_AT_P75_train_C = p75_at, n_test_subgrupo = sum(en_alta_temp),
  ganador = ganador, RMSE_ganador_MWh = met_subgan$RMSE,
  RMSE_baseline_MWh = met_subbase$RMSE,
  mejora_RMSE_pct = 100 * (met_subbase$RMSE - met_subgan$RMSE) / met_subbase$RMSE)
salvar_csv(tabla_subgrupo, "15_TEST_SUBGRUPO_ALTA_TEMPERATURA.csv")

limite_inf <- pred_final - q_error
limite_sup <- pred_final + q_error
cobertura_test <- mean(test_2015$TEY >= limite_inf & test_2015$TEY <= limite_sup)
tabla_intervalos <- data.frame(
  ganador = ganador, cobertura_nominal = nivel_cobertura,
  n_calibracion = n_cal, cuantil_abs_residuo_MWh = q_error,
  cobertura_empirica_2015 = cobertura_test,
  amplitud_media_MWh = mean(limite_sup - limite_inf),
  nota = "Mismo modelo desarrollo en calibración y test; aproximado por dependencia temporal")
salvar_csv(tabla_intervalos, "16_INTERVALOS_PREDICCION.csv")

cob_nominal <- nivel_cobertura * 100
cob_observada <- cobertura_test * 100
brecha_pp <- cob_observada - cob_nominal

df_cobertura <- data.frame(
  Categoria = factor(c("Nominal esperado", "Observado en 2015"), 
                     levels = c("Nominal esperado", "Observado en 2015")),
  Porcentaje = c(cob_nominal, cob_observada),
  Color = c("Nominal", "Observado")
)

texto_brecha <- sprintf("Brecha: %.2f puntos porcentuales. Los intervalos subcubren en la prueba final.", brecha_pp)
texto_brecha <- gsub("\\.", ",", texto_brecha)

p_cobertura <- ggplot2::ggplot(df_cobertura, ggplot2::aes(x = Categoria, y = Porcentaje, fill = Color)) +
  ggplot2::geom_col(width = 0.55) +
  ggplot2::geom_text(
    ggplot2::aes(label = gsub("\\.", ",", sprintf("%.2f %%", Porcentaje))),
    vjust = -0.6,
    size = 4,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(values = c("Nominal" = "#2C4D63", "Observado" = "#009688")) +
  ggplot2::scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, by = 20)) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "none",
    plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
    plot.caption = ggplot2::element_text(size = 8.5, hjust = 0, color = "#555555"),
    panel.grid.major.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(face = "bold", size = 10)
  ) +
  ggplot2::labs(
    title = "Intervalos de predicción al 90% — calibración insuficiente",
    x = NULL,
    y = "Porcentaje de cobertura (%)",
    caption = texto_brecha
  )

salvar_grafico(p_cobertura, "09a_intervalos_prediccion_cobertura.png", ancho = 7, alto = 5)

rmse_total_val <- met_ganador$RMSE
rmse_subgrupo_val <- met_subgan$RMSE

df_subgrupo_at <- data.frame(
  Grupo = factor(c("Total 2015", sprintf("Subgrupo cálido\nAT > %.2f °C", p75_at)), 
                 levels = c("Total 2015", sprintf("Subgrupo cálido\nAT > %.2f °C", p75_at))),
  RMSE = c(rmse_total_val, rmse_subgrupo_val),
  Color = c("Total", "Subgrupo")
)
# Reemplazar punto por coma para las etiquetas del gráfico
df_subgrupo_at$Grupo <- gsub("\\.", ",", df_subgrupo_at$Grupo)

p_subgrupo <- ggplot2::ggplot(df_subgrupo_at, ggplot2::aes(x = Grupo, y = RMSE, fill = Color)) +
  ggplot2::geom_col(width = 0.55) +
  ggplot2::geom_text(
    ggplot2::aes(label = gsub("\\.", ",", sprintf("%.4f", RMSE))),
    vjust = -0.6,
    size = 4,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(values = c("Total" = "#2C4D63", "Subgrupo" = "#009688")) +
  ggplot2::scale_y_continuous(limits = c(0, max(df_subgrupo_at$RMSE) * 1.25), breaks = seq(0, 2.0, by = 0.25)) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "none",
    plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
    plot.caption = ggplot2::element_text(size = 8, hjust = 0.5, color = "#555555"),
    panel.grid.major.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(face = "bold", size = 9.5)
  ) +
  ggplot2::labs(
    title = "Error en condiciones de mayor temperatura ambiental",
    x = NULL,
    y = "RMSE (MWh)",
    caption = "El tamaño n del subgrupo no figura en el informe; esta comparación no demuestra mayor robustez por sí sola."
  )

salvar_grafico(p_subgrupo, "09b_subgrupo_alta_temperatura_rmse.png", ancho = 7, alto = 5)


predicciones_test <- data.frame(
  Year = test_2015$Year, orden_en_anio = test_2015$orden_en_anio,
  TEY_real_MWh = test_2015$TEY, TEY_predicha_MWh = pred_final,
  pred_baseline_MWh = pred_baseline_test,
  residuo_MWh = test_2015$TEY - pred_final,
  AT_C = test_2015$AT, alta_temperatura = en_alta_temp,
  limite_inferior_90 = limite_inf, limite_superior_90 = limite_sup)
salvar_csv(predicciones_test, "17_PREDICCIONES_TEST_2015.csv")

residuos_test <- predicciones_test$residuo_MWh
p_resfit <- ggplot2::ggplot(predicciones_test,
                            ggplot2::aes(x = TEY_predicha_MWh, y = residuo_MWh)) +
  ggplot2::geom_point(alpha = .25, color = "#355C7D") +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "Residuos vs valores predichos - test 2015",
                x = "TEY predicha (MWh)", y = "Residuo: real - predicha (MWh)")
salvar_grafico(p_resfit, "06_residuos_vs_predichos_test.png")

p_hist_res <- ggplot2::ggplot(predicciones_test,
                              ggplot2::aes(x = residuo_MWh)) +
  ggplot2::geom_histogram(bins = 45, fill = "#355C7D", color = "white") +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "Distribución de residuos - test 2015",
                x = "Residuo (MWh)", y = "Frecuencia")
salvar_grafico(p_hist_res, "07_distribucion_residuos_test.png")

p_test <- ggplot2::ggplot(predicciones_test,
                          ggplot2::aes(x = TEY_real_MWh, y = TEY_predicha_MWh)) +
  ggplot2::geom_point(alpha = .25, color = "#355C7D") +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = paste("Test 2015: observado vs predicho -", ganador),
                x = "TEY real (MWh)", y = "TEY predicha (MWh)")
salvar_grafico(p_test, "08_real_vs_predicho_test.png")

p_cv_res <- ggplot2::ggplot(dplyr::filter(oof, Modelo == ganador),
                            ggplot2::aes(x = Fold, y = residuo)) +
  ggplot2::geom_boxplot(fill = "#85C1CC") + ggplot2::theme_minimal() +
  ggplot2::labs(title = "Residuos fuera de muestra por bloque temporal",
                x = "Fold de validación", y = "Residuo (MWh)")
salvar_grafico(p_cv_res, "09_residuos_por_periodo_validado.png")

png(file.path(carpeta_salida, "10_acf_residuos_test.png"),
    width = 1400, height = 900, res = 150)
stats::acf(residuos_test, lag.max = 72, na.action = na.pass,
           main = "ACF de residuos, test 2015 (orden cronológico)")
grDevices::dev.off()

cuantiles_pred <- unique(as.numeric(stats::quantile(pred_final,
                                                    probs = seq(0, 1, .25), names = FALSE)))
if (length(cuantiles_pred) >= 3L) {
  grupo <- cut(pred_final, breaks = cuantiles_pred, include.lowest = TRUE)
  heteroc <- dplyr::as_tibble(predicciones_test) |>
    dplyr::mutate(tramo_predicho = grupo) |>
    dplyr::group_by(tramo_predicho) |>
    dplyr::summarise(n = dplyr::n(), RMSE = sqrt(mean(residuo_MWh^2)),
                     sd_residuo = stats::sd(residuo_MWh), .groups = "drop")
  salvar_csv(heteroc, "18_diagnostico_error_por_nivel_predicho.csv")
}

# PASO 14. IMPORTANCIA Y PDP/ICE -----------------------------------------------
fold_imp <- folds[[3]]
modelo_imp <- ajustar_modelo(ganador, cfg_ganador, fold_imp$train,
                             semilla = SEMILLA + fold_imp$n)
set.seed(SEMILLA + 303L)
idx_imp <- sample.int(nrow(fold_imp$valid), min(2000L, nrow(fold_imp$valid)))
valid_imp <- fold_imp$valid[idx_imp, , drop = FALSE]
rmse_ref <- metricas(valid_imp$TEY,
                     predecir_modelo(modelo_imp, valid_imp))$RMSE
imp <- lapply(seq_along(predictores), function(j) {
  variable <- predictores[j]
  incrementos <- vapply(seq_len(3L), function(replica) {
    tmp <- valid_imp
    set.seed(SEMILLA + 100L * j + replica)
    tmp[[variable]] <- sample(tmp[[variable]])
    metricas(valid_imp$TEY,
             predecir_modelo(modelo_imp, tmp))$RMSE - rmse_ref
  }, numeric(1))
  data.frame(Variable = variable, Incremento_RMSE_MWh = mean(incrementos),
             sd_3_permutaciones = stats::sd(incrementos))
})
importancia <- dplyr::bind_rows(imp) |>
  dplyr::arrange(dplyr::desc(Incremento_RMSE_MWh))
salvar_csv(importancia, "19_importancia_permutacion_validacion.csv")

p_imp <- ggplot2::ggplot(importancia,
                         ggplot2::aes(x = reorder(Variable, Incremento_RMSE_MWh),
                                      y = Incremento_RMSE_MWh)) +
  ggplot2::geom_col(fill = "#286B81") + ggplot2::coord_flip() +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "Importancia predictiva por permutación (2014 inicial)",
                subtitle = "Aumento del RMSE al mezclar una variable; no es causalidad",
                x = "Variable", y = "Incremento de RMSE (MWh)")
salvar_grafico(p_imp, "11_importancia_permutacion.png")

# CORRECCIÓN EN REJILLA: se usa 'extremos[2]' en lugar de 'extremes[2]'
set.seed(SEMILLA + 400L)
idx_pdp <- sample.int(nrow(desarrollo), min(200L, nrow(desarrollo)))
base_pdp <- desarrollo[idx_pdp, , drop = FALSE]
variables_pdp <- c("CDP", "GTEP", "AT")
pdp_registros <- list()
ice_registros <- list()
for (v in variables_pdp) {
  extremos <- stats::quantile(desarrollo[[v]], c(.05, .95), names = FALSE)
  rejilla <- seq(extremos[1], extremos[2], length.out = 15L)
  for (valor in rejilla) {
    datos_variados <- base_pdp
    datos_variados[[v]] <- valor
    predicciones <- predecir_modelo(modelo_final, datos_variados)
    pdp_registros[[length(pdp_registros) + 1L]] <- data.frame(
      Variable = v, Valor = valor, PDP_MWh = mean(predicciones))
    ice_registros[[length(ice_registros) + 1L]] <- data.frame(
      Variable = v, Valor = valor, id_ICE = seq_len(min(25L, nrow(base_pdp))),
      ICE_MWh = head(predicciones, 25L))
  }
}
pdp_df <- dplyr::bind_rows(pdp_registros)
ice_df <- dplyr::bind_rows(ice_registros)
salvar_csv(pdp_df, "20_pdp.csv")
salvar_csv(ice_df, "21_ice.csv")

p_pdp <- ggplot2::ggplot() +
  ggplot2::geom_line(data = ice_df,
                     ggplot2::aes(x = Valor, y = ICE_MWh, group = id_ICE),
                     linewidth = .25, alpha = .25, color = "#999999") +
  ggplot2::geom_line(data = pdp_df,
                     ggplot2::aes(x = Valor, y = PDP_MWh), linewidth = 1.1, color = "#1B536E") +
  ggplot2::facet_wrap(~ Variable, scales = "free_x") +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "PDP (azul) e ICE (gris): variables operacionales",
                subtitle = "Efectos del modelo, no relaciones causales",
                x = "Valor de la variable", y = "TEY estimada (MWh)")
salvar_grafico(p_pdp, "12_pdp_ice.png", ancho = 12, alto = 5)

info_ejecucion <- data.frame(
  semilla = SEMILLA, sistema_ganador = ganador,
  version_script = "final_8_predictoras_2026_11_05", arboles_random_forest = ARBOLES_RF,
  profundidad_max_random_forest = PROFUNDIDAD_MAX_RF,
  hilos_maximos_utilizados = HILOS,
  anios_desarrollo = "2011-2013 y primer 75% 2014",
  periodo_calibracion = "ultimo 25% 2014",
  anio_test = 2015L, folds = length(folds),
  n_desarrollo = nrow(desarrollo), n_calibracion = nrow(calibracion),
  n_train_total = nrow(train_total), n_train_modelo_evaluado = nrow(desarrollo),
  n_test = nrow(test_2015),
  umbral_alta_temp_AT_C = p75_at,
  calibracion_nivel = nivel_cobertura,
  umbral_correlacion_XX = UMBRAL_COR,
  regla_depuracion = "NO se eliminan variables; las ocho X se conservan en todos los modelos",
  criterio_intervalos = "modelo desarrollo sin reajuste sobre calibración; sin garantía con dependencia temporal",
  modelo_averaging_pesos = "RF=0.50;XGB=0.30;GAM=0.20",
  tipo_prediccion = "Estimación contemporánea, no pronóstico de una hora",
  unidades = "TEY MWh por período horario")
salvar_csv(info_ejecucion, "22_CONFIGURACION_EJECUCION.csv")

capture.output(utils::sessionInfo(),
               file = file.path(carpeta_salida, "23_versiones_R_y_paquetes.txt"))
capture.output(str(cfg_ganador),
               file = file.path(carpeta_salida, "24_hiperparametros_ganador.txt"))
saveRDS(modelo_final, file.path(carpeta_salida, "25_modelo_evaluado_calibrado.rds"))

md5_final <- unname(tools::md5sum(rutas_originales))
stopifnot(identical(md5_original, md5_final))

cat("\n===============================================================\n")
cat("FIN: script FINAL de regresión ejecutado exitosamente\n")
cat("Ganador según CV: ", ganador, "\n", sep = "")
cat("RMSE test 2015: ", round(met_ganador$RMSE, 4), " MWh\n", sep = "")
cat("Salidas guardadas en: ", normalizePath(carpeta_salida), "\n", sep = "")
cat("===============================================================\n")
