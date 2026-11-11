# Caso de Estudio N° 1: Regresión — Predicción de la Generación Eléctrica de una Turbina de Gas Mediante Regresión Lineal

**Universidad de Santiago de Chile (USACH) — Estadística Aplicada**

## 📂 Estructura del Repositorio
Para garantizar la ejecución limpia y reproducible del pipeline, el repositorio se organiza de la siguiente manera:
- 📁 `dataset/` : Contiene los 5 archivos CSV oficiales de los años 2011 al 2015 (`gt_2011.csv`, `gt_2012.csv`, `gt_2013.csv`, `gt_2014.csv`, `gt_2015.csv`).
- 📄 `código.R` : Script principal en R que automatiza todo el análisis exploratorio, modelado, validación temporal y exportación de resultados.
- 📁 `resultados_finales/` : Carpeta generada automáticamente por el script que almacena las tablas CSV y las 12 figuras oficiales utilizadas en el informe y la presentación.

---

## ⚙️ Modo de Uso y Ejecución del Código

Para reproducir por completo el análisis, las métricas, las tablas de validación y las figuras de la presentación, siga estos pasos:

1. **Requisito de Ubicación:** 
   Asegúrese de que el archivo **`codigo.R`** y la carpeta **`dataset/`** (con los archivos CSV adentro) se encuentren ubicados **exactamente en el mismo directorio de trabajo (misma carpeta)** en su computador.

2. **Entorno de R / RStudio:**
   Abra RStudio y asegúrese de tener instaladas las librerías necesarias para el modelado y la manipulación de datos (ej. `tidyverse`, `caret`, `glmnet`, `xgboost`, `mgcv`, `car`).

3. **Ejecución del Script:**
   Establezca su directorio de trabajo en la carpeta del proyecto y ejecute el script completo mediante el comando:
   ```R
   source("script_regresion_OPTIMIZADO.R")
