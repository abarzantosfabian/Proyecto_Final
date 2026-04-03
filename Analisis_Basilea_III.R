# =============================================================================
# EVALUACIÓN DE IMPACTO BASILEA III — SISTEMA BANCARIO CHILENO
# =============================================================================
#
# PROPÓSITO
#   Estimar el efecto causal de la implementación gradual de Basilea III sobre
#   el stock de colocaciones de los bancos chilenos mediante Diferencias-en-
#   Diferencias (DiD) con panel bancario mensual 2009-2025.
#
# ESTRATEGIA DE IDENTIFICACIÓN
#   Dado que Basilea III afectó a toda la industria, el grupo de tratamiento se
#   define por exposición relativa pre-reforma. Tres criterios de clasificación:
#     1. Holgura de capital (CET1/APR − umbral regulatorio) < mediana (P50)
#     2. Holgura de capital < percentil 25 (robustez)
#     3. K-Means en 2 etapas sobre variables pre-reforma (robustez alternativa)
#
# PIPELINE
#   [1] Datos sintéticos representativos
#   [2] Construcción del panel
#   [3] Umbrales regulatorios y holgura de capital
#   [4] Clasificación PCA + K-Means en 2 etapas
#   [5] Bootstrap de estabilidad del clustering
#   [6] Variables de tratamiento y dummies de hitos
#   [7] Modelos DiD (acumulado y escalonado)
#   [8] Bootstrap del parámetro DiD
#   [9] Tablas de resultados
#  [10] Figuras
#
# REPRODUCIBILIDAD
#   Semillas fijadas con set.seed() en cada paso estocástico (SEED_KMEANS=123).
#
# DEPENDENCIAS
#   dplyr, tidyr, ggplot2, lubridate, plm, lmtest, sandwich,
#   knitr, kableExtra, scales, patchwork, purrr, stringr, ggrepel
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(plm)
  library(lmtest)
  library(sandwich)
  library(knitr)
  library(kableExtra)
  library(scales)
  library(patchwork)
  library(purrr)
  library(stringr)
  library(ggrepel)
})

options(OutDec = ",")
SEED_KMEANS <- 123


# =============================================================================
# [1] DATOS: MUESTRA SINTÉTICA REPRESENTATIVA
# =============================================================================
# Se construye un panel sintético que replica la estructura de la data original:
# 12 bancos × 201 meses (enero 2009 – septiembre 2025). Las distribuciones
# imitan los momentos reales del sistema bancario chileno (escala MM CLP,
# CET1/APR entre 8% y 20%, colocaciones correlacionadas con activos totales).
# =============================================================================

set.seed(2024)

BANCOS_INFO <- tibble::tribble(
  ~INS_COD, ~Bancos,          ~at_base,  ~cet1_mean, ~cet1_sd,
  1L,       "Chile",           2.5e7,     10.5,        1.2,
  9L,       "Internacional",   2.0e5,      9.0,        2.0,
  12L,      "Estado",          3.0e7,      8.0,        1.5,
  14L,      "Scotiabank",      3.5e6,      9.5,        1.8,
  16L,      "BCI",             1.8e7,      9.8,        1.0,
  28L,      "Bice",            2.5e6,      9.2,        1.3,
  37L,      "Santander",       2.2e7,     10.2,        1.1,
  39L,      "Itau",            5.0e6,     10.8,        2.2,
  49L,      "Security",        1.5e6,      9.0,        1.5,
  51L,      "Falabella",       4.5e5,     15.0,        5.0,
  53L,      "Ripley",          3.0e5,     18.0,        4.0,
  55L,      "Consorcio",       6.0e5,     12.0,        8.0
)

PERIODOS    <- seq.Date(as.Date("2009-01-01"), as.Date("2025-09-01"), by = "month")
N_T         <- length(PERIODOS)
macro_trend <- cumsum(c(0, rnorm(N_T - 1, mean = 0.003, sd = 0.008)))
covid_shock <- ifelse(PERIODOS >= as.Date("2020-03-01") &
                        PERIODOS <= as.Date("2021-06-01"),
                      rnorm(N_T, -0.04, 0.02), 0)

BANCA <- map_dfr(seq_len(nrow(BANCOS_INFO)), function(i) {
  b  <- BANCOS_INFO[i, ]
  at <- b$at_base * exp(macro_trend + rnorm(N_T, 0, 0.015))
  tibble(
    INS_COD         = b$INS_COD,
    Bancos          = b$Bancos,
    PERIODO         = PERIODOS,
    AT              = at,
    APR             = at * runif(N_T, 0.55, 0.75),
    CET1_APR        = pmax(b$cet1_mean + cumsum(rnorm(N_T, 0.01, 0.08)), 5),
    Coloc_Total     = at * runif(N_T, 0.50, 0.72),
    Coloc_Consumo   = at * runif(N_T, 0.08, 0.18),
    Coloc_Comercial = at * runif(N_T, 0.25, 0.45),
    Coloc_Vivienda  = at * runif(N_T, 0.10, 0.22),
    depositos       = at * runif(N_T, 0.45, 0.65),
    efectivo        = at * runif(N_T, 0.03, 0.08),
    utilidad        = at * runif(N_T, 0.005, 0.015),
    patrimonio      = at * runif(N_T, 0.07, 0.12)
  )
})

IMACEC    <- tibble(PERIODO  = PERIODOS,
                    imacec   = 100 * exp(macro_trend + covid_shock))
DESEMPLEO <- tibble(PERIODO  = PERIODOS,
                    desempleo = pmax(5, 8 + cumsum(rnorm(N_T, 0, 0.05)) +
                                      ifelse(PERIODOS >= as.Date("2020-03-01") &
                                               PERIODOS <= as.Date("2021-06-01"), 3, 0)))

cat("✓ [1] Datos sintéticos:", comma(nrow(BANCA)),
    "observaciones banco-mes generadas\n")


# =============================================================================
# [2] CONSTRUCCIÓN DEL PANEL
# =============================================================================
# Se une BANCA con las variables macro y se filtra a los 12 bancos de la
# muestra de regresión. La unidad de análisis es banco-mes en base consolidada.
# =============================================================================

bancos_regresion <- c(1, 9, 12, 14, 16, 28, 37, 39, 49, 55, 51, 53)

BANCA_PANEL <- BANCA %>%
  left_join(IMACEC,    by = "PERIODO") %>%
  left_join(DESEMPLEO, by = "PERIODO") %>%
  filter(INS_COD %in% bancos_regresion) %>%
  mutate(PERIODO = as.Date(PERIODO), ANIO = year(PERIODO))

cat("✓ [2] Panel:", n_distinct(BANCA_PANEL$INS_COD), "bancos ×",
    n_distinct(BANCA_PANEL$PERIODO), "períodos =",
    comma(nrow(BANCA_PANEL)), "obs.\n")


# =============================================================================
# [3] UMBRALES REGULATORIOS Y HOLGURA DE CAPITAL
# =============================================================================
# Cada banco tiene un umbral CET1/APR diferenciado por año bajo la transición
# de Basilea III (2021-2024). La holgura se define como:
#   Holgura_it = max(CET1/APR_it - Umbral_it, 0)
# Un banco con holgura baja pre-reforma es más vulnerable a los nuevos
# requerimientos → mayor probabilidad de pertenecer al grupo tratado.
# =============================================================================

UMBRALES <- data.frame(
  ANIO = 2009:2025,
  `1`  = c(10,10,10,10,10,10,10,10,10,10,10,10,10, 9.25, 9.25, 9.75, 9.38),
  `9`  = c( 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8.00, 8.00, 9.25, 9.25),
  `12` = c( 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9.25, 9.25,10.00, 9.50),
  `14` = c( 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9.25, 9.25,10.25, 9.50),
  `16` = c( 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9.50, 9.75, 9.75, 9.50),
  `28` = c( 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8.00, 8.00, 8.50, 8.90),
  `37` = c(11,11,11,11,11,11,11,11,11,11,11,11,11, 9.50, 9.50, 9.50, 9.63),
  `39` = c( 8, 8, 8, 8, 8, 8,10,10,10,10,10,10,10, 9.00, 9.00, 9.00, 9.00),
  `49` = c( 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8.00, 8.00, 9.25, 9.25),
  `51` = c( 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8.00, 8.00, 8.00, 8.00),
  `53` = c( 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8.00, 8.00, 8.00, 8.00),
  `55` = c( 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8.00, 8.00, 8.75,10.50),
  check.names = FALSE
)

UMBRALES_long <- UMBRALES %>%
  pivot_longer(-ANIO, names_to = "INS_COD", values_to = "UMBRAL") %>%
  mutate(INS_COD = as.integer(INS_COD))

BANCA_PANEL <- BANCA_PANEL %>%
  left_join(UMBRALES_long, by = c("INS_COD", "ANIO")) %>%
  mutate(HOLGURA = pmax(CET1_APR - UMBRAL, 0))

cat("✓ [3] Holgura calculada. Promedio industria pre-oct.2018:",
    round(mean(BANCA_PANEL$HOLGURA[BANCA_PANEL$PERIODO < as.Date("2018-10-01")],
               na.rm = TRUE), 2), "pp\n")


# =============================================================================
# [4] CLASIFICACIÓN K-MEANS EN 2 ETAPAS
# =============================================================================
# Etapa 1 — Tipología estructural:
#   K-Means (k=2) sobre variables de tamaño y composición de cartera
#   (log_AT, Con_AT, Com_AT, Viv_AT, RWA_density) estandarizadas.
#   Identifica dos tipologías de negocio: banca retail vs. banca corporativa.
#
# Etapa 2 — Exposición relativa dentro de cada tipología:
#   Dentro de cada cluster estructural, un banco es "más expuesto" si su
#   holgura promedio pre-reforma es >= mediana de su tipología.
#   Esto evita comparar bancos estructuralmente diferentes solo por holgura.
#
# Variables de clustering (calculadas en período pre-reforma):
#   log_AT      : log(Activos Totales) — proxy de tamaño
#   Con_AT      : Colocaciones consumo / AT — mix de cartera
#   Com_AT      : Colocaciones comerciales / AT
#   Viv_AT      : Colocaciones vivienda / AT
#   RWA_density : APR / AT — densidad de activos ponderados por riesgo
#   Exposure    : Holgura de capital promedio pre-reforma
# =============================================================================

hito_0      <- as.Date("2018-10-01")
hito_2      <- as.Date("2020-12-01")
cutoff_date <- as.Date("2021-12-01")
hito_4      <- as.Date("2022-12-01")
hito_5      <- as.Date("2023-12-01")
hito_6      <- as.Date("2024-12-01")

feature_stage1 <- c("log_AT", "Con_AT", "Com_AT", "Viv_AT", "RWA_density")

BANCA_CLUSTER_BASE <- BANCA_PANEL %>%
  filter(PERIODO <= hito_0) %>%
  mutate(
    INS_COD     = as.character(INS_COD),
    Bancos      = as.character(Bancos),
    log_AT      = log(as.numeric(AT) + 1),
    Con_AT      = if_else(is.finite(Coloc_Consumo / AT),   Coloc_Consumo / AT,   NA_real_),
    Com_AT      = if_else(is.finite(Coloc_Comercial / AT), Coloc_Comercial / AT, NA_real_),
    Viv_AT      = if_else(is.finite(Coloc_Vivienda / AT),  Coloc_Vivienda / AT,  NA_real_),
    RWA_density = if_else(is.finite(APR / AT), APR / AT, NA_real_),
    Exposure    = HOLGURA
  )

run_kmeans_two_stage <- function(df_pre,
                                  structural_vars = feature_stage1,
                                  exposure_var    = "Exposure",
                                  nstart          = 200,
                                  seed            = 123,
                                  centers         = 2) {
  bank_profiles <- df_pre %>%
    group_by(INS_COD, Bancos) %>%
    summarise(across(all_of(c(structural_vars, exposure_var)),
                     ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    filter(if_all(all_of(c(structural_vars, exposure_var)), is.finite))

  if (nrow(bank_profiles) < centers)
    stop("No hay suficientes bancos para clustering en 2 etapas.")

  X_struct <- scale(bank_profiles[, structural_vars])
  set.seed(seed)
  km1 <- kmeans(X_struct, centers = centers, nstart = nstart)

  out <- bank_profiles %>% mutate(cluster_estructural = km1$cluster)

  thresholds <- out %>%
    group_by(cluster_estructural) %>%
    summarise(mediana_exposure = median(.data[[exposure_var]], na.rm = TRUE),
              n_bancos = n(), .groups = "drop")

  out <- out %>%
    left_join(thresholds, by = "cluster_estructural") %>%
    mutate(
      expuesto_2stage = if_else(.data[[exposure_var]] >= mediana_exposure, 1L, 0L),
      grupo_2stage    = if_else(expuesto_2stage == 1L, "Más expuesto", "Menos expuesto")
    ) %>%
    arrange(cluster_estructural, desc(.data[[exposure_var]]), INS_COD) %>%
    group_by(cluster_estructural) %>%
    mutate(ranking_tipologia = row_number()) %>%
    ungroup()

  attr(out, "km_stage1")         <- km1
  attr(out, "stage2_thresholds") <- thresholds
  out
}

kmeans_2stage <- run_kmeans_two_stage(
  df_pre = BANCA_CLUSTER_BASE, nstart = 500, seed = SEED_KMEANS
)

cat("✓ [4] K-Means en 2 etapas:\n")
cat("  Cluster 1:", sum(kmeans_2stage$cluster_estructural == 1), "bancos\n")
cat("  Cluster 2:", sum(kmeans_2stage$cluster_estructural == 2), "bancos\n")
cat("  Más expuestos:",
    paste(kmeans_2stage$Bancos[kmeans_2stage$expuesto_2stage == 1], collapse = ", "), "\n")


# =============================================================================
# [5] BOOTSTRAP DE ESTABILIDAD DEL CLUSTERING
# =============================================================================
# Bootstrap temporal: se remuestrean los períodos pre-reforma con reemplazamiento
# y se re-estima la clasificación K-Means en 2 etapas en cada réplica.
#
# Indicador clave: P(banco i clasificado como "más expuesto")
#   ≈ 1.0 → asignación muy estable al grupo tratado
#   ≈ 0.0 → asignación muy estable al grupo control
#   ≈ 0.5 → banco en la frontera; clasificación sensible al período muestral
#
# Estabilidad: Alta si P ≥ 0.90 ó P ≤ 0.10 | Media si ≥0.75 ó ≤0.25 | Baja si entre.
# =============================================================================

B_BOOT <- 499

boot_kmeans_membership <- function(df_pre, B = 499,
                                    structural_vars = feature_stage1,
                                    seed = 123) {
  fechas_pre <- sort(unique(df_pre$PERIODO))

  map_dfr(seq_len(B), function(b) {
    set.seed(seed + b)
    boot_df <- map_dfr(
      sample(fechas_pre, size = length(fechas_pre), replace = TRUE),
      ~ df_pre %>% filter(PERIODO == .x)
    )

    km_b <- tryCatch(
      run_kmeans_two_stage(boot_df, structural_vars, nstart = 100, seed = seed + b),
      error = function(e) NULL
    )

    if (is.null(km_b)) {
      return(df_pre %>% distinct(INS_COD, Bancos) %>%
               mutate(boot_id = b, expuesto_boot = NA_integer_))
    }

    km_b %>% transmute(boot_id = b, INS_COD, Bancos,
                        expuesto_boot = expuesto_2stage)
  })
}

boot_membership_long <- boot_kmeans_membership(BANCA_CLUSTER_BASE, B_BOOT, seed = SEED_KMEANS)

bank_labels_base <- kmeans_2stage %>%
  transmute(INS_COD = as.character(INS_COD), Bancos,
            expuesto_base = as.integer(expuesto_2stage))

stability_bank <- boot_membership_long %>%
  group_by(INS_COD, Bancos) %>%
  summarise(prob_expuesto_boot = mean(expuesto_boot, na.rm = TRUE), .groups = "drop") %>%
  mutate(INS_COD = as.character(INS_COD)) %>%
  left_join(bank_labels_base, by = c("INS_COD", "Bancos")) %>%
  mutate(
    estabilidad = case_when(
      prob_expuesto_boot >= 0.90 | prob_expuesto_boot <= 0.10 ~ "Alta",
      prob_expuesto_boot >= 0.75 | prob_expuesto_boot <= 0.25 ~ "Media",
      TRUE ~ "Baja"
    )
  ) %>%
  arrange(desc(prob_expuesto_boot))

cat("✓ [5] Bootstrap de estabilidad (", B_BOOT, "réplicas)\n")
cat("  Alta estabilidad:", sum(stability_bank$estabilidad == "Alta"), "bancos\n")


# =============================================================================
# [6] VARIABLES DE TRATAMIENTO Y DUMMIES DE HITOS
# =============================================================================
# Variables de tratamiento (asignación fija por banco, definida pre-reforma):
#   dummy_holg     : 1 si holgura promedio pre-oct.2018 < P50 industria
#   dummy_holg25   : 1 si holgura promedio pre-oct.2018 < P25 industria (robustez)
#   dummy_holg_km  : 1 si banco clasificado como "más expuesto" por K-Means en 2 etapas
#
# Dummies de hitos (varían en el tiempo):
#   post_h0 : 1 si PERIODO >= oct.2018 (Ley 21.130)
#   post_h2 : 1 si PERIODO >= dic.2020 (inicio implementación — fecha base DiD)
#   post    : 1 si PERIODO >= dic.2021
#   post_h4–h6: dic.2022, dic.2023, dic.2024
# =============================================================================

avrg_pre <- BANCA_PANEL %>%
  filter(PERIODO < hito_0) %>%
  group_by(INS_COD) %>%
  summarise(avrg_holg = mean(HOLGURA, na.rm = TRUE), .groups = "drop")

q_p50 <- quantile(BANCA_PANEL$HOLGURA[BANCA_PANEL$PERIODO < hito_0], 0.50, na.rm = TRUE)
q_p25 <- quantile(BANCA_PANEL$HOLGURA[BANCA_PANEL$PERIODO < hito_0], 0.25, na.rm = TRUE)

bancos_mas_expuestos_km <- kmeans_2stage %>%
  filter(expuesto_2stage == 1L) %>% pull(INS_COD) %>% as.character()

treat_bank <- avrg_pre %>%
  mutate(
    INS_COD        = as.character(INS_COD),
    dummy_holg     = as.integer(avrg_holg <= q_p50),
    dummy_holg25   = as.integer(avrg_holg <= q_p25),
    dummy_holg_km  = as.integer(INS_COD %in% bancos_mas_expuestos_km)
  )

# Formateo numérico para tablas
fmt_num <- function(x, d = 2) formatC(x, format = "f", digits = d,
                                       decimal.mark = ",", big.mark = ".")
fmt_int <- function(x)        formatC(x, format = "d", big.mark = ".")
.stars  <- function(p) ifelse(is.na(p), "",
                               ifelse(p < 0.01, "***", ifelse(p < 0.05, "**",
                                                               ifelse(p < 0.10, "*", ""))))

panel_df <- BANCA_PANEL %>%
  mutate(
    PERIODO   = as.Date(PERIODO),
    INS_COD   = as.character(INS_COD),
    post_h0   = as.integer(PERIODO >= hito_0),
    post_h2   = as.integer(PERIODO >= hito_2),
    post      = as.integer(PERIODO >= cutoff_date),
    post_h4   = as.integer(PERIODO >= hito_4),
    post_h5   = as.integer(PERIODO >= hito_5),
    post_h6   = as.integer(PERIODO >= hito_6),
    across(c(Coloc_Total, AT, depositos, efectivo, utilidad, patrimonio), as.numeric),
    ln_coloc1  = log(Coloc_Total + 1),
    ln_at      = log(AT + 1),
    ln_imacec  = log(imacec),
    cet1_apr   = CET1_APR,
    holgura    = HOLGURA,
    deposit_at = depositos / AT,
    cash_at    = efectivo  / AT,
    ROE        = utilidad  / patrimonio
  ) %>%
  left_join(treat_bank, by = "INS_COD") %>%
  arrange(INS_COD, PERIODO)

cat("✓ [6] Panel final:", comma(nrow(panel_df)), "obs. listos para regresión\n")
cat("  Bancos tratados (Holgura P50):", sum(treat_bank$dummy_holg), "\n")
cat("  Bancos tratados (K-Means 2 etapas):", sum(treat_bank$dummy_holg_km), "\n")


# =============================================================================
# [7] MODELOS DIFF-IN-DIFF
# =============================================================================
# Especificación 1 — Efecto acumulado:
#   ln(Coloc+1)_it = αi + θt + γ·Macro_t + δ·X_it + β·(Post_t × Treat_i) + ε_it
#
# Especificación 2 — Efecto escalonado (hitos regulatorios Dτ):
#   ln(Coloc+1)_it = αi + θt + γ·Macro_t + δ·X_it + Σ βτ·(Dτ_t × Treat_i) + ε_it
#
# Efectos fijos de doble vía (twoways): αi banco + θt tiempo
# Errores estándar: robustos HC1, clusterizados por banco
#
# Controles: ln_at, deposit_at, cash_at, ROE, holgura, desempleo, ln_imacec
# =============================================================================

pdata <- pdata.frame(as.data.frame(panel_df), index = c("INS_COD", "PERIODO"))

.extract_robust <- function(model, term) {
  V  <- vcovHC(model, type = "HC1", cluster = "group")
  ct <- coeftest(model, vcov = V)
  if (!(term %in% rownames(ct))) return(NULL)
  list(est = ct[term, 1], se = ct[term, 2],
       t   = ct[term, 3], p  = ct[term, 4])
}

controles <- "ln_at + deposit_at + cash_at + ROE + holgura + desempleo + ln_imacec"

# Holgura P50 — Especificación principal
M1B_TBFE <- plm(as.formula(paste(
  "ln_coloc1 ~ dummy_holg + post_h2 + dummy_holg:post_h2 +", controles)),
  data = pdata, model = "within", effect = "twoways")

M2B_TBFE <- plm(as.formula(paste(
  "ln_coloc1 ~ dummy_holg +
    post_h0 + dummy_holg:post_h0 +
    post_h2 + dummy_holg:post_h2 +
    post    + dummy_holg:post    +
    post_h4 + dummy_holg:post_h4 +
    post_h5 + dummy_holg:post_h5 +
    post_h6 + dummy_holg:post_h6 +", controles)),
  data = pdata, model = "within", effect = "twoways")

# Holgura P25 — Robustez
M1B_TBFE25 <- plm(as.formula(paste(
  "ln_coloc1 ~ dummy_holg25 + post_h2 + dummy_holg25:post_h2 +", controles)),
  data = pdata, model = "within", effect = "twoways")

M2B_TBFE25 <- plm(as.formula(paste(
  "ln_coloc1 ~ dummy_holg25 +
    post_h0 + dummy_holg25:post_h0 +
    post_h2 + dummy_holg25:post_h2 +
    post    + dummy_holg25:post    +
    post_h4 + dummy_holg25:post_h4 +
    post_h5 + dummy_holg25:post_h5 +
    post_h6 + dummy_holg25:post_h6 +", controles)),
  data = pdata, model = "within", effect = "twoways")

# K-Means 2 etapas — Robustez alternativa
M1B_TBFE_KM <- plm(as.formula(paste(
  "ln_coloc1 ~ dummy_holg_km + post_h2 + dummy_holg_km:post_h2 +",
  gsub("holgura", "holgura", controles))),
  data = pdata, model = "within", effect = "twoways")

M2B_TBFE_KM <- plm(as.formula(paste(
  "ln_coloc1 ~ dummy_holg_km +
    post_h0 + dummy_holg_km:post_h0 +
    post_h2 + dummy_holg_km:post_h2 +
    post    + dummy_holg_km:post    +
    post_h4 + dummy_holg_km:post_h4 +
    post_h5 + dummy_holg_km:post_h5 +
    post_h6 + dummy_holg_km:post_h6 +", controles)),
  data = pdata, model = "within", effect = "twoways")

cat("✓ [7] Modelos DiD estimados\n")
r1 <- .extract_robust(M1B_TBFE, "dummy_holg:post_h2")
cat(sprintf("  Holgura P50 | Acumulado | FE dobles: β = %.3f (SE = %.3f)%s\n",
            r1$est, r1$se, .stars(r1$p)))
r2 <- .extract_robust(M2B_TBFE, "dummy_holg:post_h2")
cat(sprintf("  Holgura P50 | Escalonado D2020: β = %.3f (SE = %.3f)%s\n",
            r2$est, r2$se, .stars(r2$p)))


# =============================================================================
# [8] BOOTSTRAP DEL PARÁMETRO DiD (K-Means en 2 etapas)
# =============================================================================
# Para la clasificación K-Means en 2 etapas, se evalúa si el coeficiente DiD
# es robusto a cambios en la asignación tratado/control. En cada réplica
# bootstrap del clustering (sección [5]) se re-estima el modelo DiD y se
# registra el coeficiente de la interacción tratado × post_dic2020.
#
# Indicadores del resumen:
#   β mediana bootstrap : estimación central bajo reclasificación
#   Fracción mismo signo: estabilidad cualitativa del resultado
#   Fracción significativa al 5%: robustez de la inferencia estadística
# =============================================================================

estimate_did_km <- function(panel_df_in, treated_vec) {
  df_m <- panel_df_in %>%
    mutate(dummy_holg_km = as.integer(INS_COD %in% as.character(treated_vec)))
  pd <- pdata.frame(as.data.frame(df_m), index = c("INS_COD", "PERIODO"))
  fit <- tryCatch(
    plm(as.formula(paste(
      "ln_coloc1 ~ dummy_holg_km + post_h2 + dummy_holg_km:post_h2 +", controles)),
      data = pd, model = "within", effect = "twoways"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(tibble(estimate = NA_real_, p.value = NA_real_))
  ct   <- coeftest(fit, vcov = vcovHC(fit, type = "HC1", cluster = "group"))
  term <- "dummy_holg_km:post_h2"
  if (!(term %in% rownames(ct))) return(tibble(estimate = NA_real_, p.value = NA_real_))
  tibble(estimate = ct[term, 1], p.value = ct[term, 4])
}

beta_base_km <- estimate_did_km(panel_df, bancos_mas_expuestos_km)$estimate

boot_did_km <- map_dfr(sort(unique(boot_membership_long$boot_id)), function(bid) {
  assign_b  <- boot_membership_long %>%
    filter(boot_id == bid, !is.na(expuesto_boot)) %>%
    mutate(INS_COD = as.character(INS_COD))
  treated_b <- assign_b %>% filter(expuesto_boot == 1L) %>% pull(INS_COD)
  n_changed <- assign_b %>%
    left_join(bank_labels_base %>% select(INS_COD, expuesto_base), by = "INS_COD") %>%
    summarise(n = sum(expuesto_boot != expuesto_base, na.rm = TRUE)) %>% pull(n)
  res <- if (length(treated_b) == 0)
    tibble(estimate = NA_real_, p.value = NA_real_)
  else estimate_did_km(panel_df, treated_b)
  tibble(boot_id = bid, estimate = res$estimate, p.value = res$p.value,
         n_treated = length(treated_b), n_changed = n_changed,
         share_changed = n_changed / n_distinct(assign_b$INS_COD))
})

cat("✓ [8] Bootstrap DiD (K-Means 2 etapas)\n")
cat(sprintf("  Beta base = %.3f | Beta mediana = %.3f | mismo signo = %.1f%% | signif. 5%% = %.1f%%\n",
            beta_base_km,
            median(boot_did_km$estimate, na.rm = TRUE),
            mean(sign(boot_did_km$estimate) == sign(beta_base_km), na.rm = TRUE) * 100,
            mean(boot_did_km$p.value < 0.05, na.rm = TRUE) * 100))


# =============================================================================
# [9] TABLAS DE RESULTADOS
# =============================================================================

# --- Tabla 1: Estadística descriptiva por grupo y período ---
tabla_desc <- panel_df %>%
  mutate(grupo   = if_else(dummy_holg == 1L, "Tratados", "Controles"),
         periodo = if_else(PERIODO < hito_2, "Pre", "Post")) %>%
  group_by(periodo, grupo) %>%
  summarise(
    `Log(Coloc.)`  = mean(ln_coloc1,  na.rm = TRUE),
    `CET1/APR`     = mean(cet1_apr,   na.rm = TRUE),
    `Log(AT)`      = mean(ln_at,      na.rm = TRUE),
    `Depósitos/AT` = mean(deposit_at, na.rm = TRUE),
    `Efectivo/AT`  = mean(cash_at,    na.rm = TRUE),
    `ROE`          = mean(ROE,        na.rm = TRUE),
    `N Bancos`     = n_distinct(INS_COD),
    .groups = "drop"
  )

# --- Tabla 2: Regresiones principales (Holgura P50) ---
make_reg_table <- function(models_list, terms_list, row_labels, model_names) {
  
  # Bloque de coeficientes
  rows <- map2_dfr(terms_list, row_labels, function(term, label) {
    vals <- map_chr(models_list, function(m) {
      r <- .extract_robust(m, term)
      if (is.null(r)) return("")
      paste0(fmt_num(r$est, 3), .stars(r$p), "\n(", fmt_num(r$se, 3), ")")
    })
    
    tibble(
      Variable = label,
      !!!setNames(as.list(vals), model_names)
    )
  })
  
  # Bloque de metadatos
  meta <- tibble(
    Variable = c("Obs.", "R²", "Time FE", "Banco FE")
  )
  
  for (i in seq_along(models_list)) {
    m <- models_list[[i]]
    meta[[model_names[i]]] <- c(
      fmt_int(nobs(m)),
      fmt_num(as.numeric(summary(m)$r.squared["rsq"]), 2),
      "Sí",
      "Sí"
    )
  }
  
  bind_rows(rows, meta)
}

terms_holg_acum <- list("dummy_holg:post_h2")
terms_holg_esc  <- list("dummy_holg:post_h2", "dummy_holg:post",
                         "dummy_holg:post_h4", "dummy_holg:post_h5",
                         "dummy_holg:post_h6")
labels_acum     <- "D1 × T (dic. 2020)"
labels_esc      <- c("D1 × T (dic. 2020)", "D2 × T (dic. 2021)",
                      "D3 × T (dic. 2022)", "D4 × T (dic. 2023)",
                      "D5 × T (dic. 2024)")

tabla_reg_holg <- make_reg_table(
  models_list  = list(M1B_TBFE, M2B_TBFE),
  terms_list   = c(terms_holg_acum, terms_holg_esc),
  row_labels   = c(labels_acum, labels_esc),
  model_names  = c("Acumulado", "Escalonado")
)

# --- Tabla 3: Estabilidad bootstrap del clustering ---
tabla_bootstrap <- stability_bank %>%
  transmute(
    Banco             = Bancos,
    `Grupo base`      = if_else(expuesto_base == 1L, "Más expuesto", "Menos expuesto"),
    `P(más expuesto)` = round(prob_expuesto_boot, 2),
    Estabilidad       = estabilidad
  )

# --- Tabla 4: Resumen bootstrap DiD ---
tabla_boot_did <- tibble(
  Indicador = c("Réplicas válidas", "Beta base",
                 "Beta mediana bootstrap", "Intervalo P5-P95",
                 "Fracción mismo signo", "Fracción signif. al 5%"),
  Valor = c(
    fmt_int(sum(is.finite(boot_did_km$estimate))),
    fmt_num(beta_base_km, 3),
    fmt_num(median(boot_did_km$estimate, na.rm = TRUE), 3),
    paste0(fmt_num(quantile(boot_did_km$estimate, 0.05, na.rm = TRUE), 3),
           " a ",
           fmt_num(quantile(boot_did_km$estimate, 0.95, na.rm = TRUE), 3)),
    paste0(fmt_num(mean(sign(boot_did_km$estimate) == sign(beta_base_km),
                        na.rm = TRUE) * 100, 1), "%"),
    paste0(fmt_num(mean(boot_did_km$p.value < 0.05, na.rm = TRUE) * 100, 1), "%")
  )
)

cat("✓ [9] Tablas generadas\n")


# =============================================================================
# [10] FIGURAS
# =============================================================================

CMF_PURPLE      <- "#5B2C83"
CMF_PURPLE_PALE <- "#D4B8E8"
CMF_TEAL        <- "#20B2AA"
CMF_GRAY        <- "#6E6E6E"
CMF_GRAY_LIGHT  <- "#D0D0D0"

theme_cmf <- theme_minimal(base_size = 11) +
  theme(
    text             = element_text(color = CMF_GRAY),
    axis.text        = element_text(size = 9,  color = CMF_GRAY),
    axis.title       = element_text(size = 10, color = CMF_GRAY),
    panel.grid.major = element_line(color = CMF_GRAY_LIGHT, linewidth = 0.3),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(size = 9, color = CMF_GRAY),
    legend.title     = element_blank(),
    legend.text      = element_text(size = 9),
    plot.title       = element_text(size = 12, face = "bold",
                                    hjust = 0.5, color = CMF_PURPLE),
    plot.subtitle    = element_text(size = 10, hjust = 0.5, color = CMF_GRAY),
    plot.margin      = margin(8, 12, 8, 12)
  )

vlines_hitos <- list(
  geom_vline(xintercept = hito_0, linetype = "dashed", color = CMF_GRAY, linewidth = 0.4),
  geom_vline(xintercept = hito_2, linetype = "dashed", color = CMF_GRAY, linewidth = 0.4),
  annotate("text", x = hito_0, y = Inf, label = "OCT-2018",
           angle = 90, vjust = -0.3, hjust = 1.1, size = 2.5, color = CMF_GRAY),
  annotate("text", x = hito_2, y = Inf, label = "DIC-2020",
           angle = 90, vjust = -0.3, hjust = 1.1, size = 2.5, color = CMF_GRAY)
)

# Figura 1 — Evolución CET1/APR industria
fig1 <- panel_df %>%
  group_by(PERIODO) %>%
  summarise(cet1_med = mean(cet1_apr, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = PERIODO, y = cet1_med)) +
  geom_line(color = CMF_TEAL, linewidth = 0.8) +
  vlines_hitos +
  labs(title    = "Evolución del Capital Regulatorio — Industria",
       subtitle = "CET1/APR promedio (12 bancos en muestra)",
       x = NULL, y = "CET1/APR (%)") +
  theme_cmf

# Figura 2 — Tendencias paralelas: Log(Colocaciones) por grupo
fig2 <- panel_df %>%
  mutate(grupo = if_else(dummy_holg == 1L, "Tratados", "Control")) %>%
  group_by(PERIODO, grupo) %>%
  summarise(lc = mean(ln_coloc1, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = PERIODO, y = lc, color = grupo)) +
  geom_line(linewidth = 0.8) +
  vlines_hitos +
  scale_color_manual(values = c("Tratados" = CMF_PURPLE, "Control" = CMF_TEAL)) +
  labs(title    = "Evolución de Colocaciones por Grupo",
       subtitle = "Media de Log(Colocaciones+1) — Holgura P50",
       x = NULL, y = "Log(Colocaciones+1)") +
  theme_cmf + theme(legend.position = "bottom")

# Figura 3 — Clasificación K-Means en 2 etapas
fig3 <- ggplot(kmeans_2stage,
               aes(x = log_AT, y = RWA_density, color = grupo_2stage,
                   shape = factor(cluster_estructural))) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(aes(label = Bancos), size = 3, max.overlaps = 15, seed = 42) +
  scale_color_manual(values = c("Más expuesto"  = CMF_PURPLE,
                                 "Menos expuesto" = CMF_TEAL)) +
  labs(title    = "Clasificación K-Means en 2 Etapas",
       subtitle = "Forma: tipología estructural (Etapa 1) | Color: exposición relativa (Etapa 2)",
       x = "Log(Activos Totales)", y = "Densidad APR (APR/AT)") +
  theme_cmf + theme(legend.position = "bottom")

# Figura 4 — Bootstrap DiD: distribución del coeficiente
fig4 <- boot_did_km %>%
  filter(is.finite(estimate)) %>%
  ggplot(aes(x = estimate)) +
  geom_histogram(bins = 30, fill = CMF_PURPLE_PALE, color = CMF_PURPLE) +
  geom_vline(xintercept = beta_base_km,
             linetype = "dashed", color = CMF_TEAL, linewidth = 0.8) +
  labs(title    = "Distribución Bootstrap del Coeficiente DiD",
       subtitle = "K-Means 2 etapas | Línea: estimación base",
       x = "β (tratado × post dic. 2020)", y = "Réplicas bootstrap") +
  theme_cmf

# Figura 5 — Bootstrap DiD: cambio de grupos vs cambio del coeficiente
fig5 <- boot_did_km %>%
  filter(is.finite(estimate)) %>%
  ggplot(aes(x = share_changed, y = estimate)) +
  geom_point(alpha = 0.55, color = CMF_PURPLE) +
  geom_smooth(method = "lm", se = FALSE, color = CMF_TEAL, linewidth = 0.7) +
  geom_hline(yintercept = beta_base_km,
             linetype = "dashed", color = CMF_GRAY, linewidth = 0.5) +
  labs(title    = "Sensibilidad: Cambio de Grupos vs. Cambio del Coeficiente",
       subtitle = "Cada punto es una réplica bootstrap",
       x = "Proporción de bancos que cambia de grupo",
       y = "β estimado") +
  theme_cmf

cat("✓ [10] Figuras generadas: fig1, fig2, fig3, fig4, fig5\n")

# =============================================================================
# EXPORTAR PARA EL QMD DE PRESENTACIÓN
# =============================================================================
# El QMD hace source("analysis_basilea_iii.R") y usa directamente:
#   panel_df, kmeans_2stage, stability_bank, boot_did_km, beta_base_km
#   tabla_desc, tabla_reg_holg, tabla_bootstrap, tabla_boot_did
#   fig1, fig2, fig3, fig4, fig5
#   theme_cmf, CMF_PURPLE, CMF_TEAL, CMF_GRAY, CMF_GRAY_LIGHT, CMF_PURPLE_PALE
#   hito_0, hito_2, cutoff_date, hito_4, hito_5, hito_6
# =============================================================================

cat("\n", strrep("=", 65), "\n")
cat(" PIPELINE COMPLETADO — Objetos listos para la presentación\n")
cat(strrep("=", 65), "\n")
