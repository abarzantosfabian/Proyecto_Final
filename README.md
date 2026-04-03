# Evaluación de Impacto de Basilea III en Colocaciones Bancarias

## 1. Introducción
Este repositorio contiene un pipeline reproducible en R para evaluar el impacto de la implementación de Basilea III sobre las colocaciones bancarias en Chile. El análisis se basa en un enfoque econométrico de diferencias-en-diferencias (DiD) aplicado a datos de panel.

---

## 2. Problema de Identificación
La regulación de Basilea III afecta a toda la industria bancaria, lo que implica la ausencia de un grupo de control natural. Para resolver este problema, se construyen grupos de tratamiento y control en base a exposición relativa pre-reforma.

---

## 3. Estrategia Empírica

### Definición de Tratamiento
Se consideran tres enfoques:

- **Baseline:** bancos con menor holgura de capital (bajo la mediana)
- **Robustez 1:** percentil 25 de holgura
- **Robustez 2:** clustering K-Means en 2 etapas

---

## 4. Modelo Econométrico

Se estima el siguiente modelo:

ln(Colocaciones_it) = α_i + θ_t + β (Tratamiento_i × Post_t) + γX_it + ε_it

Donde:

- α_i: efectos fijos por banco  
- θ_t: efectos fijos temporales  
- X_it: controles (tamaño, liquidez, rentabilidad, macro)  

---

## 5. Clustering en 2 Etapas

El enfoque propuesto separa:

### Etapa 1: Tipología estructural
Clasificación de bancos según:
- tamaño (log activos)
- composición de cartera
- densidad de riesgo

### Etapa 2: Exposición relativa
Asignación de tratamiento dentro de cada tipología según holgura.

**Ventaja clave:** evita comparar bancos estructuralmente distintos.

---

## 6. Pipeline

El flujo del análisis es:

1. Construcción del panel banco-tiempo  
2. Cálculo de holgura regulatoria  
3. Clasificación (K-Means 2 etapas)  
4. Validación (bootstrap)  
5. Estimación DiD  
6. Generación de tablas y gráficos  
7. Renderización en Quarto (Beamer)

---

## 7. Validación

Se implementan dos estrategias:

- **Bootstrap de clustering:** estabilidad de asignación  
- **Bootstrap de coeficiente DiD:** sensibilidad a clasificación  

---

## 8. Resultados Principales

- No se observan efectos estadísticamente significativos en colocaciones  
- Magnitud del efecto cercana a cero  
- Resultados robustos a distintas definiciones de tratamiento  

---

## 9. Limitaciones

- Datos sintéticos (en versión actual)  
- Tamaño muestral reducido  
- Supuesto de tendencias paralelas  

---

## 10. Reproducibilidad

El proyecto está diseñado para ser completamente reproducible:

- Código en R  
- Pipeline modular  
- Integración con Quarto para reportes automáticos  

---

## 11. Uso

Para ejecutar el análisis:

```r
source("analysis_basilea_iii.R")
```

Para generar presentación:

```bash
quarto render presentacion.qmd
```

---

## 12. Conclusión

La evidencia sugiere que la implementación de Basilea III no tuvo un impacto inmediato significativo en las colocaciones bancarias. Los efectos, de existir, serían graduales y heterogéneos.

---

## Autor
Fabián Abarza
