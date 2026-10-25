# Exportación.

write_csv(tabla_cohorte, file.path(TAB_DIR, "01_caracteristicas_cohorte_segun_intrusion.csv"), na = "")

write_csv(datos_faltantes, file.path(TAB_DIR, "01b_datos_faltantes_tabla_cohorte.csv"), na = "")

write_csv(resumen_ventanas, file.path(TAB_DIR, "02_frecuencia_intrusion_por_ventanas.csv"), na = "")

write_csv(tabla_concordancia, file.path(TAB_DIR, "03_concordancia_intrusion_oficial_PM10.csv"), na = "")

write_csv(env_cluster, file.path(TAB_DIR, "04_exposicion_por_cluster_isla_fecha.csv"), na = "")

write_csv(
  correlaciones_ambientales,
  file.path(TAB_DIR, "05_correlaciones_intrusion_PM10_cluster.csv"),
  na = ""
)

write_csv(alpha_pm10_10d, file.path(TAB_DIR, "06_alpha_PM10_D1_D10.csv"), na = "")

write_csv(alpha_pm10_x_sexo, file.path(TAB_DIR, "07_alpha_PM10_interaccion_sexo.csv"), na = "")

write_csv(alpha_pm10_sexo_estrat, file.path(TAB_DIR, "08_alpha_PM10_estratificado_sexo.csv"), na = "")

write_csv(
  predicciones_sexo,
  file.path(TAB_DIR, "09_estimaciones_ajustadas_intrusion_sexo.csv"),
  na = ""
)

write_csv(alpha_pm10_180d, file.path(TAB_DIR, "10_alpha_PM10_180d.csv"), na = "")

write_csv(beta_exposiciones, file.path(TAB_DIR, "11_beta_exposiciones_principales.csv"), na = "")

write_csv(permanova_todas, file.path(TAB_DIR, "12_PERMANOVA_todos_terminos.csv"), na = "")

write_csv(
  permanova_bray_miteco,
  file.path(TAB_DIR, "13_PERMANOVA_Bray_covariables_intrusion.csv"),
  na = ""
)

write_csv(
  permanova_bray_pm10,
  file.path(TAB_DIR, "14_PERMANOVA_Bray_covariables_PM10_media.csv"),
  na = ""
)

write_csv(betadisper_miteco, file.path(TAB_DIR, "15_PERMDISP_Bray_intrusion.csv"), na = "")

write_csv(beta_cluster, file.path(TAB_DIR, "16_sensibilidad_Bray_cluster_isla_fecha.csv"), na = "")

write_csv(exac_pm10_180d, file.path(TAB_DIR, "17_PM10_180d_y_fenotipo_exacerbador.csv"), na = "")

write_csv(
  comparacion_metodos_completa,
  file.path(TAB_DIR, "18_ANCOMBC2_DESeq2_comparacion_completa.csv"),
  na = ""
)

write_csv(
  comparacion_metodos_seleccion,
  file.path(TAB_DIR, "19_ANCOMBC2_DESeq2_senales_seleccionadas.csv"),
  na = ""
)

write_csv(resumen_metodos, file.path(TAB_DIR, "20_ANCOMBC2_DESeq2_resumen.csv"), na = "")

write_csv(top15_genus, file.path(TAB_DIR, "21_top15_generos.csv"), na = "")

write_csv(prevalencia_plot, file.path(TAB_DIR, "22_prevalencia_y_core_generos.csv"), na = "")

write_csv(resumen_core, file.path(TAB_DIR, "23_resumen_microbioma_nucleo.csv"), na = "")

write_csv(
  tabla_indicadores_ambientales,
  file.path(TAB_DIR, "25_tabla_integrada_indicadores_ambientales.csv"),
  na = ""
)

write_csv(qc_constancia_cluster, file.path(TAB_DIR, "26_QC_constancia_exposicion_cluster.csv"), na = "")

manifest <- registrar(devolver = TRUE)

write_csv(manifest, file.path(TAB_DIR, "27_manifest_ejecucion.csv"), na = "")

# Libro Excel de apoyo para la redacción.
writexl::write_xlsx(
  list(
    Cohorte = tabla_cohorte,
    Faltantes = datos_faltantes,
    Ventanas = resumen_ventanas,
    Concordancia = tabla_concordancia,
    Correlaciones = correlaciones_ambientales,
    Alpha_PM10_10d = alpha_pm10_10d,
    PM10_x_sexo = alpha_pm10_x_sexo,
    PM10_sexo_estrat = alpha_pm10_sexo_estrat,
    Ajustadas_sexo = predicciones_sexo,
    Alpha_PM10_180d = alpha_pm10_180d,
    Beta_exposiciones = beta_exposiciones,
    PERMANOVA_Bray_intr = permanova_bray_miteco,
    PERMANOVA_Bray_PM10 = permanova_bray_pm10,
    PERMDISP = betadisper_miteco,
    Beta_cluster = beta_cluster,
    PM10_180_exac = exac_pm10_180d,
    Taxa_seleccionados = comparacion_metodos_seleccion,
    Metodos_resumen = resumen_metodos,
    Top15_generos = top15_genus,
    Core_resumen = resumen_core,
    Indicadores = tabla_indicadores_ambientales,
    Manifest = manifest
  ),
  file.path(TAB_DIR, "28_RESUMEN_CIERRE_RESULTADOS_TFM.xlsx")
)

saveRDS(
  list(
    meta = meta,
    env_cluster = env_cluster,
    distancias = distancias,
    alpha_pm10_10d = alpha_pm10_10d,
    alpha_pm10_x_sexo = alpha_pm10_x_sexo,
    alpha_pm10_180d = alpha_pm10_180d,
    permanova_todas = permanova_todas,
    comparacion_metodos = comparacion_metodos_completa,
    predicciones_sexo = predicciones_sexo
  ),
  file.path(OBJ_DIR, "objetos_cierre_resultados_TFM.rds")
)

capture.output(sessionInfo(), file = file.path(LOG_DIR, "SESSION_INFO.txt"))

message("Tablas y complementos PM10 guardados.")
