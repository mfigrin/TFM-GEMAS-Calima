# Guardar las tablas del análisis MITECO.
# No se exportan resultados reemplazados ni estados comunitarios ajenos al TFM final.

write_csv(metricas, file.path(ruta_tablas, "01_dataset_metricas_MITECO_analitico.csv"), na = "")

write_csv(alpha_global, file.path(ruta_tablas, "02_alpha_global_MITECO.csv"), na = "")

write_csv(alpha_factor, file.path(ruta_tablas, "03_alpha_grado_patron_MITECO.csv"), na = "")

write_csv(alpha_int_sexo, file.path(ruta_tablas, "04_alpha_interaccion_MITECO_sexo.csv"), na = "")

write_csv(alpha_sexo_estrat, file.path(ruta_tablas, "05_alpha_estratificado_sexo.csv"), na = "")

write_csv(alpha_int_isla, file.path(ruta_tablas, "06_alpha_interaccion_MITECO_isla.csv"), na = "")

write_csv(alpha_int_exac, file.path(ruta_tablas, "07_alpha_interaccion_MITECO_exac6m.csv"), na = "")

write_csv(
  alpha_int_regimen,
  file.path(ruta_tablas, "08_alpha_interaccion_MITECO_winter_summer.csv"),
  na = ""
)

write_csv(alpha_ventanas, file.path(ruta_tablas, "09_alpha_ventanas_D1_10_D11_30_D31_90.csv"), na = "")

write_csv(
  alpha_ventanas_sexo,
  file.path(ruta_tablas, "10_alpha_ventanas_interaccion_sexo.csv"),
  na = ""
)

write_csv(spline_res, file.path(ruta_tablas, "11_alpha_dosis_respuesta_splines_MITECO.csv"), na = "")

write_csv(beta_global, file.path(ruta_tablas, "12_beta_global_MITECO.csv"), na = "")

write_csv(beta_int_sexo, file.path(ruta_tablas, "13_beta_interaccion_MITECO_sexo.csv"), na = "")

write_csv(beta_int_isla, file.path(ruta_tablas, "14_beta_interaccion_MITECO_isla.csv"), na = "")

write_csv(beta_int_exac, file.path(ruta_tablas, "15_beta_interaccion_MITECO_exac6m.csv"), na = "")

write_csv(
  beta_int_regimen,
  file.path(ruta_tablas, "16_beta_interaccion_MITECO_winter_summer.csv"),
  na = ""
)

write_csv(beta_ventanas, file.path(ruta_tablas, "17_beta_ventanas_temporales_MITECO.csv"), na = "")

write_csv(core_global, file.path(ruta_tablas, "20_core_microbiome_global_MITECO.csv"), na = "")

write_csv(core_int_sexo, file.path(ruta_tablas, "21_core_microbiome_interaccion_sexo.csv"), na = "")

write_csv(comp_resumen, file.path(ruta_tablas, "22_composicion_top15_generos_MITECO_sexo.csv"), na = "")

write_csv(ancom_genus_all, file.path(ruta_tablas, "23_ANCOMBC2_Genus_todos.csv"), na = "")

write_csv(ancom_genus_resumen, file.path(ruta_tablas, "24_ANCOMBC2_Genus_resumen.csv"), na = "")

write_csv(
  ancom_genus_all %>%
    filter(robusto),
  file.path(ruta_tablas, "25_ANCOMBC2_Genus_robustos.csv"),
  na = ""
)

write_csv(
  ancom_sens_all,
  file.path(ruta_tablas, "26_ANCOMBC2_sensibilidad_taxonomica_todos.csv"),
  na = ""
)

write_csv(
  ancom_sens_resumen,
  file.path(ruta_tablas, "27_ANCOMBC2_sensibilidad_taxonomica_resumen.csv"),
  na = ""
)

write_csv(exac_alpha, file.path(ruta_tablas, "28_exacerbacion6M_alpha_ajustada_MITECO.csv"), na = "")

write_csv(exac_miteco180, file.path(ruta_tablas, "29_exacerbacion6M_carga_MITECO180.csv"), na = "")

registrar("Diversidad alfa", "OK", "Global, interacciones, estratos y ventanas.")
registrar("Beta secundaria", "OK", "999 permutaciones; principal en 06.")
registrar("Microbioma núcleo", "OK", "Umbrales 90 % y 95 %.")
registrar("Exacerbación 6M", "OK", "Diversidad y carga ambiental a 180 días.")
write_csv(manifest, file.path(ruta_tablas, "32_manifest_ejecucion.csv"), na = "")
if (any(manifest$estado == "ERROR")) stop("Un bloque registró ERROR. Revisar 32_manifest_ejecucion.csv.")
capture.output(sessionInfo(), file = file.path(ruta_salida, "sessionInfo.txt"))
