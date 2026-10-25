# Microbioma núcleo y composición por intrusión y sexo.

ps_saliva_g <- tax_glom(ps_saliva_u, taxrank = "Genus", NArm = TRUE)

otu_sal_g <- otu_taxa_rows(ps_saliva_g)

prev_core <- rowMeans(otu_sal_g > 0)

core90_ids <- names(prev_core)[
  prev_core >= 0.90
]

core95_ids <- names(prev_core)[
  prev_core >= 0.95
]

metricas_core <- function(
  ps_obj,
  core_ids,
  nombre_core
) {
  pg <- tax_glom(ps_obj, taxrank = "Genus", NArm = TRUE)

  ids <- intersect(core_ids, taxa_names(pg))

  o <- otu_taxa_rows(pg)

  if (
    length(ids) == 0
  ) {
    stop("No hay taxa core tras intersectar para ", nombre_core)
  }

  present <- colSums(o[ ids, , drop = FALSE ] > 0)

  core_reads <- colSums(o[ ids, , drop = FALSE ])

  total <- colSums(o)

  tibble(
    id_muestra = colnames(o),
    core = nombre_core,
    n_core_total = length(ids),
    n_core_present = as.numeric(present),
    frac_core_present = as.numeric(present) /
      length(ids),
    prop_reads_core = as.numeric(core_reads) /
      as.numeric(total)
  )
}

core_metrics <- bind_rows(
  metricas_core(ps_amb_u_m, core90_ids, "core90"),
  metricas_core(ps_amb_u_m, core95_ids, "core95")
) %>%
  left_join(meta_base, by = "id_muestra")

core_outcomes <- c("n_core_present", "prop_reads_core")

ajustar_core <- function(
  d,
  outcome,
  exposure
) {
  z <- ajustar_alpha_cont(outcome, exposure, data = d, covars = covars_alpha)

  if (is.null(z)) {
    return(NULL)
  }

  z
}

core_global <- bind_rows(
  lapply(
    c("core90", "core95"),
    function(cx) {
      d <- core_metrics %>%
        filter(core == cx)

      bind_rows(
        lapply(
          c("miteco10", "miteco_days10_z"),
          function(e) {
            bind_rows(
              lapply(
                core_outcomes,
                function(o) {
                  z <- ajustar_core(d, o, e)

                  if (is.null(z)) {
                    return(NULL)
                  }

                  z %>%
                    mutate(core = cx, .before = 1)
                }
              )
            )
          }
        )
      )
    }
  )
) %>%
  group_by(core, exposure) %>%
  mutate(p_BH_outcomes = p.adjust(p, method = "BH")) %>%
  ungroup()

core_int_sexo <- bind_rows(
  lapply(
    c("core90", "core95"),
    function(cx) {
      d <- core_metrics %>%
        filter(core == cx)

      bind_rows(
        lapply(
          core_outcomes,
          function(o) {
            z <- ajustar_alpha_interaccion(o, "miteco10", "sex", data = d)

            if (is.null(z)) {
              return(NULL)
            }

            z %>%
              mutate(core = cx, .before = 1)
          }
        )
      )
    }
  )
) %>%
  group_by(core) %>%
  mutate(p_BH_outcomes = p.adjust(p_interaccion, method = "BH")) %>%
  ungroup()

# Composicion taxonomica descriptiva.

ps_gen_desc <- tax_glom(ps_amb_u_m, taxrank = "Genus", NArm = TRUE)

otu_desc <- otu_taxa_rows(ps_gen_desc)

tt_desc <- tax_df(ps_gen_desc)

rel_desc <- sweep(otu_desc, 2, colSums(otu_desc), "/")

genus_names <- tt_desc$Genus
genus_names[
  is.na(genus_names)
] <- "Genus_no_asignado"

genus_names <- make.unique(genus_names, sep = "__")

rownames(rel_desc) <- genus_names

ra_mean <- rowMeans(rel_desc)

top15 <- names(sort(ra_mean, decreasing = TRUE))[
  seq_len(min(15, length(ra_mean)))
]

comp_long <- as.data.frame(t(rel_desc[ top15, , drop = FALSE ]))

comp_long$id_muestra <- rownames(comp_long)

comp_long <- comp_long %>%
  left_join(meta_base %>% select(id_muestra, miteco10_factor, sex), by = "id_muestra") %>%
  pivot_longer(cols = all_of(top15), names_to = "Genus", values_to = "abundancia_relativa")

comp_resumen <- comp_long %>%
  group_by(miteco10_factor, sex, Genus) %>%
  summarise(
    n = n(),
    media = mean(abundancia_relativa, na.rm = TRUE),
    mediana = median(abundancia_relativa, na.rm = TRUE),
    .groups = "drop"
  )
