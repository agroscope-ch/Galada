#Load libraries needed

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(gridExtra)


# ============================================================
# Figure 6: Ergosterol by treatment and pathogen
# Four independent biological replicates
# Primary analysis: prespecified Wilcoxon rank-sum tests
# Multiple testing: Holm adjustment within each pathogen
# ============================================================


# Output folder (created in the working directory)
output_dir <- "Figure"
dir.create(output_dir, showWarnings = FALSE)

# 1. Import data ---------------------------------

# ------------------------------------------------------------
# Data source: Zenodo (https://zenodo.org/records/24052059)
# One wide table: chitin, glucan and ergosterol per Number x Treatment x Pathogen.
# Empty cells mean "not measured" for that variable.
# ------------------------------------------------------------
zenodo_url  <- "https://zenodo.org/records/24052059/files/2026_chitin_glucan_ergosterol.csv?download=1"
zenodo_file <- tempfile(fileext = ".csv")
download.file(zenodo_url, zenodo_file, mode = "wb", quiet = TRUE)
raw_data_zenodo <- read.table(zenodo_file, header = TRUE, sep = ";", dec = ".")


raw_data_figure_6 <- raw_data_zenodo |>
  select(Treatment, Pathogen, `Ergosterol..ug.mg.DW.`) |>
  rename(
    Ergosterol = `Ergosterol..ug.mg.DW.`
  ) |>
  mutate(
    Treatment = trimws(as.character(Treatment)),
    Pathogen  = trimws(as.character(Pathogen)),
    Ergosterol = as.numeric(Ergosterol)
  )

# Define a consistent order for plotting and reporting.

treatment_order <- c(
  "PDB",
  "PDB+MeOH",
  "Frangula",
  "Rheum",
  "Curcuma"
)

pathogen_order <- c("C. coccodes", "H. solani")


raw_data_figure_6 <- raw_data_figure_6 |>
  filter(!is.na(Ergosterol)) |>
  mutate(
    Treatment = factor(Treatment, levels = treatment_order),
    Pathogen = factor(Pathogen, levels = pathogen_order)
  ) |>
  droplevels()


# 2. Descriptive summaries ------------------------------------

summary_table <- raw_data_figure_6 |>
  group_by(Pathogen, Treatment) |>
  summarise(
    n = n(),
    mean = mean(Ergosterol),
    sd = sd(Ergosterol),
    median = median(Ergosterol),
    q1 = quantile(Ergosterol, 0.25, names = FALSE),
    q3 = quantile(Ergosterol, 0.75, names = FALSE),
    min = min(Ergosterol),
    max = max(Ergosterol),
    .groups = "drop"
  )

print(summary_table)

# 3. Prespecified comparisons ---------------------------------

planned_comparisons <- list(
  c("PDB", "PDB+MeOH"),
  c("PDB", "Frangula"),
  c("PDB+MeOH", "Rheum"),
  c("PDB+MeOH", "Curcuma")
)

# Run independent two-sample Wilcoxon rank-sum tests.

run_test <- function(raw_data_figure_6, pathogen_name, group1, group2) {
  
  x <- raw_data_figure_6 |>
    filter(
      Pathogen == pathogen_name,
      Treatment == group1
    ) |>
    pull(Ergosterol)
  
  y <- raw_data_figure_6 |>
    filter(
      Pathogen == pathogen_name,
      Treatment == group2
    ) |>
    pull(Ergosterol)
  
  if (length(x) == 0 || length(y) == 0) {
    return(tibble(
      Pathogen = pathogen_name,
      Comparison = paste(group1, "vs", group2),
      n1 = length(x),
      n2 = length(y),
      W = NA_real_,
      p_value = NA_real_
    ))
  }
  
  test <- wilcox.test(
    x, y,
    alternative = "two.sided",
    exact = NULL
  )
  
  tibble(
    Pathogen = pathogen_name,
    Comparison = paste(group1, "vs", group2),
    n1 = length(x),
    n2 = length(y),
    W = unname(test$statistic),
    p_value = test$p.value
  )
}

pathogens_present <- unique(as.character(raw_data_figure_6$Pathogen))

test_results <- bind_rows(
  lapply(pathogens_present, function(pathogen_name) {
    bind_rows(
      lapply(planned_comparisons, function(comparison) {
        run_test(
          raw_data_figure_6,
          pathogen_name,
          comparison[1],
          comparison[2]
        )
      })
    )
  })
) |>
  group_by(Pathogen) |>
  mutate(
    p_holm = p.adjust(p_value, method = "holm")
  ) |>
  ungroup()

print(test_results)


# ------------------------------------------------------------
# Common theme figure
# ------------------------------------------------------------

theme_ergo <- theme_classic(base_size = 14) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    
    axis.line = element_line(colour = "black", linewidth = 0.7),
    axis.ticks = element_line(colour = "black", linewidth = 0.7),
    
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 14),
    
    plot.title = element_text(
      face = "italic",
      size = 16,
      hjust = 0.5,
      colour = "black"
    ),
    
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )

# ------------------------------------------------------------
# Prepare individual observations
# ------------------------------------------------------------

raw_plot_data <- raw_data_figure_6 |>
  filter(!is.na(Ergosterol)) |>
  mutate(
    Treatment = factor(
      Treatment,
      levels = c(
        "PDB",
        "PDB+MeOH",
        "Frangula",
        "Rheum",
        "Curcuma"
      )
    )
  ) # we just relevel, without droping any level

# ------------------------------------------------------------
# C. coccodes — Panel A
# ------------------------------------------------------------

plot_ergo_cc <- ggplot(
  data = filter(summary_table, Pathogen == "C. coccodes"),
  aes(x = Treatment, y = mean)
) +
  
  # Bars
  geom_col(
    fill = "#599AD3",
    width = 0.70
  ) +
  
  # Standard deviation
  geom_errorbar(
    aes(
      ymin = mean - sd,
      ymax = mean + sd
    ),
    width = 0.20,
    linewidth = 0.9,
    colour = "#4D4D4D"
  ) +
  
  # Labels
  labs(
    title = "",
    tag = "A",
    x = "Treatment",
    y = "Ergosterol (µg/mg DW)"
  ) +
  
  # Same scale for both pathogens
  scale_x_discrete(
    labels = c(
      "PDB" = "PDB",
      "PDB+MeOH" = "PDB+MeOH",
      "Frangula" = expression(italic("Frangula")),
      "Rheum" = expression(italic("Rheum")),
      "Curcuma" = expression(italic("Curcuma"))
    )
  ) +
scale_y_continuous(
    limits = c(0, 150),
    breaks = seq(0, 150, 25),
    expand = c(0, 0)
  ) +
  
  theme_ergo  +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 16,
      colour = "black"
    ),
    plot.tag.position = c(0, 1),
    plot.title.position = "plot"
  )

# ------------------------------------------------------------
# H. solani — Panel B
# ------------------------------------------------------------

plot_ergo_hs <- ggplot(
  data = filter(summary_table, Pathogen == "H. solani"),
  aes(x = Treatment, y = mean)
) +
  
  # Bars
  geom_col(
    fill = "#599AD3",
    width = 0.70
  ) +
  
  # Standard deviation
  geom_errorbar(
    aes(
      ymin = mean - sd,
      ymax = mean + sd
    ),
    width = 0.20,
    linewidth = 0.9,
    colour = "#4D4D4D"
  ) +
  
  # Labels
  labs(
    title = "",
    tag = "B",
    x = "Treatment",
    y = "Ergosterol (µg/mg DW)"
  ) +
  scale_x_discrete(
    labels = c(
      "PDB" = "PDB",
      "PDB+MeOH" = "PDB+MeOH",
      "Frangula" = expression(italic("Frangula")),
      "Rheum" = expression(italic("Rheum")),
      "Curcuma" = expression(italic("Curcuma"))
    )
  ) +
  
  # Same scale for both pathogens
  scale_y_continuous(
    limits = c(0, 150),
    breaks = seq(0, 150, 25),
    expand = c(0, 0)
  ) +
  
  theme_ergo +
  
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 16,
      colour = "black"
    ),
    plot.tag.position = c(0, 1),
    plot.title.position = "plot"
  )


# ------------------------------------------------------------
# Combine figures
# ------------------------------------------------------------

figure6 <- grid.arrange(
  plot_ergo_cc,
  plot_ergo_hs,
  ncol = 2
)

ggsave(
  file.path(output_dir, "figure_6.png"),
  plot = figure6,
  width = 12,
  height = 5.5,
  units = "in",
  dpi = 400
)


# ============================================================
# FIGURE 7 — Chitin and glucan equivalents
# ============================================================

# ------------------------------------------------------------
# 1. Import and prepare data
# ------------------------------------------------------------

raw_data_figure_7 <- raw_data_zenodo |>
  # keep only rows with a chitin and/or glucan measurement
  filter(!(is.na(chitin.ug.mg) & is.na(glucan.ug.mg))) |>
  mutate(
    Pathogen = trimws(as.character(Pathogen)),
    Treatment = trimws(as.character(Treatment))
  )

# Standardize pathogen names before converting to factors
raw_data_figure_7 <- raw_data_figure_7 |>
  mutate(
    Pathogen = recode(
      Pathogen,
      "C coccodes" = "C. coccodes",
      "H solani" = "H. solani"
    )
  )

treatment_order_figure_7 <- c(
  "PDB", "PDB+MeOH", "Frangula", "Rheum", "Curcuma"
)

pathogen_order_figure_7 <- c(
  "C. coccodes", "H. solani"
)

raw_data_figure_7 <- raw_data_figure_7 |>
  mutate(
    Treatment = factor(
      Treatment,
      levels = treatment_order_figure_7
    ),
    Pathogen = factor(
      Pathogen,
      levels = pathogen_order_figure_7
    )
  )

raw_data_figure_7 <- raw_data_figure_7 |>
  mutate(
    glucan.ug.mg = as.numeric(glucan.ug.mg),
    chitin.ug.mg = as.numeric(chitin.ug.mg)
  )


# ------------------------------------------------------------
# 2. Select observations for the four panels
# ------------------------------------------------------------

# C. coccodes: 
data_figure_7_cc <- raw_data_figure_7 |>
  filter(Pathogen == "C. coccodes")

# H. solani:
data_figure_7_hs <- raw_data_figure_7 |>
  filter(
    Pathogen == "H. solani" )

# ------------------------------------------------------------
# 3. Calculate means and standard deviations
# ------------------------------------------------------------

summarise_figure_7 <- function(data) {
  data |>
    group_by(Treatment, .drop = FALSE) |>
    summarise(
      n_glucan = sum(!is.na(glucan.ug.mg)),
      mean_glucan = if (n_glucan > 0) {
        mean(glucan.ug.mg, na.rm = TRUE)
      } else {
        NA_real_
      },
      sd_glucan = if (n_glucan > 1) {
        sd(glucan.ug.mg, na.rm = TRUE)
      } else {
        NA_real_
      },
      n_chitin = sum(!is.na(chitin.ug.mg)),
      mean_chitin = if (n_chitin > 0) {
        mean(chitin.ug.mg, na.rm = TRUE)
      } else {
        NA_real_
      },
      sd_chitin = if (n_chitin > 1) {
        sd(chitin.ug.mg, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    )
}

summary_figure_7_cc <- summarise_figure_7(data_figure_7_cc)

summary_figure_7_hs <- summarise_figure_7(data_figure_7_hs)


# ============================================================
# 4. Statistical analysis — Figure 7
# ============================================================

# Four planned comparisons, matching Figure 6
planned_comparisons_figure_7 <- list(
  "PDB vs PDB+MeOH" = c("PDB", "PDB+MeOH"),
  "PDB vs Frangula" = c("PDB", "Frangula"),
  "PDB+MeOH vs Rheum" = c("PDB+MeOH", "Rheum"),
  "PDB+MeOH vs Curcuma" = c("PDB+MeOH", "Curcuma")
)

# Function to perform planned Wilcoxon rank-sum tests
run_wilcoxon_figure_7 <- function(
    data,
    outcome,
    pathogen_label,
    outcome_label) {
  
  results <- lapply(
    names(planned_comparisons_figure_7),
    function(comparison_name) {
      
      groups <- planned_comparisons_figure_7[[comparison_name]]
      
      group_1 <- data |>
        filter(Treatment == groups[1]) |>
        pull(all_of(outcome)) |>
        na.omit()
      
      group_2 <- data |>
        filter(Treatment == groups[2]) |>
        pull(all_of(outcome)) |>
        na.omit()
      
      # Record sample sizes and test only when both groups contain data
      p_value <- if (length(group_1) > 0 && length(group_2) > 0) {
        wilcox.test(
          group_1,
          group_2,
          exact = FALSE
        )$p.value
      } else {
        NA_real_
      }
      
      tibble(
        Pathogen = pathogen_label,
        Outcome = outcome_label,
        Comparison = comparison_name,
        n1 = length(group_1),
        n2 = length(group_2),
        p_value = p_value
      )
    }
  ) |>
    bind_rows() |>
    mutate(
      p_adjusted = p.adjust(p_value, method = "holm")
    )
  
  results
}

# Run tests for glucan and chitin in C. coccodes
stats_glucan_cc_figure_7 <- run_wilcoxon_figure_7(
  data = data_figure_7_cc,
  outcome = "glucan.ug.mg",
  pathogen_label = "C. coccodes",
  outcome_label = "Glucan"
)

stats_chitin_cc_figure_7 <- run_wilcoxon_figure_7(
  data = data_figure_7_cc,
  outcome = "chitin.ug.mg",
  pathogen_label = "C. coccodes",
  outcome_label = "Chitin"
)

# Run tests for glucan and chitin in H. solani spores
stats_glucan_hs_figure_7 <- run_wilcoxon_figure_7(
  data = data_figure_7_hs,
  outcome = "glucan.ug.mg",
  pathogen_label = "H. solani (spores)",
  outcome_label = "Glucan"
)

stats_chitin_hs_figure_7 <- run_wilcoxon_figure_7(
  data = data_figure_7_hs,
  outcome = "chitin.ug.mg",
  pathogen_label = "H. solani (spores)",
  outcome_label = "Chitin"
)

# Combine all statistical results
stats_figure_7 <- bind_rows(
  stats_glucan_cc_figure_7,
  stats_chitin_cc_figure_7,
  stats_glucan_hs_figure_7,
  stats_chitin_hs_figure_7
)

# Print results
print(stats_figure_7, n = Inf)

# ------------------------------------------------------------
# 5. Common figure theme
# ------------------------------------------------------------

theme_figure_7 <- theme_classic(base_size = 14) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    
    axis.line = element_line(
      colour = "black",
      linewidth = 0.7
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.7
    ),
    axis.text = element_text(
      colour = "black",
      size = 11
    ),
    axis.title = element_text(
      colour = "black",
      size = 12
    ),
    
    plot.title = element_text(
      face = "italic",
      size = 14,
      hjust = 0.5,
      colour = "black"
    ),
    
    plot.tag = element_text(
      face = "bold",
      size = 16,
      colour = "black"
    ),
    plot.tag.position = c(0, 1),
    plot.title.position = "plot",
    
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )

treatment_labels_figure_7 <- c(
  "PDB" = "PDB",
  "PDB+MeOH" = "PDB+MeOH",
  "Frangula" = expression(italic("Frangula")),
  "Rheum" = expression(italic("Rheum")),
  "Curcuma" = expression(italic("Curcuma"))
)

# ------------------------------------------------------------
# Significance labels for C. coccodes glucan only
# ------------------------------------------------------------


significance_glucan_cc_figure_7 <- stats_glucan_cc_figure_7 |>
  filter(Comparison %in% c(
    "PDB vs Frangula",
    "PDB+MeOH vs Rheum",
    "PDB+MeOH vs Curcuma"
  )) |>
  mutate(
    Treatment = case_when(
      Comparison == "PDB vs Frangula" ~ "Frangula",
      Comparison == "PDB+MeOH vs Rheum" ~ "Rheum",
      Comparison == "PDB+MeOH vs Curcuma" ~ "Curcuma"
    ),
    significant_glc = !is.na(p_adjusted) & p_adjusted < 0.05
  )

summary_figure_7_cc <- summary_figure_7_cc |>
  left_join(
    significance_glucan_cc_figure_7 |>
      select(Treatment, significant_glc),
    by = "Treatment"
  ) |>
  mutate(
    significant_glc = replace_na(significant_glc, FALSE)
  )

# ------------------------------------------------------------
# 6. Function to create the figure
# ------------------------------------------------------------

make_plot_figure_7 <- function(
    summary_data,
    raw_data,
    mean_col,
    sd_col,
    value_col,
    y_label,
    pathogen_title,
    panel_letter,
    show_significance = FALSE
) {
  
  # Significance annotations (only if available)
  if (show_significance &&
      "significant_glc" %in% names(summary_data)) {
    
    significance_data <- summary_data |>
      filter(significant_glc) |>
      mutate(
        star_y = .data[[mean_col]] + .data[[sd_col]] +
          max(summary_data[[sd_col]], na.rm = TRUE) * 0.25
      )
    
  } else {
    significance_data <- summary_data[0, ]
    significance_data$star_y <- numeric(0)
  }
  
  # Determine suitable y-axis range from the raw observations
  max_value <- max(
    c(
      raw_data[[value_col]],
      summary_data[[mean_col]] + summary_data[[sd_col]]
    ),
    na.rm = TRUE
  )
  
  y_upper <- ceiling(max_value * 1.20)
  
  ggplot(
    summary_data,
    aes(
      x = Treatment,
      y = .data[[mean_col]]
    )
  ) +
    
    geom_col(
      fill = "#599AD3",
      width = 0.70,
      na.rm = TRUE
    ) +
    
    geom_errorbar(
      aes(
        ymin = .data[[mean_col]] - .data[[sd_col]],
        ymax = .data[[mean_col]] + .data[[sd_col]]
      ),
      width = 0.20,
      linewidth = 0.8,
      colour = "#4D4D4D",
      na.rm = TRUE
    ) +
    
   
    # Significance stars, only when requested
    {
      if (show_significance) {
        geom_text(
          data = significance_data,
          aes(
            x = Treatment,
            y = star_y,
            label = "*"
          ),
          inherit.aes = FALSE,
          size = 6,
          fontface = "bold",
          na.rm = TRUE
        )
      }
    } +
    
    scale_x_discrete(
      limits = treatment_order_figure_7,
      labels = treatment_labels_figure_7
    ) +
    
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 5),
      expand = expansion(mult = c(0, 0.03))
    ) +
    
    coord_cartesian(
      ylim = c(0, y_upper),
      clip = "off"
    ) +
    
    labs(
      title = pathogen_title,
      tag = panel_letter,
      x = "Treatment",
      y = y_label
    ) +
    
    theme_figure_7
}


# ------------------------------------------------------------
# 7. Create the four panels
# ------------------------------------------------------------

# A — Glucan, C. coccodes
plot_glucan_cc_figure_7 <- make_plot_figure_7(
  summary_data = summary_figure_7_cc,
  raw_data = data_figure_7_cc,
  mean_col = "mean_glucan",
  sd_col = "sd_glucan",
  value_col = "glucan.ug.mg",
  y_label = "Glucan content (µg mg\u207B\u00B9)",
  pathogen_title = "",
  panel_letter = "A",
  show_significance = TRUE
)

# B — Glucan, H. solani
plot_glucan_hs_figure_7 <- make_plot_figure_7(
  summary_data = summary_figure_7_hs,
  raw_data = data_figure_7_hs,
  mean_col = "mean_glucan",
  sd_col = "sd_glucan",
  value_col = "glucan.ug.mg",
  y_label = "Glucan content (µg mg\u207B\u00B9)",
  pathogen_title = "",
  panel_letter = "B"
)

# C — Chitin, C. coccodes
plot_chitin_cc_figure_7 <- make_plot_figure_7(
  summary_data = summary_figure_7_cc,
  raw_data = data_figure_7_cc,
  mean_col = "mean_chitin",
  sd_col = "sd_chitin",
  value_col = "chitin.ug.mg",
  y_label = "Chitin content (µg mg\u207B\u00B9)",
  pathogen_title = "",
  panel_letter = "C"
)

# D — Chitin, H. solani
plot_chitin_hs_figure_7 <- make_plot_figure_7(
  summary_data = summary_figure_7_hs,
  raw_data = data_figure_7_hs,
  mean_col = "mean_chitin",
  sd_col = "sd_chitin",
  value_col = "chitin.ug.mg",
  y_label = "Chitin content (µg mg\u207B\u00B9)",
  pathogen_title = "",
  panel_letter = "D"
)



# ------------------------------------------------------------
# 8. Combine and export Figure 7
# ------------------------------------------------------------

figure_7 <- grid.arrange(
  plot_glucan_cc_figure_7,
  plot_glucan_hs_figure_7,
  plot_chitin_cc_figure_7,
  plot_chitin_hs_figure_7,
  ncol = 2,
  nrow = 2
)

ggsave(
  file.path(output_dir, "figure_7.png"),
  plot = figure_7,
  width = 12,
  height = 9,
  dpi = 400,
  bg = "white"
)


# ============================================================
# FIGURE 8: Post-harvest disease severity
# ============================================================


# ------------------------------------------------------------
# 1. Import data
# ------------------------------------------------------------

# Data source: Zenodo (https://zenodo.org/records/24111700)
zenodo_url_figure_8  <- "https://zenodo.org/records/24111700/files/postharvest_data_analysis.csv?download=1"
zenodo_file_figure_8 <- tempfile(fileext = ".csv")
download.file(zenodo_url_figure_8, zenodo_file_figure_8, mode = "wb", quiet = TRUE)

raw_data_001 <- read.table(
  zenodo_file_figure_8,
  header = TRUE,
  sep = ";",
  dec = ".",
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

# Trim whitespace from text columns
raw_data_001 <- raw_data_001 |>
  mutate(across(where(is.character), stringr::str_trim))


# ------------------------------------------------------------
# 2. Validate the data
# ------------------------------------------------------------

required_columns <- c(
  "Unit.number", "Treatment", "Product", "Year",
  "Cc.0", "Cc.1", "Cc.2", "Cc.3", "Cc.4",
  "Hs.0", "Hs.1", "Hs.2", "Hs.3", "Hs.4"
)

missing_columns <- setdiff(
  required_columns,
  names(raw_data_001)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing columns in CSV: ",
    paste(missing_columns, collapse = ", ")
  )
}

product_order_figure_8 <- c(
  "Water", "Amistar", "Curcuma", "Frangula", "Rheum"
)

unexpected_products <- setdiff(
  unique(na.omit(raw_data_001$Product)),
  product_order_figure_8
)

if (length(unexpected_products) > 0) {
  stop(
    "Unexpected Product labels: ",
    paste(unexpected_products, collapse = ", ")
  )
}

if (anyNA(raw_data_001$Unit.number) ||
    anyNA(raw_data_001$Product) ||
    anyNA(raw_data_001$Year)) {
  stop("Unit.number, Product and Year must not contain missing values.")
}

# Confirm that the counts are numeric
count_columns_figure_8 <- c(
  "Cc.0", "Cc.1", "Cc.2", "Cc.3", "Cc.4",
  "Hs.0", "Hs.1", "Hs.2", "Hs.3", "Hs.4"
)

raw_data_001 <- raw_data_001 |>
  mutate(
    across(
      all_of(count_columns_figure_8),
      as.numeric
    )
  )

if (anyNA(raw_data_001[count_columns_figure_8])) {
  stop("Disease-category counts contain missing or non-numeric values.")
}

if (any(as.matrix(raw_data_001[count_columns_figure_8]) < 0)) {
  stop("Disease-category counts cannot be negative.")
}

# ------------------------------------------------------------
# 3. Pool counts across observers within each unit
# ------------------------------------------------------------

# Each observer assessed different tubers.
# Therefore, sum their category counts within each unit,
# product, treatment and year before calculating severity.

data_units_figure_8 <- raw_data_001 |>
  group_by(Unit.number, Treatment, Product, Year) |>
  summarise(
    across(
      all_of(count_columns_figure_8),
      ~ sum(.x)
    ),
    .groups = "drop"
  ) |>
  mutate(
    Product = factor(
      Product,
      levels = product_order_figure_8
    ),
    
    # Total tubers assessed for each pathogen
    Cc.total = Cc.0 + Cc.1 + Cc.2 + Cc.3 + Cc.4,
    Hs.total = Hs.0 + Hs.1 + Hs.2 + Hs.3 + Hs.4,
    
    # Disease severity (%), using category scores:
    # 0%, 10%, 25%, 50%, 75%
    Cc.Tax = if_else(
      Cc.total > 0,
      (Cc.1 * 10 +
         Cc.2 * 25 +
         Cc.3 * 50 +
         Cc.4 * 75) / Cc.total,
      NA_real_
    ),
    
    Hs.Tax = if_else(
      Hs.total > 0,
      (Hs.1 * 10 +
         Hs.2 * 25 +
         Hs.3 * 50 +
         Hs.4 * 75) / Hs.total,
      NA_real_
    )
  )

# Inspect the resulting unit-level dataset
print(data_units_figure_8)

# Number of independent units per product and year
replication_figure_8 <- data_units_figure_8 |>
  group_by(Product, Year) |>
  summarise(
    n_units = n(),
    .groups = "drop"
  )

print(replication_figure_8)

# ------------------------------------------------------------
# 4. Prepare the two pathogen datasets
# ------------------------------------------------------------

# C. coccodes:
data_figure_8_cc <- data_units_figure_8 |>
  filter(!is.na(Cc.Tax))

# H. solani: retain all available years
data_figure_8_hs <- data_units_figure_8 |>
  filter(!is.na(Hs.Tax))

# ------------------------------------------------------------
# 5. Descriptive statistics at the unit level
# ------------------------------------------------------------

summarise_figure_8 <- function(data, outcome) {
  data |>
    group_by(Product) |>
    summarise(
      n_units = sum(!is.na(.data[[outcome]])),
      mean = mean(.data[[outcome]], na.rm = TRUE),
      sd = sd(.data[[outcome]], na.rm = TRUE),
      median = median(.data[[outcome]], na.rm = TRUE),
      .groups = "drop"
    )
}

summary_figure_8_cc <- summarise_figure_8(
  data_figure_8_cc, "Cc.Tax"
)

summary_figure_8_hs <- summarise_figure_8(
  data_figure_8_hs, "Hs.Tax"
)

print(summary_figure_8_cc)
print(summary_figure_8_hs)


# --------------------------------------------------------------------------
# 6. Linear model with Product and Year as fixed effects and Holm correction
# --------------------------------------------------------------------------

library(emmeans)
library(dplyr)

# Ensure year and product are factors
data_figure_8_cc <- data_figure_8_cc |>
  mutate(
    Product = factor(
      Product,
      levels = c("Water", "Amistar", "Curcuma", "Frangula", "Rheum")
    ),
    Year = factor(Year)
  )

data_figure_8_hs <- data_figure_8_hs |>
  mutate(
    Product = factor(
      Product,
      levels = c("Water", "Amistar", "Curcuma", "Frangula", "Rheum")
    ),
    Year = factor(Year)
  )

# Additive models: product effect adjusted for year
model_figure_8_cc <- lm(
  Cc.Tax ~ Product + Year,
  data = data_figure_8_cc
)

model_figure_8_hs <- lm(
  Hs.Tax ~ Product + Year,
  data = data_figure_8_hs
)

# Estimated product means, adjusted for year
emm_cc <- emmeans::emmeans(model_figure_8_cc, ~ Product)
emm_hs <- emmeans::emmeans(model_figure_8_hs, ~ Product)

# Planned comparisons against Water, with Holm correction
contrasts_vs_water <- list(
  "Amistar vs Water"  = c(-1, 1, 0, 0, 0),
  "Curcuma vs Water"  = c(-1, 0, 1, 0, 0),
  "Frangula vs Water" = c(-1, 0, 0, 1, 0),
  "Rheum vs Water"    = c(-1, 0, 0, 0, 1)
)

stats_cc <- emmeans::contrast(
  emm_cc,
  method = contrasts_vs_water,
  adjust = "holm"
)

stats_hs <- emmeans::contrast(
  emm_hs,
  method = contrasts_vs_water,
  adjust = "holm"
)

summary(stats_cc, infer = TRUE)
summary(stats_hs, infer = TRUE)

#Diagnostics plots:

par(mfrow = c(2, 2))
plot(model_figure_8_cc)

par(mfrow = c(2, 2))
plot(model_figure_8_hs)

par(mfrow = c(1, 1))

# ------------------------------------------------------------
# 7. Plot settings
# ------------------------------------------------------------

theme_figure_8 <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    ),
    plot.title = element_text(
      size = 13,
      hjust = 0.5
    ),
    plot.tag = element_text(
      size = 15,
      face = "bold"
    ),
    plot.tag.position = c(0.02, 0.98),
    plot.margin = margin(10, 12, 10, 10)
  )

shapes_year_figure_8 <- c(
  "2017" = 21,  # circle
  "2018" = 22,  # square
  "2019" = 24   # triangle
)
colours_figure_8 <- c(
  "Water"    = "#FFFFFF",
  "Amistar"  = "#CCCCCC",
  "Curcuma"  = brewer.pal(n = 3, name = "Blues")[1],
  "Frangula" = brewer.pal(n = 3, name = "Blues")[2],
  "Rheum"    = brewer.pal(n = 3, name = "Blues")[3]
)

treatment_labels_figure_8 <- c(
  "Water" = "Water",
  "Amistar" = "Amistar",
  "Frangula" = expression(italic("Frangula")),
  "Rheum" = expression(italic("Rheum")),
  "Curcuma" = expression(italic("Curcuma"))
)

# ------------------------------------------------------------
# 8. Function to create each panel
# ------------------------------------------------------------

make_plot_figure_8 <- function(
    data,
    outcome,
    y_label,
    pathogen_title,
    panel_letter,
    stats = NULL,
    show_significance = FALSE,
    alpha_levels = c("***" = 0.001, "**" = 0.01, "*" = 0.05)) {

  # Significance annotations are optional. Figure 8 as published shows
  # none (show_significance = FALSE). When TRUE, stars are derived only from
  # the Holm-adjusted p-values in `stats` (the emmeans contrasts vs Water);
  # the product of each contrast is the level carrying the +1 coefficient in
  # `contrasts_vs_water`, so no product names or labels are typed here.
  if (show_significance) {

    p_values <- as.data.frame(summary(stats))$p.value

    # Label of the smallest threshold the p-value falls below; "" if none
    stars <- vapply(
      p_values,
      function(p) {
        passed <- names(alpha_levels)[p < alpha_levels]
        if (length(passed) == 0) "" else passed[which.min(alpha_levels[passed])]
      },
      character(1)
    )

    significance_data <- data.frame(
      Product = levels(data$Product)[
        vapply(contrasts_vs_water, which.max, integer(1))
      ],
      significance = stars,
      stringsAsFactors = FALSE
    ) |>
      dplyr::filter(significance != "")

    # Place stars above the highest observed unit value
    star_positions <- data |>
      dplyr::group_by(Product) |>
      dplyr::summarise(
        star_y = max(.data[[outcome]], na.rm = TRUE) + 5,
        .groups = "drop"
      ) |>
      dplyr::mutate(Product = as.character(Product))

    significance_data <- significance_data |>
      dplyr::left_join(star_positions, by = "Product") |>
      dplyr::mutate(
        Product = factor(
          Product,
          levels = product_order_figure_8
        )
      )

  } else {
    significance_data <- data.frame(
      Product = factor(character(0), levels = product_order_figure_8),
      significance = character(0),
      star_y = numeric(0)
    )
  }

  ggplot2::ggplot(
    data,
    ggplot2::aes(x = Product, y = .data[[outcome]])
  ) +
    
    # Boxplots show the distribution within each treatment
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = Product),
      width = 0.62,
      outlier.shape = NA,
      colour = "black",
      linewidth = 0.6
    ) +
    
    ggplot2::scale_fill_manual(
      values = colours_figure_8,
      guide = "none"
    ) +
    
    # Individual independent storage units; shape indicates year
    # A fixed seed makes the horizontal jitter reproducible across runs
    ggplot2::geom_jitter(
      ggplot2::aes(shape = factor(Year)),
      position = ggplot2::position_jitter(
        width = 0.09,
        height = 0,
        seed = 2026
      ),
      size = 2.5,
      fill = "black",
      colour = "black",
      alpha = 0.8
    ) +
    
    ggplot2::scale_shape_manual(
      name = "Year",
      values = shapes_year_figure_8,
      breaks = names(shapes_year_figure_8),
      drop = FALSE
    ) +
    
    # Significance annotations
    ggplot2::geom_text(
      data = significance_data,
      ggplot2::aes(
        x = Product,
        y = star_y,
        label = significance
      ),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold",
      na.rm = TRUE
    ) +
    
    ggplot2::scale_x_discrete(
      limits = product_order_figure_8,
      labels = treatment_labels_figure_8,
      drop = FALSE
    ) +
    
    ggplot2::scale_y_continuous(
      breaks = seq(0, 80, by = 20),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    
    ggplot2::coord_cartesian(
      ylim = c(0, 79),
      clip = "off"
    ) +
    
    ggplot2::labs(
      title = pathogen_title,
      tag = panel_letter,
      x = "Post-harvest treatment",
      y = y_label
    ) +
    
    theme_figure_8
}

# ------------------------------------------------------------
# 9. Create panels A and B
# ------------------------------------------------------------

plot_figure_8_cc <- make_plot_figure_8(
  data = data_figure_8_cc,
  outcome = "Cc.Tax",
  y_label = "Black dot disease severity (%)",
  pathogen_title = expression(italic("")),
  panel_letter = "A",
  stats = stats_cc,
  show_significance = TRUE
)

plot_figure_8_hs <- make_plot_figure_8(
  data = data_figure_8_hs,
  outcome = "Hs.Tax",
  y_label = "Silver scurf disease severity (%)",
  pathogen_title = expression(italic("")),
  panel_letter = "B",
  stats = stats_hs,
  show_significance = TRUE
)

# Display the panels
figure_8 <- gridExtra::arrangeGrob(
  plot_figure_8_cc,
  plot_figure_8_hs,
  ncol = 2
)

grid::grid.newpage()
grid::grid.draw(figure_8)


# ------------------------------------------------------------
# Table S3a: Descriptive statistics by product and pathogen
# Each row in the input data must represent one storage unit
# ------------------------------------------------------------

summarise_figure_8 <- function(data, outcome, pathogen_label) {
  data |>
    group_by(Product) |>
    summarise(
      n_units = sum(!is.na(.data[[outcome]])),
      mean = mean(.data[[outcome]], na.rm = TRUE),
      sd = sd(.data[[outcome]], na.rm = TRUE),
      median = median(.data[[outcome]], na.rm = TRUE),
      q1 = quantile(.data[[outcome]], 0.25, na.rm = TRUE),
      q3 = quantile(.data[[outcome]], 0.75, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      Pathogen = pathogen_label,
      mean_sd = sprintf("%.2f ± %.2f", mean, sd),
      median_iqr = sprintf("%.2f (%.2f–%.2f)", median, q1, q3)
    ) |>
    select(
      Pathogen, Product, n_units,
      mean_sd, median_iqr,
      mean, sd, median, q1, q3
    )
}

table_S3a <- bind_rows(
  summarise_figure_8(
    data_figure_8_cc, "Cc.Tax", "Black dot"
  ),
  summarise_figure_8(
    data_figure_8_hs, "Hs.Tax", "Silver scurf"
  )
) |>
  mutate(
    Product = factor(
      Product,
      levels = product_order_figure_8
    )
  ) |>
  arrange(Pathogen, Product) |>
  mutate(Product = as.character(Product))

# ------------------------------------------------------------
# Table S3b: Year-adjusted comparisons against Water
# ------------------------------------------------------------

extract_model_comparisons <- function(stats, pathogen_label) {
  as.data.frame(summary(stats, infer = TRUE)) |>
    mutate(
      Pathogen = pathogen_label,
      comparison = contrast,
      estimate = estimate,
      SE = SE,
      df = df,
      CI_lower = lower.CL,
      CI_upper = upper.CL,
      p_Holm = p.value
    ) |>
    select(
      Pathogen, comparison, estimate, SE, df,
      CI_lower, CI_upper, p_Holm
    )
}

table_S3b <- bind_rows(
  extract_model_comparisons(stats_cc, "Black dot"),
  extract_model_comparisons(stats_hs, "Silver scurf")
) |>
  mutate(
    across(
      c(estimate, SE, CI_lower, CI_upper),
      ~ round(.x, 2)
    ),
    p_Holm = signif(p_Holm, 3)
  )

# ------------------------------------------------------------
# 10. Export the figure and tables
# ------------------------------------------------------------

ggsave(
  file.path(output_dir, "figure_8.png"),
  plot = figure_8,
  width = 12,
  height = 5.5,
  dpi = 600,
  bg = "white"
)

write.csv(
  table_S3a,
  file.path(output_dir, "Supplementary_Table_S3a_Descriptive.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  table_S3b,
  file.path(output_dir, "Supplementary_Table_S3b_Model_comparisons.csv"),
  row.names = FALSE,
  na = ""
)
