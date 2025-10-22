## Introduction

This document performs spatially explicit pixel-by-pixel model agreement
analysis to assess algorithm-specific uncertainty in presence-only
distribution modeling for *Bactericera cockerelli*. The analysis
compares **MaxEnt** and **Random Forest (RF)** predictions using
**global-scale binarized models (Bco-Glob)** across three climate
scenarios.

The comparison quantifies prediction uncertainty by systematically
evaluating algorithmic differences, dividing forecasting concordance
into three categories:

-   **Stable Area** (blue): Both algorithms consistently predict
    establishment potential
-   **New Area ↑Risk** (red): RF predicts larger suitable area than
    MaxEnt (increased risk under RF)
-   **Area Loss ↓Risk** (green): MaxEnt predicts larger suitable area
    than RF (reduced risk under RF, heightened under MaxEnt)

This multi-algorithm validation approach provides robust uncertainty
bounds for *B. cockerelli*-transmitted disease risk assessments,
supporting integrated pest management and phytosanitary surveillance
strategies.

------------------------------------------------------------------------

## Environment Configuration

Configure parallel processing and memory optimization for comparative
spatial analysis.

    library(future)
    library(furrr)
    library(parallel)

    # Clear memory and enable efficient raster processing
    gc()

    ##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells  5333426 284.9   17670570  943.8  11421366  610.0
    ## Vcells 56087761 428.0  149388046 1139.8 149388046 1139.8

    options(terra.blocks = TRUE)

    # Configure parallel processing with 24 cores
    plan(multisession, workers = 24)

## Load Required Libraries

Import essential packages for spatial analysis, visualization, and
algorithm comparison.

    library(terra)
    library(SDMtune)
    library(readxl)
    library(dplyr)
    library(tidyverse)
    library(raster)
    library(geodata)
    library(predicts)
    library(tidyterra)
    library(cluster)
    library(png)
    library(grid)
    library(rgl)
    library(ggplot2)

## Import Binary Predictions - MaxEnt

Load global-scale binarized MaxEnt predictions (Bco-Glob) using omission
rate thresholds.

    # Load global-scale MaxEnt binary predictions
    A <- readRDS("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Paper_comnplejo_punta_morada/Maxent/Prediction_Bactericera_Global_binarizado.rds")

## Import Binary Predictions - Random Forest

Load global-scale binarized Random Forest predictions (Bco-Glob) using
omission rate thresholds.

    # Load global-scale Random Forest binary predictions
    B <- readRDS("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Paper_comnplejo_punta_morada/Random Forest/Predicion_Bactericera_Global_Binarizada.rds")

## Import Geographic Boundaries

Load Colombian department polygons for spatial context and visualization
masking.

    # Import Colombia administrative boundaries (WGS84 projection)
    colombia_shp <- vect("C:/Users/manuel.cortes/Downloads/Departamentos2010WGS84.shp")

    # Convert to dataframe for ggplot2 compatibility
    colombia_df <- as.data.frame(geom(colombia_shp), geom = "XY")

## Data Structure Verification

Inspect dataset structure and verify scenario integrity for both
algorithms before comparison.

    print("Verifying data structure...")

    ## [1] "Verifying data structure..."

    print("=== OBJECT A (MAXENT) ===")

    ## [1] "=== OBJECT A (MAXENT) ==="

    print(paste("Total records in A$data:", nrow(A$data)))

    ## [1] "Total records in A$data: 163014"

    print("Unique scenarios in A:")

    ## [1] "Unique scenarios in A:"

    print(table(A$data$escenario))

    ## 
    ##     A     B     C 
    ## 54338 54338 54338

    print("Unique values in lyr1 (A):")

    ## [1] "Unique values in lyr1 (A):"

    print(table(A$data$lyr1))

    ## 
    ##  Ausencia Presencia 
    ##    145269     17745

    print("=== OBJECT B (RANDOM FOREST) ===")

    ## [1] "=== OBJECT B (RANDOM FOREST) ==="

    print(paste("Total records in B$data:", nrow(B$data)))

    ## [1] "Total records in B$data: 163014"

    print("Unique scenarios in B:")

    ## [1] "Unique scenarios in B:"

    print(table(B$data$escenario))

    ## 
    ##     A     B     C 
    ## 54338 54338 54338

    print("Unique values in lyr1 (B):")

    ## [1] "Unique values in lyr1 (B):"

    print(table(B$data$lyr1))

    ## 
    ##  Ausencia Presencia 
    ##    127543     35471

## Extract Available Scenarios

Identify available climate scenarios in both MaxEnt and Random Forest
datasets.

    # Extract unique scenario names from both objects
    escenarios_A <- unique(A$data$escenario)
    escenarios_B <- unique(B$data$escenario)

    print("Available scenarios:")

    ## [1] "Available scenarios:"

    print(paste("A (Maxent):", paste(escenarios_A, collapse = ", ")))

    ## [1] "A (Maxent): A, B, C"

    print(paste("B (Random Forest):", paste(escenarios_B, collapse = ", ")))

    ## [1] "B (Random Forest): A, B, C"

## Define Data Extraction Function

Create standardized function to extract and prepare scenario-specific
data for both algorithms.

    # Function to extract data by scenario and convert to binary presence/absence
    extract_scenario_data <- function(data_obj, scenario_name) {
      scenario_data <- data_obj$data[data_obj$data$escenario == scenario_name, ]
      
      # Convert categorical labels to binary (1 = presence, 0 = absence)
      scenario_data$presence <- ifelse(scenario_data$lyr1 == "Presencia", 1, 0)
      
      return(scenario_data[, c("x", "y", "presence")])
    }

## Define Algorithm Comparison Function

Create function to perform pixel-by-pixel comparison between MaxEnt and
Random Forest predictions.

    # Function to compare MaxEnt vs Random Forest predictions
    compare_methods <- function(data_maxent, data_rf, comparison_name) {
      # Create complete coordinate grid from both datasets
      all_coords <- unique(rbind(data_maxent[, c("x", "y")], data_rf[, c("x", "y")]))
      
      # Merge with complete grid to ensure spatial alignment
      merged_maxent <- merge(all_coords, data_maxent, by = c("x", "y"), all.x = TRUE)
      merged_rf <- merge(all_coords, data_rf, by = c("x", "y"), all.x = TRUE)
      
      # Replace NA with 0 (absence where no prediction exists)
      merged_maxent$presence[is.na(merged_maxent$presence)] <- 0
      merged_rf$presence[is.na(merged_rf$presence)] <- 0
      
      # Combine both datasets with suffixes for differentiation
      comparison_data <- merge(merged_maxent, merged_rf, by = c("x", "y"), 
                               suffixes = c("_maxent", "_rf"))
      
      # Classify agreement types between algorithms
      comparison_data$agreement_type <- ifelse(
        comparison_data$presence_maxent == 1 & comparison_data$presence_rf == 1, "Coinciden_Presencia",
        ifelse(comparison_data$presence_maxent == 0 & comparison_data$presence_rf == 0, "Coinciden_Ausencia",
               ifelse(comparison_data$presence_maxent == 0 & comparison_data$presence_rf == 1, "RF_Mayor_Riesgo",
                      ifelse(comparison_data$presence_maxent == 1 & comparison_data$presence_rf == 0, "Maxent_Mayor_Riesgo", "Error"))))
      
      # Add scenario identifier
      comparison_data$escenario <- comparison_name
      
      return(comparison_data)
    }

------------------------------------------------------------------------

## Scenario 1: Current Climate Comparison

Compare MaxEnt and Random Forest predictions under current climate
conditions (first scenario position).

    print("Comparing algorithms by scenario position...")

    ## [1] "Comparing algorithms by scenario position..."

    # Initialize empty dataframe for all comparisons
    all_comparisons <- data.frame()

    # FIRST SCENARIO (CURRENT CLIMATE)
    if (length(escenarios_A) >= 1 & length(escenarios_B) >= 1) {
      data_maxent_1 <- extract_scenario_data(A, escenarios_A[1])
      data_rf_1 <- extract_scenario_data(B, escenarios_B[1])
      
      comparison_1 <- compare_methods(data_maxent_1, data_rf_1, "a")
      all_comparisons <- rbind(all_comparisons, comparison_1)
    }

------------------------------------------------------------------------

## Scenario 2: SSP2-4.5 Comparison

Compare MaxEnt and Random Forest predictions under moderate mitigation
scenario (SSP2-4.5, 2021-2040).

    # SECOND SCENARIO (SSP2-4.5)
    if (length(escenarios_A) >= 2 & length(escenarios_B) >= 2) {
      data_maxent_2 <- extract_scenario_data(A, escenarios_A[2])
      data_rf_2 <- extract_scenario_data(B, escenarios_B[2])
      
      comparison_2 <- compare_methods(data_maxent_2, data_rf_2, "b")
      all_comparisons <- rbind(all_comparisons, comparison_2)
    }

------------------------------------------------------------------------

## Scenario 3: SSP5-8.5 Comparison

Compare MaxEnt and Random Forest predictions under high emissions
scenario (SSP5-8.5, 2021-2040).

    # THIRD SCENARIO (SSP5-8.5)
    if (length(escenarios_A) >= 3 & length(escenarios_B) >= 3) {
      data_maxent_3 <- extract_scenario_data(A, escenarios_A[3])
      data_rf_3 <- extract_scenario_data(B, escenarios_B[3])
      
      comparison_3 <- compare_methods(data_maxent_3, data_rf_3, "c")
      all_comparisons <- rbind(all_comparisons, comparison_3)
    }

------------------------------------------------------------------------

## Prepare Data for Visualization

Consolidate all algorithm comparisons and filter relevant disagreements
for mapping.

    print("Preparing data for visualization...")

    ## [1] "Preparing data for visualization..."

    print(paste("Total comparisons:", nrow(all_comparisons)))

    ## [1] "Total comparisons: 163014"

    # Display summary of agreement types by scenario
    print("Summary of algorithm agreement:")

    ## [1] "Summary of algorithm agreement:"

    print(table(all_comparisons$escenario, all_comparisons$agreement_type))

    ##    
    ##     Coinciden_Ausencia Coinciden_Presencia Maxent_Mayor_Riesgo RF_Mayor_Riesgo
    ##   a              42562                5579                 112            6085
    ##   b              42370                5881                  62            6025
    ##   c              42378                6052                  59            5849

    # Filter only relevant pixels (exclude mutual absence)
    relevant_comparisons <- all_comparisons[all_comparisons$agreement_type != "Coinciden_Ausencia", ]

    print(paste("Relevant pixels for visualization:", nrow(relevant_comparisons)))

    ## [1] "Relevant pixels for visualization: 35704"

------------------------------------------------------------------------

## Generate Algorithm Comparison Maps

Create faceted visualization showing algorithm-specific uncertainty
across three climate scenarios.

    if(nrow(relevant_comparisons) > 0) {
      
      # Convert agreement types to categorical labels with forced factor levels
      relevant_comparisons$lyr1 <- factor(relevant_comparisons$agreement_type,
                                          levels = c("Coinciden_Presencia", "RF_Mayor_Riesgo", "Maxent_Mayor_Riesgo"),
                                          labels = c("Stable Area", 
                                                     "New Area (↑ Risk)", 
                                                     "Area Loss (↓ Risk)"))
      
      # Force all levels to be present (even if no data exists)
      relevant_comparisons$lyr1 <- factor(relevant_comparisons$lyr1, 
                                          levels = c("New Area (↑ Risk)", 
                                                     "Area Loss (↓ Risk)", 
                                                     "Stable Area"))
      
      # Define algorithm comparison color palette
      method_colors <- c(
        "New Area (↑ Risk)" = "red",      # RF predicts higher risk than MaxEnt
        "Stable Area" = "blue",            # Both algorithms agree
        "Area Loss (↓ Risk)" = "green"    # MaxEnt predicts higher risk than RF
      )
      
      # Create faceted map with Colombian geographic boundaries
      comparison_map <- ggplot() +
        # Background layer: Colombian territory
        geom_polygon(data = colombia_df, aes(x = x, y = y, group = geom), 
                     fill = "#449c8c") +
        # Agreement layer: pixel-level algorithm comparison
        geom_tile(data = relevant_comparisons, aes(x = x, y = y, fill = lyr1)) +
        scale_fill_manual(values = method_colors, name = "", 
                          limits = c("New Area (↑ Risk)", "Area Loss (↓ Risk)", "Stable Area"),
                          drop = FALSE,
                          guide = guide_legend(override.aes = list(color = c("red", "green", "blue"),
                                                                   fill = c("red", "green", "blue")))) +
        # Boundary overlay (transparent)
        geom_polygon(data = colombia_df, aes(x = x, y = y), 
                     fill = NA, color = NA) +
        coord_equal() +
        facet_wrap(~escenario, nrow = 1) +
        theme_minimal() +
        theme(
          panel.background = element_rect(fill = "grey90", color = NA),
          panel.grid.major = element_line(color = "white", linewidth = 0.5),
          panel.grid.minor = element_blank(),
          plot.title = element_blank(),
          legend.position = "right",
          legend.title = element_text(face = "bold"),
          strip.text = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 9)
        )
      
      # Display final map
      print(comparison_map)
      
      # Print interpretation guide
      cat("\n=== ANALYSIS COMPLETED SUCCESSFULLY ===\n")
      cat("🔴 RED: New Area (RF predicts increased risk vs MaxEnt)\n")
      cat("🟢 GREEN: Area Loss (MaxEnt predicts increased risk vs RF)\n") 
      cat("🔵 BLUE: Stable Area (both algorithms agree on suitability)\n")
      
    } else {
      cat("ERROR: No significant differences found between algorithms.\n")
      cat("Possible causes:\n")
      cat("1. Algorithms produce very similar predictions\n")
      cat("2. Data structure inconsistency\n")
      cat("3. Universal absence in both algorithm predictions\n")
    }

![](C:\Users\manuel.cortes\Desktop\PMP2024\Bactericera\Paper_comnplejo_punta_morada\MaxEnt%20VS%20RF_files/figure-markdown_strict/generate_maps-1.png)

    ## 
    ## === ANALYSIS COMPLETED SUCCESSFULLY ===
    ## 🔴 RED: New Area (RF predicts increased risk vs MaxEnt)
    ## 🟢 GREEN: Area Loss (MaxEnt predicts increased risk vs RF)
    ## 🔵 BLUE: Stable Area (both algorithms agree on suitability)

------------------------------------------------------------------------

## Display Final Map

Render the comparative visualization showing algorithmic concordance and
disagreement.

    plot(comparison_map)

![](C:\Users\manuel.cortes\Desktop\PMP2024\Bactericera\Paper_comnplejo_punta_morada\MaxEnt%20VS%20RF_files/figure-markdown_strict/display_map-1.png)

------------------------------------------------------------------------

## Summary Statistics

Quantify the extent of algorithm agreement and disagreement across
climate scenarios.

    if(nrow(relevant_comparisons) > 0) {
      # Calculate pixel counts by scenario and agreement type
      summary_table <- relevant_comparisons %>%
        group_by(escenario, lyr1) %>%
        summarise(pixel_count = n(), .groups = 'drop') %>%
        pivot_wider(names_from = lyr1, values_from = pixel_count, values_fill = 0)
      
      print(summary_table)
      
      # Calculate percentage agreement
      total_pixels <- relevant_comparisons %>%
        group_by(escenario) %>%
        summarise(total = n(), .groups = 'drop')
      
      agreement_percentage <- relevant_comparisons %>%
        filter(lyr1 == "Stable Area") %>%
        group_by(escenario) %>%
        summarise(stable_pixels = n(), .groups = 'drop') %>%
        left_join(total_pixels, by = "escenario") %>%
        mutate(agreement_percent = (stable_pixels / total) * 100)
      
      print("\nPercentage of algorithm agreement by scenario:")
      print(agreement_percentage)
    }

    ## # A tibble: 3 × 4
    ##   escenario `New Area (↑ Risk)` `Area Loss (↓ Risk)` `Stable Area`
    ##   <chr>                   <int>                <int>         <int>
    ## 1 a                        6085                  112          5579
    ## 2 b                        6025                   62          5881
    ## 3 c                        5849                   59          6052
    ## [1] "\nPercentage of algorithm agreement by scenario:"
    ## # A tibble: 3 × 4
    ##   escenario stable_pixels total agreement_percent
    ##   <chr>             <int> <int>             <dbl>
    ## 1 a                  5579 11776              47.4
    ## 2 b                  5881 11968              49.1
    ## 3 c                  6052 11960              50.6

------------------------------------------------------------------------

## Conclusions

This algorithm comparison analysis reveals **spatially explicit patterns
of prediction uncertainty** in *B. cockerelli* distribution modeling at
the global scale (Bco-Glob). The results demonstrate:

### Key Findings:

1.  **High-Confidence Zones** (blue pixels): Areas where both MaxEnt and
    RF consistently predict establishment potential represent
    high-priority surveillance targets with robust algorithmic support.

2.  **RF-Specific Expansions** (red pixels): Regions where Random Forest
    predicts broader suitable habitat than MaxEnt indicate **potential
    colonization fronts** requiring proactive monitoring, particularly
    in:

    -   Pacific coastal regions (Chocó, Valle del Cauca)
    -   Caribbean lowland departments
    -   Eastern piedmont zones (Meta, Casanare)

3.  **MaxEnt-Specific Predictions** (green pixels): Areas where MaxEnt
    forecasts higher risk than RF suggest **conservative establishment
    zones** where field validation is critical to resolve algorithmic
    uncertainty.

### Methodological Implications:

The differential responses between MaxEnt (conservative, entropy-based)
and RF (ensemble-based, broader gradients) reflect **fundamental
algorithmic differences** in species-environment modeling:

-   **MaxEnt**: Focuses on environmental combinations closest to known
    presence conditions → more conservative predictions
-   **Random Forest**: Identifies broader environmental gradients
    through ensemble learning → potential expansion forecasts

### Management Recommendations:

Areas where **both algorithms converge** (stable zones) warrant
**focused surveillance and management interventions**. Regions of
**model disagreement** (red/green pixels) flag **emerging risk areas
requiring proactive field validation** to resolve prediction uncertainty
and optimize resource allocation for phytosanitary surveillance.

These findings align with the manuscript’s conclusion that
algorithm-specific uncertainty provides valuable insights for risk
stratification, supporting **robust spatial intelligence** for
integrated pest management under climate change scenarios.
