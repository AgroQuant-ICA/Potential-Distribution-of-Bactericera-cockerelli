## Introduction

This document performs spatially explicit pixel-by-pixel comparative
analysis to quantify distributional stability and epidemiological risk
transitions for *Bactericera cockerelli* under climate change scenarios.
Three pairwise comparisons are analyzed: Current vs SSP2-4.5, Current vs
SSP5-8.5, and SSP2-4.5 vs SSP5-8.5 (2021-2040).

The analysis identifies three critical risk categories:

-   **Stable Areas** (blue): Regions with consistent establishment
    potential across scenarios
-   **New Areas ↑Risk** (red): Zones with increased establishment
    potential under future climate
-   **Area Loss ↓Risk** (green): Zones with reduced suitability in
    projected scenarios

------------------------------------------------------------------------

## Environment Configuration

Configure parallel processing and memory management for efficient
spatial analysis.

    library(future)
    library(furrr)
    library(parallel)

    # Clear memory and optimize raster processing
    gc()

    ##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells  5379629 287.4   14142814  755.4  11166829  596.4
    ## Vcells 55018268 419.8  136311769 1040.0 136311769 1040.0

    options(terra.blocks = TRUE)

    # Configure parallel processing with 24 cores
    plan(multisession, workers = 24)

## Load Required Libraries

Import essential packages for spatial analysis and visualization.

    library(terra)
    library(SDMtune)
    library(readxl)
    library(dplyr)
    library(tidyverse)
    library(raster)
    library(geodata)
    library(predicts)
    library(cluster)
    library(png)
    library(ggplot2)

## Import Binary Prediction Data

Load binarized habitat suitability projections derived from MaxEnt
models using omission rate thresholds.

    # Load binary predictions for B. cockerelli across climate scenarios
    A <- readRDS("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Paper_comnplejo_punta_morada/Maxent/Prediction_Bactericera_Colombia_binarizado.rds")
    print(A)

![](C:\Users\manuel.cortes\Desktop\PMP2024\Bactericera\Paper_comnplejo_punta_morada\Spatial%20Risk%20Bactericera%20climate%20change_files/figure-markdown_strict/load_data-1.png)

## Import Geographic Boundaries

Load Colombian department polygons for spatial masking and
visualization.

    # Import Colombia administrative boundaries (WGS84 projection)
    colombia_shp <- vect("C:/Users/manuel.cortes/Downloads/Departamentos2010WGS84.shp")

    # Convert to dataframe for ggplot2 compatibility
    colombia_df <- as.data.frame(geom(colombia_shp), geom = "XY")

## Data Structure Verification

Inspect dataset structure and verify scenario integrity before analysis.

    # Verify data structure and completeness
    print("Data structure verification:")

    ## [1] "Data structure verification:"

    print(paste("Total records:", nrow(A$data)))

    ## [1] "Total records: 163014"

    print("\nScenario distribution:")

    ## [1] "\nScenario distribution:"

    print(table(A$data$escenario))

    ## 
    ##     A     B     C 
    ## 54338 54338 54338

    print("\nPresence/Absence distribution:")

    ## [1] "\nPresence/Absence distribution:"

    print(table(A$data$lyr1))

    ## 
    ##  Ausencia Presencia 
    ##    153386      9628

## Extract Scenario-Specific Data

Partition dataset into three climate scenarios for pairwise comparison.

    # Extract data by climate scenario
    data_A <- A$data[A$data$escenario == "A", ]  # Current climate
    data_B <- A$data[A$data$escenario == "B", ]  # SSP2-4.5 (2021-2040)
    data_C <- A$data[A$data$escenario == "C", ]  # SSP5-8.5 (2021-2040)

    print(paste("Current climate (A):", nrow(data_A), "pixels"))

    ## [1] "Current climate (A): 54338 pixels"

    print(paste("SSP2-4.5 (B):", nrow(data_B), "pixels"))

    ## [1] "SSP2-4.5 (B): 54338 pixels"

    print(paste("SSP5-8.5 (C):", nrow(data_C), "pixels"))

    ## [1] "SSP5-8.5 (C): 54338 pixels"

## Convert Categorical to Binary

Transform categorical presence labels into numeric binary format (0/1)
for computational analysis.

    # Convert categorical labels to binary (1 = presence, 0 = absence)
    data_A$presence <- ifelse(data_A$lyr1 == "Presencia", 1, 0)
    data_B$presence <- ifelse(data_B$lyr1 == "Presencia", 1, 0)
    data_C$presence <- ifelse(data_C$lyr1 == "Presencia", 1, 0)

## Define Comparison Function

Create standardized function to align spatial grids and perform
pixel-by-pixel comparisons between scenarios.

    # Function to prepare datasets for pixel-by-pixel comparison
    prepare_comparison <- function(data1, data2, name1, name2) {
      # Create complete coordinate grid from both datasets
      all_coords <- unique(rbind(data1[, c("x", "y")], data2[, c("x", "y")]))
      
      # Merge with complete grid to ensure spatial alignment
      merged_1 <- merge(all_coords, data1[, c("x", "y", "presence")], 
                        by = c("x", "y"), all.x = TRUE)
      merged_2 <- merge(all_coords, data2[, c("x", "y", "presence")], 
                        by = c("x", "y"), all.x = TRUE)
      
      # Replace NA with 0 (absence where no prediction exists)
      merged_1$presence[is.na(merged_1$presence)] <- 0
      merged_2$presence[is.na(merged_2$presence)] <- 0
      
      # Combine both datasets with suffixes for differentiation
      final_data <- merge(merged_1, merged_2, by = c("x", "y"), 
                          suffixes = c(paste0("_", name1), paste0("_", name2)))
      
      return(final_data)
    }

------------------------------------------------------------------------

## Analysis 1: Current vs SSP2-4.5

Quantify distributional changes between current climate and moderate
mitigation scenario (SSP2-4.5).

    print("Analyzing Current vs SSP2-4.5...")

    ## [1] "Analyzing Current vs SSP2-4.5..."

    coords_AB <- prepare_comparison(data_A, data_B, "A", "B")

    # Classify pixel-level change types
    coords_AB$change_type <- ifelse(coords_AB$presence_A == 0 & coords_AB$presence_B == 0, "Sin_Presencia",
                                    ifelse(coords_AB$presence_A == 1 & coords_AB$presence_B == 0, "Perdida",
                                           ifelse(coords_AB$presence_A == 0 & coords_AB$presence_B == 1, "Nueva",
                                                  ifelse(coords_AB$presence_A == 1 & coords_AB$presence_B == 1, "Estable", "Error"))))

    # Filter relevant changes (exclude non-presence pixels)
    risk_AB <- coords_AB[coords_AB$change_type != "Sin_Presencia", ]
    print(paste("Pixels with changes (A vs B):", nrow(risk_AB)))

    ## [1] "Pixels with changes (A vs B): 3409"

    # Assign categorical labels for visualization
    if(nrow(risk_AB) > 0) {
      risk_AB$lyr1 <- factor(risk_AB$change_type,
                             levels = c("Nueva", "Perdida", "Estable"),
                             labels = c("New Area (↑ Risk)", 
                                        "Area Loss (↓ Risk)", 
                                        "Stable Area"))
      risk_AB$escenario <- "a"
    }

------------------------------------------------------------------------

## Analysis 2: Current vs SSP5-8.5

Quantify distributional changes between current climate and high
emissions scenario (SSP5-8.5).

    print("Analyzing Current vs SSP5-8.5...")

    ## [1] "Analyzing Current vs SSP5-8.5..."

    coords_AC <- prepare_comparison(data_A, data_C, "A", "C")

    # Classify pixel-level change types
    coords_AC$change_type <- ifelse(coords_AC$presence_A == 0 & coords_AC$presence_C == 0, "Sin_Presencia",
                                    ifelse(coords_AC$presence_A == 1 & coords_AC$presence_C == 0, "Perdida",
                                           ifelse(coords_AC$presence_A == 0 & coords_AC$presence_C == 1, "Nueva",
                                                  ifelse(coords_AC$presence_A == 1 & coords_AC$presence_C == 1, "Estable", "Error"))))

    # Filter relevant changes
    risk_AC <- coords_AC[coords_AC$change_type != "Sin_Presencia", ]
    print(paste("Pixels with changes (A vs C):", nrow(risk_AC)))

    ## [1] "Pixels with changes (A vs C): 3485"

    # Assign categorical labels
    if(nrow(risk_AC) > 0) {
      risk_AC$lyr1 <- factor(risk_AC$change_type,
                             levels = c("Nueva", "Perdida", "Estable"),
                             labels = c("New Area (↑ Risk)", 
                                        "Area Loss (↓ Risk)", 
                                        "Stable Area"))
      risk_AC$escenario <- "b"
    }

------------------------------------------------------------------------

## Analysis 3: SSP2-4.5 vs SSP5-8.5

Compare distributional differences between moderate and high emissions
scenarios.

    print("Analyzing SSP2-4.5 vs SSP5-8.5...")

    ## [1] "Analyzing SSP2-4.5 vs SSP5-8.5..."

    coords_BC <- prepare_comparison(data_B, data_C, "B", "C")

    # Classify pixel-level change types
    coords_BC$change_type <- ifelse(coords_BC$presence_B == 0 & coords_BC$presence_C == 0, "Sin_Presencia",
                                    ifelse(coords_BC$presence_B == 1 & coords_BC$presence_C == 0, "Perdida",
                                           ifelse(coords_BC$presence_B == 0 & coords_BC$presence_C == 1, "Nueva",
                                                  ifelse(coords_BC$presence_B == 1 & coords_BC$presence_C == 1, "Estable", "Error"))))

    # Filter relevant changes
    risk_BC <- coords_BC[coords_BC$change_type != "Sin_Presencia", ]
    print(paste("Pixels with changes (B vs C):", nrow(risk_BC)))

    ## [1] "Pixels with changes (B vs C): 3567"

    # Assign categorical labels
    if(nrow(risk_BC) > 0) {
      risk_BC$lyr1 <- factor(risk_BC$change_type,
                             levels = c("Nueva", "Perdida", "Estable"),
                             labels = c("New Area (↑ Risk)", 
                                        "Area Loss (↓ Risk)", 
                                        "Stable Area"))
      risk_BC$escenario <- "c"
    }

------------------------------------------------------------------------

## Consolidate All Analyses

Combine all three pairwise comparisons into unified dataset for
visualization.

    # Initialize empty dataframe
    all_risk_data <- data.frame()

    # Append each analysis if data exists
    if(nrow(risk_AB) > 0) {
      all_risk_data <- rbind(all_risk_data, risk_AB[, c("x", "y", "lyr1", "escenario")])
    }
    if(nrow(risk_AC) > 0) {
      all_risk_data <- rbind(all_risk_data, risk_AC[, c("x", "y", "lyr1", "escenario")])
    }
    if(nrow(risk_BC) > 0) {
      all_risk_data <- rbind(all_risk_data, risk_BC[, c("x", "y", "lyr1", "escenario")])
    }

    print(paste("Total pixels in final analysis:", nrow(all_risk_data)))

    ## [1] "Total pixels in final analysis: 10461"

------------------------------------------------------------------------

## Generate Risk Maps

Create faceted visualization showing climate-driven distributional
shifts across three pairwise scenarios.

    if(nrow(all_risk_data) > 0) {
      
      # Define risk category color palette
      risk_colors <- c(
        "New Area (↑ Risk)" = "red",      # Increased establishment risk
        "Stable Area" = "blue",            # Persistent suitability
        "Area Loss (↓ Risk)" = "green"    # Reduced establishment risk
      )
      
      # Create faceted map with Colombian geographic boundaries
      A_risk <- ggplot() +
        # Background layer: Colombian territory
        geom_polygon(data = colombia_df, aes(x = x, y = y, group = geom), 
                     fill = "#449c8c") +
        # Risk layer: pixel-level change classification
        geom_tile(data = all_risk_data, aes(x = x, y = y, fill = lyr1)) +
        scale_fill_manual(values = risk_colors, name = "") +
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
      print(A_risk)
      
      # Print interpretation guide
      cat("\n=== ANALYSIS COMPLETED SUCCESSFULLY ===\n")
      cat("🔴 RED: New Area (increased establishment risk)\n")
      cat("🟢 GREEN: Area Loss (decreased establishment risk)\n") 
      cat("🔵 BLUE: Stable Area (consistent suitability)\n")
      
    } else {
      cat("ERROR: No significant changes detected between scenarios.\n")
      cat("Possible causes:\n")
      cat("1. Scenarios exhibit minimal divergence\n")
      cat("2. Data structure inconsistency\n")
      cat("3. Universal absence across all scenarios\n")
    }

![](C:\Users\manuel.cortes\Desktop\PMP2024\Bactericera\Paper_comnplejo_punta_morada\Spatial%20Risk%20Bactericera%20climate%20change_files/figure-markdown_strict/generate_maps-1.png)

    ## 
    ## === ANALYSIS COMPLETED SUCCESSFULLY ===
    ## 🔴 RED: New Area (increased establishment risk)
    ## 🟢 GREEN: Area Loss (decreased establishment risk)
    ## 🔵 BLUE: Stable Area (consistent suitability)

------------------------------------------------------------------------

## Summary Statistics

Quantify the extent of each risk category across pairwise comparisons.

    if(nrow(all_risk_data) > 0) {
      summary_table <- all_risk_data %>%
        group_by(escenario, lyr1) %>%
        summarise(pixel_count = n(), .groups = 'drop') %>%
        pivot_wider(names_from = lyr1, values_from = pixel_count, values_fill = 0)
      
      print(summary_table)
    }

    ## # A tibble: 3 × 4
    ##   escenario `New Area (↑ Risk)` `Area Loss (↓ Risk)` `Stable Area`
    ##   <chr>                   <int>                <int>         <int>
    ## 1 a                         541                   36          2832
    ## 2 b                         617                   98          2770
    ## 3 c                         194                  180          3193

------------------------------------------------------------------------

## Conclusions

This pixel-by-pixel analysis provides spatially explicit risk
stratification for *B. cockerelli* under climate change scenarios
through 2040. The results support targeted surveillance strategies by
identifying:

1.  **High-confidence establishment zones** (stable areas across
    scenarios) requiring intensive monitoring
2.  **Emerging risk areas** (new suitable habitat under future climate)
    warranting proactive surveillance
3.  **Retreating populations** (habitat loss zones) indicating reduced
    establishment potential

These findings align with the manuscript’s conclusion that *B.
cockerelli* distribution remains largely constrained to the Eastern
Cordillera due to thermal limitations, with algorithmic uncertainty
(MaxEnt stability vs RF expansions) highlighting the need for field
validation in predicted expansion zones.
