## Introduction

This document implements the Multivariate Environmental Similarity
Surfaces (MESS) analysis for *Bactericera cockerelli* species
distribution modeling in Colombia. Following the manuscript methodology,
we use potato crop suitability maps from Colombia’s UPRA as reference
frameworks for identifying climatically analogous zones.

------------------------------------------------------------------------

## Section 1: Environment Setup and Library Loading

    library(future)
    library(furrr)
    library(parallel)
    gc()

    ##            used  (Mb) gc trigger (Mb) max used (Mb)
    ## Ncells  6888237 367.9   13836597  739 13836597  739
    ## Vcells 16905201 129.0   62385287  476 77981608  595

    options(terra.blocks = TRUE)
    plan(multisession, workers = 24) 

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
    library(grid)
    library(rgl)
    library(dismo)
    library(ggplot2)
    library(tidyterra)

------------------------------------------------------------------------

## Section 2: Load Climate Data for Current and Future Scenarios

    # Load current climate data
    bioclim_data <- rast("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/data/2.5 minutos/bioclim_data.tif")
    Envirem <- rast("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/ensambles mensuales/envirem_vars_historico.tif")
     
    names(bioclim_data) <- paste0("BIO", 1:nlyr(bioclim_data))

    climate_data <- c(bioclim_data, Envirem)
    names(climate_data)

    ##  [1] "BIO1"                     "BIO2"                     "BIO3"                    
    ##  [4] "BIO4"                     "BIO5"                     "BIO6"                    
    ##  [7] "BIO7"                     "BIO8"                     "BIO9"                    
    ## [10] "BIO10"                    "BIO11"                    "BIO12"                   
    ## [13] "BIO13"                    "BIO14"                    "BIO15"                   
    ## [16] "BIO16"                    "BIO17"                    "BIO18"                   
    ## [19] "BIO19"                    "annualPET"                "aridityIndexThornthwaite"
    ## [22] "climaticMoistureIndex"    "continentality"           "embergerQ"               
    ## [25] "growingDegDays0"          "growingDegDays5"          "maxTempColdest"          
    ## [28] "minTempWarmest"           "meanTempColdest"          "meanTempWarmest"         
    ## [31] "monthCountByTemp10"       "PETColdestQuarter"        "PETDriestQuarter"        
    ## [34] "PETseasonality"           "PETWarmestQuarter"        "PETWettestQuarter"       
    ## [37] "thermicityIndex"

    # Load SSP2-4.5 scenario
    ssp245 <- rast("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Ensamble GCMs/245/BIO_ensemble_245_BIO.tif")
    env_245 <- rast("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Ensamble GCMs/245/ENVIREM_FINAL_245.tif")

    names(ssp245) <- paste0("BIO", 1:nlyr(ssp245))
    climate_data245 <- c(ssp245, env_245)
    names(climate_data245)

    ##  [1] "BIO1"                     "BIO2"                     "BIO3"                    
    ##  [4] "BIO4"                     "BIO5"                     "BIO6"                    
    ##  [7] "BIO7"                     "BIO8"                     "BIO9"                    
    ## [10] "BIO10"                    "BIO11"                    "BIO12"                   
    ## [13] "BIO13"                    "BIO14"                    "BIO15"                   
    ## [16] "BIO16"                    "BIO17"                    "BIO18"                   
    ## [19] "BIO19"                    "annualPET"                "aridityIndexThornthwaite"
    ## [22] "climaticMoistureIndex"    "continentality"           "embergerQ"               
    ## [25] "growingDegDays0"          "growingDegDays5"          "maxTempColdest"          
    ## [28] "minTempWarmest"           "meanTempColdest"          "meanTempWarmest"         
    ## [31] "monthCountByTemp10"       "PETColdestQuarter"        "PETDriestQuarter"        
    ## [34] "PETseasonality"           "PETWarmestQuarter"        "PETWettestQuarter"       
    ## [37] "thermicityIndex"

    # Load SSP5-8.5 scenario
    ssp585 <- rast("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Ensamble GCMs/BIO_ensemble_585_BIO.tif")
    env_585 <- rast("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Ensamble GCMs/ENVIREM_FINAL_585.tif")

    names(ssp585) <- paste0("BIO", 1:nlyr(ssp585))
    climate_data585 <- c(ssp585, env_585)
    names(climate_data585)

    ##  [1] "BIO1"                     "BIO2"                     "BIO3"                    
    ##  [4] "BIO4"                     "BIO5"                     "BIO6"                    
    ##  [7] "BIO7"                     "BIO8"                     "BIO9"                    
    ## [10] "BIO10"                    "BIO11"                    "BIO12"                   
    ## [13] "BIO13"                    "BIO14"                    "BIO15"                   
    ## [16] "BIO16"                    "BIO17"                    "BIO18"                   
    ## [19] "BIO19"                    "annualPET"                "aridityIndexThornthwaite"
    ## [22] "climaticMoistureIndex"    "continentality"           "embergerQ"               
    ## [25] "growingDegDays0"          "growingDegDays5"          "maxTempColdest"          
    ## [28] "minTempWarmest"           "meanTempColdest"          "meanTempWarmest"         
    ## [31] "monthCountByTemp10"       "PETColdestQuarter"        "PETDriestQuarter"        
    ## [34] "PETseasonality"           "PETWarmestQuarter"        "PETWettestQuarter"       
    ## [37] "thermicityIndex"

------------------------------------------------------------------------

## Section 3: Variable Selection and Scaling

    # MaxEnt variable selection
    variables_interes <- c(2, 5, 6, 12, 21, 22, 26, 34)

    climate_data_subset <- subset(climate_data, variables_interes)
    climate_data_subset245 <- subset(climate_data245, variables_interes)
    climate_data_subset585 <- subset(climate_data585, variables_interes)

    climate_data <- climate_data_subset
    climate_data245 <- climate_data_subset245
    climate_data585 <- climate_data_subset585

    # Scale predictors
    climate_data <- scale(climate_data)  
    climate_data245 <- scale(climate_data245)
    climate_data585 <- scale(climate_data585)

------------------------------------------------------------------------

## Section 4: Load Presence Records

    obs_data <- read_excel("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Bactericera cockerelli Colombia Oficial surveillance ICA.xlsx")

    obs_data <- obs_data[!is.na(obs_data$longitude) & !is.na(obs_data$latitude), ]
    obs_data <- obs_data[!duplicated(obs_data[, c("longitude", "latitude")]), ]

    presence <- obs_data[, c("longitude", "latitude")]
    presence$pa <- 1
    presence <- presence %>% distinct(longitude, latitude, .keep_all = TRUE)

    presence$longitude <- as.numeric(presence$longitude)
    presence$latitude <- as.numeric(presence$latitude)

------------------------------------------------------------------------

## Section 5: Generate Background Points through Stratified Environmental Sampling

    set.seed(20210707)

    presencias_area <- vect(presence, geom = c("longitude", "latitude"), crs = "epsg:4326")

    tryCatch({
      colombia_boundary <- gadm(country = "COL", level = 0, path = tempdir())
      climate_crop <- crop(climate_data, colombia_boundary)
      climate_crop <- mask(climate_crop, colombia_boundary)
      
      if(all(is.na(values(climate_crop[[1]])))) {
        stop("Error: El recorte de Colombia resultó en valores NA")
      }
      
    }, error = function(e) {
      cat("Error al descargar límites de Colombia:", e$message, "\n")
      cat("Intentando con método alternativo...\n")
      
      colombia_ext <- ext(-82, -66, -5, 13)
      colombia_boundary <- as.polygons(colombia_ext, crs = "epsg:4326")
      climate_crop <- crop(climate_data, colombia_boundary)
      climate_crop <- mask(climate_crop, colombia_boundary)
    })

    tryCatch({
      climate_values <- values(climate_crop, na.rm = TRUE)
      climate_coords <- xyFromCell(climate_crop, which(!is.na(values(climate_crop[[1]]))))
      
      if(nrow(climate_coords) == 0) {
        stop("Error: No se encontraron coordenadas válidas")
      }
      
      cat("Píxeles válidos encontrados:", nrow(climate_coords), "\n")
      
    }, error = function(e) {
      stop("Error al extraer valores climáticos: ", e$message)
    })

    ## Píxeles válidos encontrados: 54539

    complete_cases <- complete.cases(climate_values)
    climate_values_clean <- climate_values[complete_cases, ]
    climate_coords_clean <- climate_coords[complete_cases, ]

    cat("Píxeles después de limpiar NA:", nrow(climate_coords_clean), "\n")

    ## Píxeles después de limpiar NA: 54539

    if(nrow(climate_values_clean) < nrow(presence)) {
      stop("Error: No hay suficientes píxeles válidos para el muestreo")
    }

    tryCatch({
      climate_pca <- prcomp(climate_values_clean, scale. = TRUE, center = TRUE)
      cat("PCA completado. Varianza explicada PC1-3:", 
          round(sum(climate_pca$sdev[1:3]^2)/sum(climate_pca$sdev^2)*100, 2), "%\n")
    }, error = function(e) {
      stop("Error en PCA: ", e$message)
    })

    ## PCA completado. Varianza explicada PC1-3: 88.04 %

    pca_scores <- climate_pca$x[, 1:3]

    n_clusters <- min(50, nrow(presence) * 2, nrow(climate_values_clean) %/% 10)

    cat("Número de clusters a usar:", n_clusters, "\n")

    ## Número de clusters a usar: 50

    tryCatch({
      kmeans_result <- kmeans(pca_scores, centers = n_clusters, nstart = 25, iter.max = 100)
      cat("Clustering completado\n")
    }, error = function(e) {
      cat("Error en clustering, reduciendo número de clusters...\n")
      n_clusters <- min(20, nrow(presence))
      kmeans_result <- kmeans(pca_scores, centers = n_clusters, nstart = 10)
    })

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 2726950)
    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 2726950)
    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 2726950)

    ## Clustering completado

    cluster_assignments <- kmeans_result$cluster

    cluster_counts <- table(cluster_assignments)
    samples_per_cluster <- round((cluster_counts / sum(cluster_counts)) * nrow(presence))

    if(sum(samples_per_cluster) < nrow(presence)) {
      diff_needed <- nrow(presence) - sum(samples_per_cluster)
      largest_clusters <- order(cluster_counts, decreasing = TRUE)[1:diff_needed]
      samples_per_cluster[largest_clusters] <- samples_per_cluster[largest_clusters] + 1
    }

    background_indices <- c()

    for(i in 1:n_clusters) {
      if(samples_per_cluster[i] > 0) {
        cluster_pixels <- which(cluster_assignments == i)
        if(length(cluster_pixels) >= samples_per_cluster[i]) {
          sampled_indices <- sample(cluster_pixels, samples_per_cluster[i], replace = FALSE)
        } else {
          sampled_indices <- cluster_pixels
        }
        background_indices <- c(background_indices, sampled_indices)
      }
    }

    background <- climate_coords_clean[background_indices, ]
    colnames(background) <- c("x", "y")

    background <- as.data.frame(background)

    absence <- as.data.frame(background)

    colnames(absence) <- c("longitude", "latitude")
    absence$pa <- 0

    presence$latitude <- as.numeric(presence$latitude)
    presence$longitude <- as.numeric(presence$longitude)

    all_points <- rbind(presence, absence)

------------------------------------------------------------------------

## Section 6: Extract Environmental Values Using Buffers

    rad_to_deg <- function(rad) rad * 180 / pi
    deg_to_rad <- function(deg) deg * pi / 180

    earth_radius <- 6371000

    lat_degree <- 250 / (pi / 180 * earth_radius)

    lon_degree <- 250 / (pi / 180 * earth_radius * cos(deg_to_rad(presence$latitude)))

    avg_width <- (lat_degree + mean(lon_degree)) / 2

    presence_vect <- vect(presence, geom = c("longitude", "latitude"), crs = "epsg:4326")
    presence_buffer <- buffer(presence_vect, width = avg_width)

    background_vect <- vect(absence, geom = c("longitude", "latitude"), crs = "epsg:4326")
    background_buffer <- buffer(background_vect, width = avg_width)

    presence_climate <- extract(climate_data, presence_buffer, fun = mean, na.rm = TRUE)
    background_climate <- extract(climate_data, background_buffer, fun = mean, na.rm = TRUE)

    presence_climate <- cbind(presence, presence_climate)
    background_climate <- cbind(absence, background_climate)

    points_climate <- rbind(presence_climate, background_climate)

    drop_cols <- which(colnames(points_climate) %in% c("longitude", "latitude", "ID"))
    points_climate <- points_climate[, -drop_cols]

    response <- points_climate$pa

    names(climate_data)

    ## [1] "BIO2"                     "BIO5"                     "BIO6"                    
    ## [4] "BIO12"                    "aridityIndexThornthwaite" "climaticMoistureIndex"   
    ## [7] "growingDegDays5"          "PETseasonality"

    names(climate_data245)

    ## [1] "BIO2"                     "BIO5"                     "BIO6"                    
    ## [4] "BIO12"                    "aridityIndexThornthwaite" "climaticMoistureIndex"   
    ## [7] "growingDegDays5"          "PETseasonality"

    names(climate_data585)

    ## [1] "BIO2"                     "BIO5"                     "BIO6"                    
    ## [4] "BIO12"                    "aridityIndexThornthwaite" "climaticMoistureIndex"   
    ## [7] "growingDegDays5"          "PETseasonality"

------------------------------------------------------------------------

## Section 7: Crop Climate Data to Colombia Territory

    colombia_shp <- vect("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Colombia.shp")

    if (!crs(colombia_shp) == crs(climate_data)) {
      colombia_shp <- project(colombia_shp, crs(climate_data))
    }

    if (!crs(colombia_shp) == crs(climate_data245)) {
      colombia_shp <- project(colombia_shp, crs(climate_data245))
    }

    if (!crs(colombia_shp) == crs(climate_data585)) {
      colombia_shp <- project(colombia_shp, crs(climate_data585))
    }

    cropped <- crop(climate_data, colombia_shp)
    climate_data <- mask(cropped, colombia_shp)

    cropped245 <- crop(climate_data245, colombia_shp)
    climate_data245 <- mask(cropped245, colombia_shp)

    cropped585 <- crop(climate_data585, colombia_shp)
    climate_data585 <- mask(cropped585, colombia_shp)

------------------------------------------------------------------------

## Section 8: Load Calibrated MaxEnt Model

    mejor_modelo <- readRDS("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Paper_comnplejo_punta_morada/Maxent/Bactericera_colombia_modelo.rds")

------------------------------------------------------------------------

## Section 9: Load UPRA Potato Suitability Reference Map

    mapa_papa_sp <-readRDS("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/shape_unificado_aptitud_papa/mapa_papa_sp.rds")
    library(terra)
    library(tidyterra)
    library(ggplot2)

    mapa_papa_vect <- terra::vect(mapa_papa_sp)

    ggplot() +
      geom_spatvector(data = mapa_papa_vect,
                      fill = "#1B9E77",
                      color = "#1B9E77",
                      linewidth = 0.1) +
      scale_fill_princess_c(
        palette = "america",
        na.value = "transparent",
        name = "Suitability"
      ) +
      coord_sf() +   
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "grey90", color = NA),
        panel.grid.major = element_line(color = "white", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(face = "bold")
      ) +
      labs(title = "Area of aptitud in papa cultivation")

![](C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Paper_comnplejo_punta_morada/MESS_files/figure-markdown_strict/load-upra-1.png)

------------------------------------------------------------------------

## Section 10: MESS Analysis and Predictions

### Current Climate Scenario

    climate_data1  <- stack(climate_data)
    ref_vals1 <- extract(climate_data1, mapa_papa_sp)
    ref_df1 <- na.omit(as.data.frame(ref_vals1))
    mess_result1 <- mess(x = climate_data1, v = ref_df1, full = TRUE)
    rmess1 <- mess_result1[["rmess"]]
    climate_analogous1 <- mask(climate_data1, rmess1 >= 0)
    climate_analogous <- terra::rast(climate_analogous1)
    prediction_current <- predict(mejor_modelo, climate_analogous, type = "cloglog")
    plot(rmess1, main = "Current: Analysis of MESS-Analogous Zones")

![](C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Paper_comnplejo_punta_morada/MESS_files/figure-markdown_strict/mess-current-1.png)

### SSP2-4.5 Climate Scenario

    climate_data2  <- stack(climate_data245)
    ref_vals2  <- extract(climate_data2, mapa_papa_sp)
    ref_df2 <- na.omit(as.data.frame(ref_vals2))
    mess_result2  <- mess(x = climate_data2, v = ref_df2, full = TRUE)
    rmess2   <- mess_result2[["rmess"]]
    climate_analogous2  <- mask(climate_data2, rmess2 >= 0)
    climate_analogous2 <- terra::rast(climate_analogous2)
    prediction_245  <- predict(mejor_modelo, climate_analogous2, type = "cloglog")
    plot(rmess2, main = "245: Analysis of MESS-Analogous Zones")

![](C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Paper_comnplejo_punta_morada/MESS_files/figure-markdown_strict/mess-ssp245-1.png)

### SSP5-8.5 Climate Scenario

    climate_data3 <- stack(climate_data585)
    ref_vals3 <- extract(climate_data3, mapa_papa_sp)
    ref_df3 <- na.omit(as.data.frame(ref_vals3))
    mess_result3 <- mess(x = climate_data3, v = ref_df3, full = TRUE)
    rmess3 <- mess_result3[["rmess"]]
    climate_analogous3 <- mask(climate_data3, rmess3 >= 0)
    climate_analogous3  <- terra::rast(climate_analogous3)
    prediction_585   <- predict(mejor_modelo, climate_analogous3, type = "cloglog")
    plot(rmess3, main = "585: Analysis of MESS-Analogous Zones")

![](C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Paper_comnplejo_punta_morada/MESS_files/figure-markdown_strict/mess-ssp585-1.png)

------------------------------------------------------------------------

## Section 11: Visualization with ggplot2

    pred_df <- as.data.frame(prediction_current, xy = TRUE)

    scenario_current <- ggplot(pred_df, aes(x = x, y = y, fill = lyr1)) +
      geom_tile() +
      scale_fill_princess_c(
        palette = "america",
        na.value = "transparent",
        name = "Suitability"
      ) +
      coord_equal() +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "grey90", color = NA),
        panel.grid.major = element_line(color = "white", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(face = "bold")
      )

    print(scenario_current)

![](C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Paper_comnplejo_punta_morada/MESS_files/figure-markdown_strict/visualize-predictions-current-1.png)

    pred_df1 <- as.data.frame(prediction_245, xy = TRUE)

    scenario_245 <- ggplot(pred_df1, aes(x = x, y = y, fill = lyr1)) +
      geom_tile() +
      scale_fill_princess_c(
        palette = "america",
        na.value = "transparent",
        name = "Suitability"
      ) +
      coord_equal() +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "grey90", color = NA),
        panel.grid.major = element_line(color = "white", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(face = "bold")
      )

    print(scenario_245)

![](C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Paper_comnplejo_punta_morada/MESS_files/figure-markdown_strict/visualize-predictions-245-1.png)

    pred_df2 <- as.data.frame(prediction_585, xy = TRUE)

    scenario_585 <- ggplot(pred_df2, aes(x = x, y = y, fill = lyr1)) +
      geom_tile() +
      scale_fill_princess_c(
        palette = "america",
        na.value = "transparent",
        name = "Suitability"
      ) +
      coord_equal() +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "grey90", color = NA),
        panel.grid.major = element_line(color = "white", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(face = "bold")
      )

    print(scenario_585)

![](C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Paper_comnplejo_punta_morada/MESS_files/figure-markdown_strict/visualize-predictions-585-1.png)

------------------------------------------------------------------------

## Session Information

    sessionInfo()

    ## R version 4.3.3 (2024-02-29 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19045)
    ## 
    ## Matrix products: default
    ## 
    ## 
    ## locale:
    ## [1] LC_COLLATE=Spanish_Colombia.utf8  LC_CTYPE=Spanish_Colombia.utf8   
    ## [3] LC_MONETARY=Spanish_Colombia.utf8 LC_NUMERIC=C                     
    ## [5] LC_TIME=Spanish_Colombia.utf8    
    ## 
    ## time zone: America/Bogota
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] RColorBrewer_1.1-3           sf_1.0-20                    tidyterra_0.6.1             
    ##  [4] tmap_4.0                     randomForestExplainer_0.10.1 pROC_1.18.5                 
    ##  [7] dismo_1.3-14                 randomForest_4.7-1.1         rgl_1.3.18                  
    ## [10] png_0.1-8                    cluster_2.1.8.1              predicts_0.1-11             
    ## [13] geodata_0.6-2                raster_3.6-32                sp_2.2-0                    
    ## [16] lubridate_1.9.3              forcats_1.0.0                stringr_1.5.1               
    ## [19] purrr_1.0.2                  readr_2.1.5                  tidyr_1.3.1                 
    ## [22] tibble_3.3.0                 ggplot2_3.5.1                tidyverse_2.0.0             
    ## [25] dplyr_1.1.4                  readxl_1.4.5                 SDMtune_1.3.1               
    ## [28] terra_1.8-29                 furrr_0.3.1                  future_1.34.0               
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] DBI_1.2.3               tmaptools_3.2           gridExtra_2.3          
    ##  [4] s2_1.1.6                logger_0.4.0            rlang_1.1.3            
    ##  [7] magrittr_2.0.3          e1071_1.7-16            compiler_4.3.3         
    ## [10] vctrs_0.6.5             wk_0.9.1                pkgconfig_2.0.3        
    ## [13] crayon_1.5.3            fastmap_1.2.0           lwgeom_0.2-14          
    ## [16] leafem_0.2.5            labeling_0.4.3          rmarkdown_2.29         
    ## [19] tzdb_0.5.0              spacesXYZ_1.5-1         xfun_0.53              
    ## [22] jsonlite_2.0.0          R6_2.6.1                stringi_1.8.7          
    ## [25] GGally_2.2.1            parallelly_1.39.0       stars_0.6-8            
    ## [28] cellranger_1.1.0        Rcpp_1.0.12             knitr_1.50             
    ## [31] base64enc_0.1-3         leaflet.providers_2.0.0 timechange_0.3.0       
    ## [34] tidyselect_1.2.1        rstudioapi_0.17.1       dichromat_2.0-0.1      
    ## [37] abind_1.4-8             yaml_2.3.10             viridis_0.6.5          
    ## [40] codetools_0.2-19        listenv_0.9.1           leafsync_0.1.0         
    ## [43] lattice_0.22-5          plyr_1.8.9              withr_3.0.2            
    ## [46] evaluate_1.0.5          ggstats_0.9.0           units_0.8-5            
    ## [49] proxy_0.4-27            pillar_1.11.0           KernSmooth_2.23-22     
    ## [52] DT_0.33                 generics_0.1.4          hms_1.1.3              
    ## [55] scales_1.4.0            globals_0.16.3          class_7.3-22           
    ## [58] glue_1.7.0              tools_4.3.3             leaflegend_1.2.1       
    ## [61] data.table_1.17.8       XML_3.99-0.18           crosstalk_1.2.2        
    ## [64] colorspace_2.1-1        cols4all_0.8            cli_3.6.2              
    ## [67] viridisLite_0.4.2       gtable_0.3.6            digest_0.6.35          
    ## [70] classInt_0.4-10         ggrepel_0.9.5           htmlwidgets_1.6.4      
    ## [73] farver_2.1.2            htmltools_0.5.8.1       lifecycle_1.0.4        
    ## [76] leaflet_2.2.2           microbenchmark_1.5.0
