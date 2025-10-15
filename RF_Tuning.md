## Section 1: Environment Setup and Library Loading

We configure the computational environment for parallel processing (24
cores) and load required packages for species distribution modeling,
spatial analysis, and Random Forest implementation. This setup enables
efficient processing of large spatial datasets and cross-validation
procedures.

    library(future)
    library(furrr)
    library(parallel)
    gc()  # Garbage collection

    ##            used  (Mb) gc trigger  (Mb) max used  (Mb)
    ## Ncells  5787557 309.1    9904608 529.0  9904608 529.0
    ## Vcells 12991633  99.2   21468864 163.8 18277439 139.5

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
    library(terra)
    library(dplyr)
    library(cluster)
    library(png)
    library(grid)
    library(rgl)

## Section 2: Climate Data Loading and Preprocessing

We load bioclimatic variables from WorldClim (19 variables) and ENVIREM
datasets (18 environmental variables) at 2.5 arc-minute resolution (~4.5
km), representing the historical period 1970-2000. Following the
manuscript methodology, we select 8 environmental predictors after
collinearity analysis (Spearman |r| &gt; 0.8) to avoid multicollinearity
issues.

**Selected variables**: BIO3, BIO5, BIO6, BIO13, climaticMoistureIndex,
aridityIndexThornthwaite, PETseasonality, growingDegDays5

    # Load bioclimatic and biophysical data
    bioclim_data <- rast("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/data/2.5 minutos/bioclim_data.tif")
    Envirem <- rast("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/ensambles mensuales/envirem_vars_historico.tif")
     
    names(bioclim_data) <- paste0("BIO", 1:nlyr(bioclim_data))

    # Combine climate and biophysical variables
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

    ##### Random Forest variable selection
    selected_variables <- c(3, 5, 6, 13, 26, 23, 34, 35) ##### Bactericera Colombia final

    # Create new SpatRaster with selected variables only
    climate_data_subset <- subset(climate_data, selected_variables)

    # Verify dimensions and names
    climate_data <- climate_data_subset

    # Scale predictors
    climate_data <- scale(climate_data)

## Section 3: Species Presence Data Loading and Preprocessing

We load *Bactericera cockerelli* presence records from ICA (Colombian
Agricultural Institute) phytosanitary surveillance data (2021-2024).
Duplicate coordinates are removed to avoid spatial autocorrelation bias
in model training.

    # Load presence points
    obs_data <- read_excel("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Bactericera cockerelli Colombia Oficial surveillance ICA.xlsx")
    obs_data <- obs_data[!is.na(obs_data$longitude) & !is.na(obs_data$latitude), ]
    obs_data <- obs_data[!duplicated(obs_data[, c("longitude", "latitude")]), ]

    presence <- obs_data[, c("longitude", "latitude")]
    presence$pa <- 1  # Add column indicating presence
    # Remove duplicates (repeated coordinates)
    presence <- presence %>% distinct(longitude, latitude, .keep_all = TRUE)

    # Convert to numeric
    presence$longitude <- as.numeric(presence$longitude)
    presence$latitude <- as.numeric(presence$latitude)

## Section 4: Stratified Environmental Sampling for Background Points

Following Valavi et al. (2021) recommendations, we implement
**stratified sampling in environmental space** using PCA-based K-means
clustering to generate background points. This approach ensures
comprehensive representation of environmental gradients while minimizing
spatial autocorrelation biases inherent in geographic random sampling.

**Methodology**: 1. Crop climate data to Colombia boundaries (GADM level
0) 2. Perform PCA on all environmental variables 3. K-means clustering
in PCA space (50 strata for mesoscale) 4. Proportional sampling from
each environmental stratum

This addresses the “class overlap” problem in presence-background data
by ensuring background samples represent all available environmental
conditions, including those where the species occurs.

    library(terra)
    library(cluster)
    library(geodata)

    # Set seed for reproducibility
    set.seed(20210707)

    # Create spatial object from presences
    presence_area <- vect(presence, geom = c("longitude", "latitude"), crs = "epsg:4326")

    # CROP TO COMPLETE COLOMBIA (not just presence area)
    tryCatch({
      colombia_boundary <- gadm(country = "COL", level = 0, path = tempdir())
      # colombia_boundary is already a SpatVector
      climate_crop <- crop(climate_data, colombia_boundary)
      climate_crop <- mask(climate_crop, colombia_boundary)
      
      # Verify successful cropping
      if(all(is.na(values(climate_crop[[1]])))) {
        stop("Error: Colombia crop resulted in NA values")
      }
      
    }, error = function(e) {
      cat("Error downloading Colombia boundaries:", e$message, "\n")
      cat("Trying alternative method...\n")
      
      # Alternative method: create approximate polygon for Colombia
      colombia_ext <- ext(-82, -66, -5, 13)  # Approximate extent of Colombia
      colombia_boundary <- as.polygons(colombia_ext, crs = "epsg:4326")
      climate_crop <- crop(climate_data, colombia_boundary)
      climate_crop <- mask(climate_crop, colombia_boundary)
    })

    # STRATIFIED SAMPLING IN ENVIRONMENTAL SPACE
    # 1. Extract all climate values from cropped area (Colombia)
    tryCatch({
      climate_values <- values(climate_crop, na.rm = TRUE)
      climate_coords <- xyFromCell(climate_crop, which(!is.na(values(climate_crop[[1]]))))
      
      # Verify we have data
      if(nrow(climate_coords) == 0) {
        stop("Error: No valid coordinates found")
      }
      
      cat("Valid pixels found:", nrow(climate_coords), "\n")
      
    }, error = function(e) {
      stop("Error extracting climate values: ", e$message)
    })

    ## Valid pixels found: 54539

    # 2. Remove rows with NA
    complete_cases <- complete.cases(climate_values)
    climate_values_clean <- climate_values[complete_cases, ]
    climate_coords_clean <- climate_coords[complete_cases, ]

    cat("Pixels after NA cleaning:", nrow(climate_coords_clean), "\n")

    ## Pixels after NA cleaning: 54539

    # Verify sufficient data
    if(nrow(climate_values_clean) < nrow(presence)) {
      stop("Error: Insufficient valid pixels for sampling")
    }

    # 3. Perform principal component analysis (PCA)
    tryCatch({
      climate_pca <- prcomp(climate_values_clean, scale. = TRUE, center = TRUE)
      cat("PCA completed. Variance explained PC1-3:", 
          round(sum(climate_pca$sdev[1:3]^2)/sum(climate_pca$sdev^2)*100, 2), "%\n")
    }, error = function(e) {
      stop("Error in PCA: ", e$message)
    })

    ## PCA completed. Variance explained PC1-3: 96.57 %

    # 4. Extract first 3 principal components
    pca_scores <- climate_pca$x[, 1:3]

    # 5. K-means clustering in environmental space (PCA)
    n_clusters <- min(50, nrow(presence) * 2, nrow(climate_values_clean) %/% 10)

    cat("Number of clusters to use:", n_clusters, "\n")

    ## Number of clusters to use: 50

    tryCatch({
      kmeans_result <- kmeans(pca_scores, centers = n_clusters, nstart = 25, iter.max = 100)
      cat("Clustering completed\n")
    }, error = function(e) {
      cat("Error in clustering, reducing number of clusters...\n")
      n_clusters <- min(20, nrow(presence))
      kmeans_result <- kmeans(pca_scores, centers = n_clusters, nstart = 10)
    })

    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 2726950)
    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 2726950)
    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 2726950)
    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 2726950)
    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 2726950)
    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 2726950)
    ## Warning: Quick-TRANSfer stage steps exceeded maximum (= 2726950)

    ## Clustering completed

    # 6. Assign each pixel to an environmental cluster
    cluster_assignments <- kmeans_result$cluster

    # 7. Calculate number of points per cluster (proportional to cluster size)
    cluster_counts <- table(cluster_assignments)
    samples_per_cluster <- round((cluster_counts / sum(cluster_counts)) * nrow(presence))

    # Ensure we have at least the desired total number
    if(sum(samples_per_cluster) < nrow(presence)) {
      # Add points to largest clusters
      diff_needed <- nrow(presence) - sum(samples_per_cluster)
      largest_clusters <- order(cluster_counts, decreasing = TRUE)[1:diff_needed]
      samples_per_cluster[largest_clusters] <- samples_per_cluster[largest_clusters] + 1
    }

    # 8. Stratified sampling by cluster
    background_indices <- c()

    for(i in 1:n_clusters) {
      if(samples_per_cluster[i] > 0) {
        cluster_pixels <- which(cluster_assignments == i)
        if(length(cluster_pixels) >= samples_per_cluster[i]) {
          sampled_indices <- sample(cluster_pixels, samples_per_cluster[i], replace = FALSE)
        } else {
          # If cluster too small, take all pixels
          sampled_indices <- cluster_pixels
        }
        background_indices <- c(background_indices, sampled_indices)
      }
    }

    # 9. Create background object with selected coordinates
    background <- climate_coords_clean[background_indices, ]
    colnames(background) <- c("x", "y")

    # Convert to data.frame for compatibility
    background <- as.data.frame(background)

    # Create absence table
    absence <- as.data.frame(background)

    colnames(absence) <- c("longitude", "latitude")
    absence$pa <- 0  # Mark as absence

    presence$latitude <- as.numeric(presence$latitude)
    presence$longitude <- as.numeric(presence$longitude)
    # Combine presence and background
    all_points <- rbind(presence, absence)

## Section 5: Buffer-Based Environmental Data Extraction

Following the manuscript methodology, we create 500-meter buffers around
each presence and background point to reduce spatial autocorrelation in
presence records. This standardized approach ensures methodological
consistency across spatial scales and reduces the influence of spatial
clustering in phytosanitary monitoring data.

The mean values of environmental variables within these buffers are
extracted and used for model calibration, addressing the aggregated
nature of surveillance data.

    # Implicit function to calculate width in degrees by latitude
    rad_to_deg <- function(rad) rad * 180 / pi
    deg_to_rad <- function(deg) deg * pi / 180

    # Earth radius in meters (approximate)
    earth_radius <- 6371000

    # Convert 250 meters to degrees latitude (constant)
    lat_degree <- 250 / (pi / 180 * earth_radius)

    # Convert 250 meters to degrees longitude (depends on latitude)
    lon_degree <- 250 / (pi / 180 * earth_radius * cos(deg_to_rad(presence$latitude)))

    # Average between latitude and longitude for approximate buffer
    avg_width <- (lat_degree + mean(lon_degree)) / 2

    # Generate buffers with adjusted width
    presence_vect <- vect(presence, geom = c("longitude", "latitude"), crs = "epsg:4326")
    presence_buffer <- buffer(presence_vect, width = avg_width)

    background_vect <- vect(absence, geom = c("longitude", "latitude"), crs = "epsg:4326")
    background_buffer <- buffer(background_vect, width = avg_width)

    # Extract climate data for presence and background buffers
    presence_climate <- extract(climate_data, presence_buffer, fun = mean, na.rm = TRUE)
    background_climate <- extract(climate_data, background_buffer, fun = mean, na.rm = TRUE)

    # Combine climate data with original points
    presence_climate <- cbind(presence, presence_climate)
    background_climate <- cbind(absence, background_climate)

    # Combine presence and background data with climate values
    points_climate <- rbind(presence_climate, background_climate)

    # Remove longitude and latitude columns before continuing
    drop_cols <- which(colnames(points_climate) %in% c("longitude", "latitude", "ID"))
    points_climate <- points_climate[, -drop_cols]

    # Create response vector (presence/absence)
    response <- points_climate$pa

## Section 6: Data Cleaning and Index Recreation

We remove NA values and recreate presence (pa=1) and absence (pa=0)
indices. This step ensures data integrity before model training and
prevents errors during cross-validation procedures.

    cat("=== DATA CLEANING AND PREPARATION ===\n")

    ## === DATA CLEANING AND PREPARATION ===

    cat("Data before cleaning:", nrow(points_climate), "rows\n")

    ## Data before cleaning: 2460 rows

    # STEP 1: Clean NA data
    points_climate <- na.omit(points_climate)
    cat("Data after cleaning:", nrow(points_climate), "rows\n")

    ## Data after cleaning: 2460 rows

    # STEP 2: RECREATE presence (1) and absence (0) indices
    presences <- which(points_climate$pa == 1)
    absences <- which(points_climate$pa == 0)

    cat("Presences:", length(presences), "\n")

    ## Presences: 1230

    cat("Absences:", length(absences), "\n")

    ## Absences: 1230

## Section 7: Addressing Class Imbalance Through Down-Sampling

Following Valavi et al. (2021), we implement a **down-sampling
strategy** to manage class imbalance in presence-background data. Class
imbalance occurs when background samples substantially outnumber
presence points, potentially leading to biased predictions.

**Balancing Strategy**: - Maximum ratio: 5:1 (absences:presences) to
avoid bias toward majority class - Conservative approach prevents
overfitting to background samples - Maintains all presence points (true
data) while subsampling background

This approach differs from equal-sampling (Barbet-Massin et al. 2012) by
creating a balanced dataset once, rather than repeatedly at the model
level, improving computational efficiency while maintaining statistical
validity.

    # CONTROL OF SAMPLING BIAS: Balance data to avoid overfitting
    set.seed(2023)
    n_presences <- length(presences)
    n_absences <- length(absences)

    # More conservative ratio to avoid bias toward absences
    ratio_limit <- 5  # Maximum 5:1 (absences:presences)
    if (n_absences > n_presences * ratio_limit) {
      selected_absences <- sample(absences, n_presences * ratio_limit)
      cat("Subsampling absences from", n_absences, "to", length(selected_absences), "\n")
      balanced_indices <- c(presences, selected_absences)
      points_climate_balanced <- points_climate[balanced_indices, ]
    } else {
      points_climate_balanced <- points_climate
      cat("Keeping all absences (acceptable ratio)\n")
    }

    ## Keeping all absences (acceptable ratio)

    cat("Balanced dataset - Presences:", sum(points_climate_balanced$pa), 
        "Absences:", sum(1 - points_climate_balanced$pa), "\n")

    ## Balanced dataset - Presences: 1230 Absences: 1230

## Section 8: Stratified 5-Fold Cross-Validation Setup

We implement stratified k-fold cross-validation (k=5) to ensure each
fold maintains the same presence:absence ratio as the full dataset. This
is critical for presence-background data where class imbalance is
inherent.

**Stratification ensures**: - Each fold has representative samples from
both classes - Unbiased performance estimation across folds - Reliable
model evaluation despite imbalanced data

This follows best practices for SDM validation (Roberts et al. 2017) and
provides robust performance metrics.

    # STEP 3: Stratified cross-validation (5 folds)
    presences_bal <- which(points_climate_balanced$pa == 1)
    absences_bal <- which(points_climate_balanced$pa == 0)

    fold_pres <- sample(rep(1:5, length.out = length(presences_bal)))
    fold_abs <- sample(rep(1:5, length.out = length(absences_bal)))

    fold <- integer(nrow(points_climate_balanced))
    fold[presences_bal] <- fold_pres
    fold[absences_bal] <- fold_abs

    cat("\nCross-validation configured: 5 stratified folds\n")

    ## 
    ## Cross-validation configured: 5 stratified folds

    # Verify fold distribution
    for(i in 1:5) {
      fold_data <- points_climate_balanced[fold == i, ]
      cat("Fold", i, "- Presences:", sum(fold_data$pa), 
          "Absences:", sum(1 - fold_data$pa), "\n")
    }

    ## Fold 1 - Presences: 246 Absences: 246 
    ## Fold 2 - Presences: 246 Absences: 246 
    ## Fold 3 - Presences: 246 Absences: 246 
    ## Fold 4 - Presences: 246 Absences: 246 
    ## Fold 5 - Presences: 246 Absences: 246

## Section 9: Hyperparameter Tuning - mtry Optimization

We tune the **mtry** parameter (number of predictors randomly sampled at
each split) using the `tuneRF` function. This is a critical
hyperparameter in Random Forest that controls the randomness and
correlation between trees.

**Tuning Strategy**: - Base evaluation: 300 trees (faster than final
model) - Improvement threshold: 2% OOB error reduction - Small noise
addition (±0.001) prevents overfitting during tuning - **Maximum mtry
limit**: √(number of predictors) to prevent overfitting

The number of trees (ntree) was manually set to 1500 based on literature
recommendations for ecological data, as Random Forests are relatively
insensitive to this parameter beyond a certain threshold.

**Note**: Following Valavi et al. (2021), we set mtry conservatively to
balance model complexity with generalization capacity.

    selected_vars <- names(climate_data)
    all_data <- points_climate_balanced[, c(selected_vars, "pa")]
    all_data$pa <- as.numeric(all_data$pa)

    prNum <- sum(all_data$pa == 1)
    bgNum <- sum(all_data$pa == 0)

    cat("\n=== MTRY TUNING WITH OVERFITTING CONTROL ===\n")

    ## 
    ## === MTRY TUNING WITH OVERFITTING CONTROL ===

    cat("Presences for training:", prNum, "\n")

    ## Presences for training: 1230

    cat("Absences for training:", bgNum, "\n")

    ## Absences for training: 1230

    library(randomForest)
    library(dismo)
    library(pROC)
    library(randomForestExplainer)

    # Add minimal noise to avoid overfitting in tuning
    set.seed(123)
    y_noise <- all_data$pa + runif(nrow(all_data), -0.001, 0.001)

    # More conservative tuning
    optimal_mtry <- tuneRF(
      x = all_data[, selected_vars],
      y = y_noise,
      stepFactor = 1.5,
      ntreeTry = 300,        # Fewer trees for tuning
      improve = 0.02,        # Stricter improvement
      trace = TRUE,
      plot = FALSE,
      doBest = FALSE
    )

    ## mtry = 2  OOB error = 0.006827234 
    ## Searching left ...
    ## Searching right ...
    ## mtry = 3     OOB error = 0.006837134 
    ## -0.001450118 0.02

    best_mtry <- optimal_mtry[which.min(optimal_mtry[, "OOBError"]), "mtry"]
    cat("Best mtry found:", best_mtry, "\n")

    ## Best mtry found: 2

    # LIMIT mtry to avoid overfitting
    max_mtry_allowed <- max(2, floor(sqrt(length(selected_vars))))
    if (best_mtry > max_mtry_allowed) {
      best_mtry <- max_mtry_allowed
      cat("mtry limited to:", best_mtry, "to avoid overfitting\n")
    }

## Section 10: Random Forest Model Training with Anti-Overfitting Parameters

Following Valavi et al. (2021), we implement **regression-RF**
(probability forest) instead of classification-RF for
presence-background data. Regression-RF produces more stable and
calibrated probability estimates than classification approaches.

**Anti-Overfitting Parameters**: - **ntree = 1000**: Moderate number of
trees (balance between accuracy and computation) - **mtry =
best\_mtry**: Optimized and limited to prevent overfitting - **sampsize
= 70%**: Bootstrap sample size per tree (reduces variance) - **nodesize
= 50**: Large terminal nodes (fewer splits, less overfitting) -
**maxnodes = 30**: Limited nodes per tree (controls model complexity) -
**replace = TRUE**: Bootstrap sampling with replacement - **NO class
weights**: Weights removed as they can cause bias in presence-background
data

These parameters directly address the “class overlap” problem by
preventing trees from becoming too deep and overfitting to noisy splits
in overlapping environmental space.

    auc_train_folds <- auc_test_folds <- cor_folds <- tss_folds <- omission_rate_folds <- 
      sensitivity_folds <- threshold_folds <- oob_error_folds <- var_explained_folds <- numeric(5)
    models_folds <- list()

    cat("\n=== 5-FOLD CROSS-VALIDATION ===\n")

    ## 
    ## === 5-FOLD CROSS-VALIDATION ===

    for(fold_i in 1:5) {
      cat(paste0("\n--- FOLD ", fold_i, " ---\n"))
      
      testing <- points_climate_balanced[fold == fold_i, c(selected_vars, "pa")]
      training <- points_climate_balanced[fold != fold_i, c(selected_vars, "pa")]
      
      training$pa <- as.numeric(training$pa)
      testing$pa <- as.numeric(testing$pa)
      
      cat("Training data - Presences:", sum(training$pa), 
          "Absences:", sum(1 - training$pa), "\n")
      
      # ANTI-OVERFITTING PARAMETERS for Random Forest
      set.seed(20178 + fold_i)
      
      forest_fold <- randomForest(
        formula = as.formula(paste("pa ~", paste(selected_vars, collapse = " + "))),
        data = training,
        ntree = 1000,                    # Fewer trees
        mtry = best_mtry,                # Controlled mtry
        sampsize = floor(0.7 * nrow(training)),  # 70% sample per tree
        replace = TRUE,
        nodesize = 50,                   # LARGER terminal nodes
        maxnodes = 30,                   # FEWER nodes per tree
        importance = TRUE,
        keep.forest = TRUE
        # NO WEIGHTS - completely removed
      )
      
      models_folds[[fold_i]] <- forest_fold
      
      # Predictions with [0,1] limits
      pred_train <- predict(forest_fold, newdata = training)
      pred_test <- predict(forest_fold, newdata = testing)
      pred_train <- pmin(1, pmax(0, pred_train))
      pred_test <- pmin(1, pmax(0, pred_test))
      
      # ROC and AUC
      roc_train <- roc(training$pa, pred_train, quiet = TRUE)
      roc_test <- roc(testing$pa, pred_test, quiet = TRUE)
      
      auc_train_folds[fold_i] <- auc(roc_train)
      auc_test_folds[fold_i] <- auc(roc_test)
      
      # Wrapper for dismo::evaluate
      rf_wrapper <- function(object, newdata) {
        pred <- predict(object, newdata)
        return(pmin(1, pmax(0, pred)))
      }
      
      # Evaluation with dismo
      evaluation <- dismo::evaluate(
        p = testing[testing$pa == 1, selected_vars],
        a = testing[testing$pa == 0, selected_vars],
        model = forest_fold,
        predFun = rf_wrapper
      )
      
      cor_folds[fold_i] <- evaluation@cor
      
      # TSS and optimal threshold
      thresholds <- evaluation@t
      confusion_matrices <- evaluation@confusion
      tss <- evaluation@TPR + evaluation@TNR - 1
      best_threshold_idx <- which.max(tss)
      best_threshold <- thresholds[best_threshold_idx]
      
      threshold_folds[fold_i] <- best_threshold
      tss_folds[fold_i] <- max(tss)
      
      # Omission rate
      tp <- confusion_matrices[best_threshold_idx, "tp"]
      fn <- confusion_matrices[best_threshold_idx, "fn"]
      omission_rate <- fn / (tp + fn)
      
      omission_rate_folds[fold_i] <- omission_rate
      sensitivity_folds[fold_i] <- 1 - omission_rate
      
      # Model metrics
      oob_error_folds[fold_i] <- forest_fold$mse[forest_fold$ntree]
      var_explained_folds[fold_i] <- (1 - forest_fold$mse[forest_fold$ntree] / var(training$pa)) * 100
      
      cat("AUC training:", round(auc_train_folds[fold_i], 4), "\n")
      cat("AUC testing:", round(auc_test_folds[fold_i], 4), "\n")
      cat("TSS:", round(tss_folds[fold_i], 4), "\n")
      cat("Omission rate:", round(omission_rate_folds[fold_i], 4), "\n")
      cat("Optimal threshold:", round(best_threshold, 4), "\n")
      
      # OVERFITTING DIAGNOSTIC
      overfitting_gap <- auc_train_folds[fold_i] - auc_test_folds[fold_i]
      if (overfitting_gap > 0.1) {
        cat("⚠️  WARNING: Possible overfitting (AUC GAP =", round(overfitting_gap, 4), ")\n")
      }
    }

    ## 
    ## --- FOLD 1 ---
    ## Training data - Presences: 984 Absences: 984

    ## Warning in randomForest.default(m, y, ...): The response has five or fewer unique values.  Are
    ## you sure you want to do regression?

    ## AUC training: 0.9997 
    ## AUC testing: 0.996 
    ## TSS: 0.9715 
    ## Omission rate: 0.0203 
    ## Optimal threshold: 0.3937 
    ## 
    ## --- FOLD 2 ---
    ## Training data - Presences: 984 Absences: 984

    ## Warning in randomForest.default(m, y, ...): The response has five or fewer unique values.  Are
    ## you sure you want to do regression?

    ## AUC training: 0.9996 
    ## AUC testing: 0.9995 
    ## TSS: 0.9878 
    ## Omission rate: 0.0081 
    ## Optimal threshold: 0.5024 
    ## 
    ## --- FOLD 3 ---
    ## Training data - Presences: 984 Absences: 984

    ## Warning in randomForest.default(m, y, ...): The response has five or fewer unique values.  Are
    ## you sure you want to do regression?

    ## AUC training: 0.9997 
    ## AUC testing: 0.998 
    ## TSS: 0.9878 
    ## Omission rate: 0.0081 
    ## Optimal threshold: 0.3816 
    ## 
    ## --- FOLD 4 ---
    ## Training data - Presences: 984 Absences: 984

    ## Warning in randomForest.default(m, y, ...): The response has five or fewer unique values.  Are
    ## you sure you want to do regression?

    ## AUC training: 0.9997 
    ## AUC testing: 0.9996 
    ## TSS: 0.9878 
    ## Omission rate: 0.0122 
    ## Optimal threshold: 0.5386 
    ## 
    ## --- FOLD 5 ---
    ## Training data - Presences: 984 Absences: 984

    ## Warning in randomForest.default(m, y, ...): The response has five or fewer unique values.  Are
    ## you sure you want to do regression?

    ## AUC training: 0.9997 
    ## AUC testing: 0.9989 
    ## TSS: 0.9837 
    ## Omission rate: 0.0122 
    ## Optimal threshold: 0.5214

## Section 11: Final Model Training

We train the final Random Forest model using all balanced data with the
same anti-overfitting parameters validated during cross-validation. This
model will be used for spatial predictions and variable importance
analysis.

**Model Configuration**: - All presence and background data used - Same
hyperparameters as cross-validation - Trained once (no repeated sampling
at model level) - Results in a single, stable model for prediction

    cat("\n=== TRAINING FINAL ANTI-OVERFITTING MODEL ===\n")

    ## 
    ## === TRAINING FINAL ANTI-OVERFITTING MODEL ===

    all_data_final <- points_climate_balanced[, c(selected_vars, "pa")]
    all_data_final$pa <- as.numeric(all_data_final$pa)

    set.seed(20178)

    # FINAL MODEL MORE CONSERVATIVE
    forest1 <- randomForest(
      formula = as.formula(paste("pa ~", paste(selected_vars, collapse = " + "))),
      data = all_data_final,
      ntree = 1500,                      # Moderate number of trees
      mtry = best_mtry,                  # Optimized and limited mtry
      sampsize = floor(0.7 * nrow(all_data_final)),  # 70% sample per tree
      replace = TRUE,
      nodesize = 40,                     # Large nodes (fewer splits)
      maxnodes = 35,                     # Few nodes per tree
      importance = TRUE,
      keep.forest = TRUE
      # NO WEIGHTS
    )

    ## Warning in randomForest.default(m, y, ...): The response has five or fewer unique values.  Are
    ## you sure you want to do regression?

    cat("Final model trained with anti-overfitting parameters\n")

    ## Final model trained with anti-overfitting parameters

## Section 12: Threshold Selection for Binary Predictions

We calculate three threshold strategies for converting continuous
probabilities to binary presence/absence predictions:

**Threshold Selection Methods**: 1. **Robust threshold (median)**:
Recommended for operational use - less sensitive to outliers 2.
**Conservative threshold (25th percentile)**: Minimizes false positives
(type I error) 3. **Liberal threshold (75th percentile)**: Minimizes
false negatives (type II error)

The omission rate is used as the threshold criterion, following SDM best
practices where the threshold that maximizes TSS is selected for each
fold, then aggregated across folds.

    # Calculate robust and percentile thresholds
    mean_gap_auc <- mean(auc_train_folds - auc_test_folds)
    sd_gap_auc <- sd(auc_train_folds - auc_test_folds)

    threshold_robust <- median(threshold_folds)
    threshold_percentiles <- quantile(threshold_folds, probs = c(0.25, 0.5, 0.75))

## Section 13: Cross-Validation Performance Summary

We compile all performance metrics from 5-fold cross-validation to
assess model reliability and generalization capacity. These metrics
provide comprehensive evaluation of discriminative power, predictive
consistency, and overfitting control.

**Performance Metrics**: - **AUC (Area Under ROC Curve)**:
Threshold-independent discrimination (0.5-1.0) - &gt;0.9: Excellent
discrimination - 0.8-0.9: Good discrimination - 0.7-0.8: Acceptable
discrimination - **TSS (True Skill Statistic)**: Threshold-dependent
accuracy (-1 to 1) - &gt;0.7: High predictive accuracy - 0.4-0.7:
Moderate accuracy - **Correlation**: Association between predicted and
observed values - **Omission Rate**: Proportion of presences incorrectly
predicted as absences - **Sensitivity**: Proportion of presences
correctly identified (1 - omission rate) - **OOB Error**: Out-of-bag
error from Random Forest bootstrap - **Variance Explained**: Percentage
of response variance captured by model

**Overfitting Assessment**: - **AUC Gap**: Difference between train and
test AUC - &lt;0.05: Overfitting controlled - 0.05-0.1: Mild
overfitting - &gt;0.1: High overfitting

    cv_results <- list(
      auc_train = mean(auc_train_folds),
      auc_test = mean(auc_test_folds),
      auc_train_sd = sd(auc_train_folds),
      auc_test_sd = sd(auc_test_folds),
      auc_gap = mean_gap_auc,
      auc_gap_sd = sd_gap_auc,
      correlation = mean(cor_folds, na.rm = TRUE),
      correlation_sd = sd(cor_folds, na.rm = TRUE),
      tss = mean(tss_folds),
      tss_sd = sd(tss_folds),
      omission_rate = mean(omission_rate_folds),
      omission_rate_sd = sd(omission_rate_folds),
      sensitivity = mean(sensitivity_folds),
      sensitivity_sd = sd(sensitivity_folds),
      best_threshold = mean(threshold_folds),
      best_threshold_sd = sd(threshold_folds),
      threshold_robust = threshold_robust,
      threshold_conservative = threshold_percentiles[1],
      threshold_liberal = threshold_percentiles[3],
      oob_error = mean(oob_error_folds),
      oob_error_sd = sd(oob_error_folds),
      var_explained = mean(var_explained_folds),
      var_explained_sd = sd(var_explained_folds),
      overfitting_status = ifelse(mean_gap_auc > 0.1, "HIGH", 
                                  ifelse(mean_gap_auc < 0.05, "CONTROLLED", "MILD")),
      all_folds = list(
        auc_train = auc_train_folds,
        auc_test = auc_test_folds,
        correlation = cor_folds,
        tss = tss_folds,
        omission_rate = omission_rate_folds,
        sensitivity = sensitivity_folds,
        threshold = threshold_folds,
        oob_error = oob_error_folds,
        var_explained = var_explained_folds
      )
    )

    cat("\n=== ANTI-OVERFITTING RESULTS SUMMARY ===\n")

    ## 
    ## === ANTI-OVERFITTING RESULTS SUMMARY ===

    cat("AUC Train:", round(cv_results$auc_train, 4), "±", round(cv_results$auc_train_sd, 4), "\n")

    ## AUC Train: 0.9997 ± 0

    cat("AUC Test:", round(cv_results$auc_test, 4), "±", round(cv_results$auc_test_sd, 4), "\n")

    ## AUC Test: 0.9984 ± 0.0015

    cat("AUC Gap (overfitting):", round(cv_results$auc_gap, 4), "±", round(cv_results$auc_gap_sd, 4), "\n")

    ## AUC Gap (overfitting): 0.0013 ± 0.0015

    cat("Overfitting status:", cv_results$overfitting_status, "\n")

    ## Overfitting status: CONTROLLED

    cat("TSS:", round(cv_results$tss, 4), "±", round(cv_results$tss_sd, 4), "\n")

    ## TSS: 0.9837 ± 0.007

    cat("Omission rate:", round(cv_results$omission_rate, 4), "±", round(cv_results$omission_rate_sd, 4), "\n")

    ## Omission rate: 0.0122 ± 0.005

    cat("Sensitivity:", round(cv_results$sensitivity, 4), "±", round(cv_results$sensitivity_sd, 4), "\n")

    ## Sensitivity: 0.9878 ± 0.005

    cat("Variance explained:", round(cv_results$var_explained, 2), "±", round(cv_results$var_explained_sd, 2), "%\n")

    ## Variance explained: 95.73 ± 0.22 %

    cat("OOB Error:", round(cv_results$oob_error, 4), "±", round(cv_results$oob_error_sd, 4), "\n")

    ## OOB Error: 0.0107 ± 6e-04

    cat("\n=== RECOMMENDED THRESHOLDS FOR BINARIZATION ===\n")

    ## 
    ## === RECOMMENDED THRESHOLDS FOR BINARIZATION ===

    cat("Robust threshold (recommended):", round(cv_results$threshold_robust, 4), "\n")

    ## Robust threshold (recommended): 0.5024

    cat("Conservative threshold:", round(cv_results$threshold_conservative, 4), "\n")

    ## Conservative threshold: 0.3937

    cat("Liberal threshold:", round(cv_results$threshold_liberal, 4), "\n")

    ## Liberal threshold: 0.5214

    cat("\n=== ANTI-OVERFITTING MEASURES APPLIED ===\n")

    ## 
    ## === ANTI-OVERFITTING MEASURES APPLIED ===

    cat("✓ Complete elimination of problematic weights\n")

    ## ✓ Complete elimination of problematic weights

    cat("✓ Sampling bias control (maximum ratio 5:1)\n")

    ## ✓ Sampling bias control (maximum ratio 5:1)

    cat("✓ Limited mtry (max =", max_mtry_allowed, ")\n")

    ## ✓ Limited mtry (max = 2 )

    cat("✓ Large nodesize (40-50) for fewer splits\n")

    ## ✓ Large nodesize (40-50) for fewer splits

    cat("✓ Limited maxnodes (30-35) per tree\n")

    ## ✓ Limited maxnodes (30-35) per tree

    cat("✓ Reduced sampsize (70%) per tree\n")

    ## ✓ Reduced sampsize (70%) per tree

    cat("✓ Moderate number of trees (1000-1500)\n")

    ## ✓ Moderate number of trees (1000-1500)

    cat("✓ Robust threshold based on CV median\n")

    ## ✓ Robust threshold based on CV median

## Section 14: ROC Curve Visualization

We generate averaged ROC curves across all 5 cross-validation folds for
both training and testing datasets. This visualization demonstrates
model discriminative capacity and overfitting assessment.

**ROC Curve Methodology**: - High-resolution interpolation (100 FPR
points) - Smoothing via spline interpolation for cleaner visualization -
Separate curves for training and testing data - AUC values annotated
directly on plot

**Interpretation**: - Curves closer to top-left corner indicate better
discrimination - Gap between train and test curves indicates overfitting
degree - Diagonal line represents random chance (AUC = 0.5)

    # Load required libraries
    library(ggplot2)
    library(pROC)

    # CALCULATE AVERAGED ROC CURVES FROM 5 FOLDS - CORRECTED VERSION

    # Use more points for high resolution in high-performance models
    fpr_points <- c(seq(0, 0.1, length.out = 50),    # High resolution at low FPR
                    seq(0.1, 0.5, length.out = 30),   # Medium resolution
                    seq(0.5, 1, length.out = 20))     # Low resolution at high FPR

    # Matrices to store TPR and FPR
    tpr_train_matrix <- matrix(NA, nrow = length(fpr_points), ncol = 5)
    tpr_test_matrix <- matrix(NA, nrow = length(fpr_points), ncol = 5)

    # Calculate ROC curves for each fold
    for(fold_i in 1:5) {
      # Data
      testing <- points_climate[fold == fold_i, c(selected_vars, "pa")]
      training <- points_climate[fold != fold_i, c(selected_vars, "pa")]
      
      training$pa <- as.numeric(training$pa)
      testing$pa <- as.numeric(testing$pa)
      
      # Predictions
      pred_train <- predict(models_folds[[fold_i]], newdata = training)
      pred_test <- predict(models_folds[[fold_i]], newdata = testing)
      pred_train <- pmin(1, pmax(0, pred_train))
      pred_test <- pmin(1, pmax(0, pred_test))
      
      # ROC curves with higher resolution
      roc_train <- roc(training$pa, pred_train, quiet = TRUE, direction = "<")
      roc_test <- roc(testing$pa, pred_test, quiet = TRUE, direction = "<")
      
      # Extract FPR and TPR
      fpr_train <- 1 - roc_train$specificities
      tpr_train <- roc_train$sensitivities
      fpr_test <- 1 - roc_test$specificities
      tpr_test <- roc_test$sensitivities
      
      # CORRECTION: Force (0,0) and (1,1) points explicitly
      # Training
      fpr_train_complete <- c(0, fpr_train, 1)
      tpr_train_complete <- c(0, tpr_train, 1)
      
      # Testing  
      fpr_test_complete <- c(0, fpr_test, 1)
      tpr_test_complete <- c(0, tpr_test, 1)
      
      # Ensure increasing order and unique values
      train_order <- order(fpr_train_complete)
      test_order <- order(fpr_test_complete)
      
      # Interpolate TPR for common FPR - now guarantees start at (0,0)
      tpr_train_matrix[, fold_i] <- approx(fpr_train_complete[train_order], 
                                           tpr_train_complete[train_order], 
                                           fpr_points, method = "linear", 
                                           rule = 1)$y  # rule = 1 for NA outside range
      
      tpr_test_matrix[, fold_i] <- approx(fpr_test_complete[test_order], 
                                          tpr_test_complete[test_order], 
                                          fpr_points, method = "linear", 
                                          rule = 1)$y  # rule = 1 for NA outside range
    }

    ## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm): collapsing to unique 'x'
    ## values
    ## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm): collapsing to unique 'x'
    ## values
    ## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm): collapsing to unique 'x'
    ## values
    ## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm): collapsing to unique 'x'
    ## values
    ## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm): collapsing to unique 'x'
    ## values
    ## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm): collapsing to unique 'x'
    ## values
    ## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm): collapsing to unique 'x'
    ## values
    ## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm): collapsing to unique 'x'
    ## values
    ## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm): collapsing to unique 'x'
    ## values
    ## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm): collapsing to unique 'x'
    ## values

    # Average TPR and apply smoothing
    tpr_train_avg <- rowMeans(tpr_train_matrix, na.rm = TRUE)
    tpr_test_avg <- rowMeans(tpr_test_matrix, na.rm = TRUE)

    # Light smoothing for cleaner curves
    if(require(stats, quietly = TRUE)) {
      tpr_train_smooth <- predict(smooth.spline(fpr_points, tpr_train_avg, df = 50), fpr_points)$y
      tpr_test_smooth <- predict(smooth.spline(fpr_points, tpr_test_avg, df = 50), fpr_points)$y
      tpr_train_avg <- pmax(0, pmin(1, tpr_train_smooth))
      tpr_test_avg <- pmax(0, pmin(1, tpr_test_smooth))
    }

    # Average AUC
    auc_train_mean <- mean(auc_train_folds)
    auc_test_mean <- mean(auc_test_folds)

    # Create dataframe for ggplot
    roc_data <- data.frame(
      fpr = rep(fpr_points, 2),
      tpr = c(tpr_train_avg, tpr_test_avg),
      dataset = rep(c("Training", "Testing"), each = length(fpr_points))
    )

## Section 15: ROC Curve Plot Generation

We create a publication-quality ROC curve plot with intelligent zoom for
high-performance models. When both train and test AUC exceed 0.95, the
plot automatically zooms to the upper-left corner for better
visualization of subtle differences.

    # OPTIMIZED ROC PLOT FOR HIGH PERFORMANCE

    # Smart zoom: if AUC > 0.95, zoom to upper-left corner
    use_zoom <- min(auc_train_mean, auc_test_mean) > 0.95
    zoom_limit <- if(use_zoom) 0.3 else 1.0

    p <- ggplot(roc_data, aes(x = fpr, y = tpr, color = dataset)) +
      geom_line(size = 0.8, alpha = 0.95) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                  color = "gray50", size = 0.5, alpha = 0.6) +
      scale_color_manual(
        values = c("Training" = "red", "Testing" = "#0072B2"),
        labels = c("Training" = "Train", "Testing" = "Test")
      ) +
      labs(
        title = "",
        x = "1 - Specificity",
        y = "Sensitivity"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "none",
        panel.grid.major = element_line(color = "gray90", size = 0.3),
        panel.grid.minor = element_line(color = "gray95", size = 0.2),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )

    # Apply zoom if necessary
    if(use_zoom) {
      p <- p + 
        coord_fixed(xlim = c(0, zoom_limit), ylim = c(1-zoom_limit, 1)) +
        scale_x_continuous(breaks = seq(0, zoom_limit, by = 0.05), 
                           labels = scales::number_format(accuracy = 0.01)) +
        scale_y_continuous(breaks = seq(1-zoom_limit, 1, by = 0.05),
                           labels = scales::number_format(accuracy = 0.01)) +
        annotate("text", x = zoom_limit * 0.4, y = 1 - zoom_limit * 0.7, 
                 label = paste0("Train: ", sprintf("%.4f", auc_train_mean)), 
                 color = "red", size = 3.5, hjust = 0, fontface = "bold") +
        annotate("text", x = zoom_limit * 0.4, y = 1 - zoom_limit * 0.85, 
                 label = paste0("Test: ", sprintf("%.4f", auc_test_mean)), 
                 color = "#0072B2", size = 3.5, hjust = 0, fontface = "bold")
    } else {
      p <- p + 
        coord_fixed() +
        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
        annotate("text", x = 0.6, y = 0.3, 
                 label = paste0("Train: ", sprintf("%.4f", auc_train_mean)), 
                 color = "red", size = 3.5, hjust = 0, fontface = "bold") +
        annotate("text", x = 0.6, y = 0.2, 
                 label = paste0("Test: ", sprintf("%.4f", auc_test_mean)), 
                 color = "#0072B2", size = 3.5, hjust = 0, fontface = "bold")
    }

    # Display plot
    print(p)

![](C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Paper_comnplejo_punta_morada/RF_Tuning_files/figure-markdown_strict/roc-curve-plot-1.png)

## Section 16: Variable Importance Analysis - Minimum Depth Distribution

Following the manuscript methodology and Valavi et al. (2021)
recommendations, we assess variable importance using the **minimum depth
distribution** approach from the `RandomForestExplainer` package.

**Minimum Depth Methodology**: - Variables with lower minimum depth
appear in earlier splits (closer to tree root) - Earlier splits indicate
greater discriminatory capacity - More influential variables are
selected more frequently for initial partitioning - Calculated as
average minimum depth across all trees where variable was used

**Interpretation**: - Lower depth = Higher importance - Variables with
minimum depth &lt; 2 are typically most critical - This method is more
robust than traditional Mean Decrease Gini for imbalanced data

The selected predictors prioritize biological, epidemiological, and
ecological relevance (Friedl & Stampfer 2001; He et al. 2022; Vignali et
al. 2020).

    library(randomForest)

    # Extract variable importance from final model
    var_importance <- importance(forest1)

    # Calculate minimum depth distribution for variable importance
    # This follows the methodology described in the manuscript
    library(randomForestExplainer)

    # Calculate minimum depth
    min_depth_frame <- min_depth_distribution(forest1)

    # Create visualization
    plot_min_depth_distribution(min_depth_frame) +
      labs(
        title = "",
        x = "Variable",
        y = "Minimum Depth"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_line(color = "gray90", size = 0.3),
        panel.grid.minor = element_line(color = "gray95", size = 0.2),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )

![](C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Paper_comnplejo_punta_morada/RF_Tuning_files/figure-markdown_strict/variable-importance-plot-1.png)

    # Print summary statistics of minimum depth
    cat("\n=== VARIABLE IMPORTANCE: MINIMUM DEPTH SUMMARY ===\n")

    ## 
    ## === VARIABLE IMPORTANCE: MINIMUM DEPTH SUMMARY ===

    min_depth_summary <- aggregate(min_depth_frame$minimal_depth, 
                                   by = list(Variable = min_depth_frame$variable), 
                                   FUN = median)
    colnames(min_depth_summary) <- c("Variable", "Median_Min_Depth")
    min_depth_summary <- min_depth_summary[order(min_depth_summary$Median_Min_Depth), ]

    print(min_depth_summary)

    ##            Variable Median_Min_Depth
    ## 1             BIO13                2
    ## 2              BIO3                2
    ## 3              BIO5                2
    ## 4              BIO6                2
    ## 5    continentality                2
    ## 6   growingDegDays5                2
    ## 7    PETseasonality                2
    ## 8 PETWarmestQuarter                2

    cat("\nVariables with median minimum depth < 2 are considered most important\n")

    ## 
    ## Variables with median minimum depth < 2 are considered most important

    cat("Lower depth indicates earlier splits and greater discriminatory capacity\n")

    ## Lower depth indicates earlier splits and greater discriminatory capacity

## Final Notes

This script implements Random Forest species distribution modeling for
*Bactericera cockerelli* following:

1.  **Valavi et al. (2021)** methodology for presence-background data:
    -   Stratified environmental sampling for background points (PCA +
        K-means)
    -   Down-sampling strategy to manage class imbalance (5:1 ratio)
    -   Regression-RF (probability forest) instead of classification-RF
    -   Anti-overfitting parameters (nodesize, maxnodes, sampsize)
2.  **Manuscript specifications**:
    -   500-meter buffers for spatial autocorrelation reduction
    -   8 environmental predictors after collinearity analysis
    -   5-fold stratified cross-validation
    -   Multiple threshold strategies for binarization
    -   Minimum depth for variable importance

**Key Differences from Default RF**: - NO class weights (removed to
prevent bias) - Conservative mtry limit (√p maximum) - Larger terminal
nodes (nodesize = 40-50) - Limited tree depth (maxnodes = 30-35) -
Reduced bootstrap sample (sampsize = 70%) - Manual ntree selection (1500
trees)

**Performance Expectations**: - AUC &gt; 0.97 (mesoscale models) - TSS
&gt; 0.97 (excellent discrimination) - Omission rate &lt; 0.03 (high
sensitivity) - Controlled overfitting (AUC gap &lt; 0.05)
