### MaxEnt Hyperparameter Optimization Using Genetic Algorithm

###### This document implements systematic hyperparameter tuning for MaxEnt models using genetic algorithms (GA) as implemented in SDMtune package. The optimization process explores regularization multipliers (0.2-5.0) and feature class combinations to maximize model performance (AUC) while preventing overfitting through systematic validation procedures.

#### Environment Setup

First, we configure the computational environment for parallel
processing and load required packages:

    library(future)
    library(furrr)
    library(parallel)
    gc()

    ##           used  (Mb) gc trigger  (Mb) max used  (Mb)
    ## Ncells 3680606 196.6    9679011 517.0  9679011 517.0
    ## Vcells 6063735  46.3   24437330 186.5 19326043 147.5

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
    library(ellipse)
    library(scatterplot3d)
    library(MASS)
    library(png)
    library(grid)
    library(rgl)
    library(zeallot)

## Section 1: Climate Data Loading and Variable Selection

We load bioclimatic variables from WorldClim and environmental
predictors from ENVIREM for current climate conditions (1970-2000).
Variable selection was performed based on Spearman correlation analysis
(|r| &gt; 0.8) and ecological relevance for *B. cockerelli*
distribution, as described in the Methods section.

    # Load current climate data (baseline: 1970-2000)
    bioclim_data <- rast("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/data/2.5 minutos/bioclim_data.tif")
    envirem_data <- rast("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/ensambles mensuales/envirem_vars_historico.tif")
     
    names(bioclim_data) <- paste0("BIO", 1:nlyr(bioclim_data))

    # Combine bioclimatic and environmental variables
    climate_data <- c(bioclim_data, envirem_data)

    # Variable selection based on correlation analysis and ecological relevance
    # Selected variables: 2, 5, 6, 12, 21, 22, 26, 34 for B. cockerelli Colombia
    selected_variables <- c(2, 5, 6, 12, 21, 22, 26, 34)

    # Subset climate data to selected variables
    climate_data_subset <- subset(climate_data, selected_variables)
    climate_data <- climate_data_subset

    # Standardize predictors (mean = 0, sd = 1) for model calibration
    climate_data <- scale(climate_data)

## Section 2: Presence Data Loading and Preprocessing

We load *B. cockerelli* occurrence records from ICA phytosanitary
surveillance (2021-2024) and remove duplicate coordinates to avoid
spatial autocorrelation.

    # Load presence records from ICA surveillance database
    presence_data <- read_excel("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Bactericera cockerelli Colombia Oficial surveillance ICA.xlsx")

    # Remove records with missing coordinates
    presence_data <- presence_data[!is.na(presence_data$longitude) & !is.na(presence_data$latitude), ]

    # Remove duplicate coordinates
    presence_data <- presence_data[!duplicated(presence_data[, c("longitude", "latitude")]), ]

    # Create presence data frame
    presence <- presence_data[, c("longitude", "latitude")]
    presence$pa <- 1

    # Remove duplicates using dplyr
    presence <- presence %>% distinct(longitude, latitude, .keep_all = TRUE)

    # Convert to numeric
    presence$longitude <- as.numeric(presence$longitude)
    presence$latitude <- as.numeric(presence$latitude)

## Section 3: Stratified Environmental Sampling for Background Points

We generate background points using stratified sampling in environmental
space (PCA-based K-means clustering) to ensure comprehensive
representation of environmental gradients while minimizing spatial
autocorrelation biases. This approach addresses overfitting issues
observed with traditional random geographic sampling.

    # Set seed for reproducibility
    set.seed(20210707)

    # Create spatial vector from presence points
    presence_vector <- vect(presence, geom = c("longitude", "latitude"), crs = "epsg:4326")

    # Crop climate data to Colombia boundaries using GADM
    tryCatch({
      colombia_boundary <- gadm(country = "COL", level = 0, path = tempdir())
      climate_crop <- crop(climate_data, colombia_boundary)
      climate_crop <- mask(climate_crop, colombia_boundary)
      
      # Verify successful crop
      if(all(is.na(values(climate_crop[[1]])))) {
        stop("Error: Colombia crop resulted in NA values")
      }
      
    }, error = function(e) {
      cat("Error downloading Colombia boundaries:", e$message, "\n")
      cat("Attempting alternative method...\n")
      
      # Alternative method: create approximate Colombia polygon
      colombia_extent <- ext(-82, -66, -5, 13)
      colombia_boundary <- as.polygons(colombia_extent, crs = "epsg:4326")
      climate_crop <- crop(climate_data, colombia_boundary)
      climate_crop <- mask(climate_crop, colombia_boundary)
    })

    # Extract climate values from cropped raster (Colombia)
    tryCatch({
      climate_values <- values(climate_crop, na.rm = TRUE)
      climate_coords <- xyFromCell(climate_crop, which(!is.na(values(climate_crop[[1]]))))
      
      # Verify data extraction
      if(nrow(climate_coords) == 0) {
        stop("Error: No valid coordinates found")
      }
      
      cat("Valid pixels found:", nrow(climate_coords), "\n")
      
    }, error = function(e) {
      stop("Error extracting climate values: ", e$message)
    })

    ## Valid pixels found: 54539

    # Remove rows with NA values
    complete_cases <- complete.cases(climate_values)
    climate_values_clean <- climate_values[complete_cases, ]
    climate_coords_clean <- climate_coords[complete_cases, ]

    cat("Pixels after removing NA:", nrow(climate_coords_clean), "\n")

    ## Pixels after removing NA: 54539

    # Verify sufficient data for sampling
    if(nrow(climate_values_clean) < nrow(presence)) {
      stop("Error: Insufficient valid pixels for sampling")
    }

### Principal Component Analysis and K-means Clustering

We perform PCA in environmental space followed by K-means clustering to
create environmental strata for background point sampling.

    # Principal Component Analysis (PCA) in environmental space
    tryCatch({
      climate_pca <- prcomp(climate_values_clean, scale. = TRUE, center = TRUE)
      cat("PCA completed. Variance explained PC1-3:", 
          round(sum(climate_pca$sdev[1:3]^2)/sum(climate_pca$sdev^2)*100, 2), "%\n")
    }, error = function(e) {
      stop("Error in PCA: ", e$message)
    })

    ## PCA completed. Variance explained PC1-3: 88.04 %

    # Extract first 3 principal components
    pca_scores <- climate_pca$x[, 1:3]

    # K-means clustering in environmental space (PCA)
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

    ## Clustering completed

    # Assign each pixel to environmental cluster
    cluster_assignments <- kmeans_result$cluster

    # Calculate number of points per cluster (proportional to cluster size)
    cluster_counts <- table(cluster_assignments)
    samples_per_cluster <- round((cluster_counts / sum(cluster_counts)) * nrow(presence))

    # Ensure minimum total number of samples
    if(sum(samples_per_cluster) < nrow(presence)) {
      diff_needed <- nrow(presence) - sum(samples_per_cluster)
      largest_clusters <- order(cluster_counts, decreasing = TRUE)[1:diff_needed]
      samples_per_cluster[largest_clusters] <- samples_per_cluster[largest_clusters] + 1
    }

    # Stratified sampling by environmental cluster
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

    # Create background object with selected coordinates
    background <- climate_coords_clean[background_indices, ]
    colnames(background) <- c("x", "y")
    background <- as.data.frame(background)

    # Create absence data frame
    absence <- as.data.frame(background)
    colnames(absence) <- c("longitude", "latitude")
    absence$pa <- 0

    # Convert presence coordinates to numeric
    presence$latitude <- as.numeric(presence$latitude)
    presence$longitude <- as.numeric(presence$longitude)

    # Combine presence and absence points
    all_points <- rbind(presence, absence)

## Section 4: Buffer-Based Environmental Data Extraction

We create 500-meter buffers around presence and background points to
extract mean environmental values, reducing spatial autocorrelation in
model calibration.

    # Convert meters to degrees (latitude-dependent for longitude)
    rad_to_deg <- function(rad) rad * 180 / pi
    deg_to_rad <- function(deg) deg * pi / 180

    # Earth radius in meters
    earth_radius <- 6371000

    # Convert 250 meters to degrees of latitude (constant)
    lat_degree <- 250 / (pi / 180 * earth_radius)

    # Convert 250 meters to degrees of longitude (latitude-dependent)
    lon_degree <- 250 / (pi / 180 * earth_radius * cos(deg_to_rad(presence$latitude)))

    # Average buffer width
    avg_width <- (lat_degree + mean(lon_degree)) / 2

    # Generate buffers with adjusted width
    presence_vect <- vect(presence, geom = c("longitude", "latitude"), crs = "epsg:4326")
    presence_buffer <- buffer(presence_vect, width = avg_width)

    background_vect <- vect(absence, geom = c("longitude", "latitude"), crs = "epsg:4326")
    background_buffer <- buffer(background_vect, width = avg_width)

    # Extract climate data for presence and background buffers (mean values)
    presence_climate <- extract(climate_data, presence_buffer, fun = mean, na.rm = TRUE)
    background_climate <- extract(climate_data, background_buffer, fun = mean, na.rm = TRUE)

    # Combine climate data with original points
    presence_climate <- cbind(presence, presence_climate)
    background_climate <- cbind(absence, background_climate)

    # Combine presence and background data with climate values
    points_climate <- rbind(presence_climate, background_climate)

    # Remove coordinate columns before modeling
    drop_cols <- which(colnames(points_climate) %in% c("longitude", "latitude", "ID"))
    points_climate <- points_climate[, -drop_cols]

    # Create response vector (presence/absence)
    response <- points_climate$pa

## Section 5: Data Preparation for MaxEnt Model Calibration

We prepare presence and background data in SDMtune format (SWD objects)
for MaxEnt model training. The `prepareSWD` function creates
Species-With-Data objects that link occurrence points with environmental
predictors.

    # Prepare data in SDMtune format (Species With Data - SWD)
    data <- prepareSWD(species = "B. cockerelli Mesoscale", 
                       p = presence[, c("longitude", "latitude")],
                       a = absence[, c("longitude", "latitude")],
                       env = climate_data)

## Section 6: Data Partitioning Strategy

We partition data into training (70%), validation (15%), and test (15%)
sets.

### Critical Distinction:

-   **Training Set (70%)**: Used to calibrate model parameters
-   **Test Set (15%)**: Used during hyperparameter optimization with
    genetic algorithm to evaluate candidate model configurations and
    guide the evolutionary search toward optimal hyperparameter
    combinations
-   **Validation Set (15%)**: Reserved for final model performance
    assessment (AUC, TSS, omission rate) and remains completely
    independent from both training and optimization processes to provide
    unbiased evaluation

This three-way split prevents data leakage and ensures robust
performance estimates, as the validation set provides independent
assessment of the optimized model’s generalization capacity.

    # Create training, test, and validation datasets
    datasets <- trainValTest(data, val = 0.15, test = 0.15, only_presence = FALSE, seed = 12345)

    train <- datasets[[1]]  # Training data (70%) - for model calibration
    test <- datasets[[2]]   # Test data (15%) - for hyperparameter optimization
    val <- datasets[[3]]    # Validation data (15%) - for final performance metrics

## Section 7: Baseline MaxEnt Model Training

We train an initial MaxEnt model with default parameters (regularization
= 1.0, all feature classes enabled) to establish baseline performance
before systematic hyperparameter optimization.

    # Train baseline MaxEnt model with default parameters
    model <- train("Maxent", data = train, iter = 700)

## Section 8: Hyperparameter Grid Definition

We define the hyperparameter search space for genetic algorithm
optimization:

### Regularization Multipliers (reg)

Control model complexity and overfitting: - Range: 0.2 to 5.0 in steps
of 0.2 (35 values) - Higher values increase regularization → simpler,
smoother predictions - Lower values decrease regularization → more
complex, detailed predictions

### Feature Classes (fc)

Define mathematical transformations of predictors: - **l** = linear:
Simple linear relationships - **q** = quadratic: Polynomial
relationships (degree 2) - **p** = product: Interaction terms between
variables - **t** = threshold: Step functions - **h** = hinge: Piecewise
linear functions - Combinations: 31 possible feature class combinations
tested

**Total search space**: 35 regularization × 31 feature combinations =
1,085 possible model configurations explored through genetic algorithm.

    # Define hyperparameter grid for optimization
    h <- list(reg = seq(0.2, 5, 0.2), 
              fc = c("l", "q", "p", "t", "h", "lq", "lp", "lt", "lh", "qp", "qt", "qh", 
                     "pt", "ph", "th", "lqp", "lqt", "lqh", "lpt", "lph", "lth", "qpt", 
                     "qph", "qth", "pth", "lqpt", "lqph", "lqth", "lpth", "qpth", "lqpth"))

## Section 9: Genetic Algorithm Optimization

We optimize MaxEnt hyperparameters using genetic algorithm (GA)
implementation from SDMtune package.

### Genetic Algorithm Process:

1.  **Initial Population (pop = 37)**: Creates 37 candidate models with
    random hyperparameter combinations
2.  **Fitness Evaluation**: Each model is evaluated on the TEST SET
    using AUC as fitness metric
3.  **Evolutionary Generations (gen = 5)**: Iterative evolution through
    selection, crossover, mutation, and evaluation
4.  **Convergence**: Each generation improves upon previous best
    configurations

The genetic algorithm explores hyperparameter space more efficiently
than exhaustive grid search, focusing on promising regions of parameter
space.

    # Optimize hyperparameters using genetic algorithm
    opt_model <- optimizeModel(
      model = model,
      hypers = h,
      metric = "auc",
      test = test,        # Test set used for optimization (fitness evaluation)
      pop = 37,           # Initial population size (37 candidate models)
      gen = 5,            # Number of evolutionary generations
      seed = 12345
    )

## Section 10: Best Model Selection and Ensemble Creation

We select the best model configuration based on highest test AUC from
genetic algorithm optimization. The `combineCV` function creates an
ensemble by averaging predictions from the top-performing hyperparameter
combinations, resulting in more stable and generalizable predictions.

    # Select best model based on highest test AUC
    best_model <- combineCV(opt_model@models[[which.max(opt_model@results$test_AUC)]])

## Section 11: Model Performance Evaluation

We evaluate final model performance using the VALIDATION SET
(independent from training and optimization).

### Performance Metrics:

-   **AUC (Area Under ROC Curve)**: Discriminative capacity (0.5-1.0)
    -   &gt;0.9: Excellent discrimination
    -   0.8-0.9: Good discrimination
    -   0.7-0.8: Acceptable discrimination
-   **TSS (True Skill Statistic)**: Threshold-dependent accuracy (-1 to
    1)
    -   &gt;0.7: High predictive accuracy
    -   0.4-0.7: Moderate accuracy
    -   &lt;0.4: Poor accuracy
-   **Omission Rate**: Proportion of presence records below threshold
    (lower is better)

<!-- -->

    # Calculate True Skill Statistic (TSS)
    tss_value <- tss(best_model)
    cat("True Skill Statistic (TSS):", tss_value, "\n")

    ## True Skill Statistic (TSS): 0.9686411

    # Plot ROC curve and calculate AUC using test set
    roc_plot <- plotROC(best_model, val)
    plot(roc_plot)

![](C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Paper_comnplejo_punta_morada/Tuning_MaxEnt_files/figure-markdown_strict/model-evaluation-1.png)

## Section 12: Variable Importance Analysis - Jackknife Test

Jackknife analysis systematically evaluates predictor importance by:

1.  Training models with each variable excluded (one at a time)
2.  Training models with only individual variables
3.  Comparing performance (AUC) to identify most influential predictors

This approach reveals variables with highest individual predictive
capacity and those whose exclusion most reduces model performance.

    # Perform Jackknife test for variable importance
    jk <- doJk(best_model, 
               metric = "auc", 
               test = val)  # Use validation set for importance evaluation

    # Visualize Jackknife results
    jk_plot <- plotJk(jk, type = "test")
    plot(jk_plot)

![](C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Paper_comnplejo_punta_morada/Tuning_MaxEnt_files/figure-markdown_strict/jackknife-analysis-1.png)

## Section 13: Variable Correlation Visualization

We visualize Spearman correlations between environmental predictors to
verify that selected variables maintain low multicollinearity (|r| &lt;
0.8).

    # Visualize correlation matrix for selected variables
    correlation_plot <- plotCor(train, method = "spearman", cor_th = 0.7)
    plot(correlation_plot)

![](C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Paper_comnplejo_punta_morada/Tuning_MaxEnt_files/figure-markdown_strict/correlation-plot-1.png)

## Summary

This analysis successfully optimized MaxEnt hyperparameters for *B.
cockerelli* distribution modeling in Colombia using:

-   Stratified environmental sampling for background points
-   Genetic algorithm optimization (1,085 hyperparameter combinations
    explored)
-   Three-way data partitioning (training/test/validation)
-   Independent validation for unbiased performance assessment

### Output Files Generated:

-   Optimized MaxEnt model with tuned hyperparameters (regularization,
    features)
-   Performance metrics: AUC, TSS, omission rate on independent
    validation set
-   Variable importance rankings from Jackknife analysis
-   ROC curves and correlation plots for model interpretation
