## Introduction

This document implements a weighted ensemble approach to reduce
uncertainty in climate projections from 10 CMIP6 General Circulation
Models (GCMs). Each GCM contains 19 annual mean bioclimatic layers
(WorldClim variables) for the 2021-2040 period under SSP2-4.5 and
SSP5-8.5 scenarios.

Ensemble weights are calculated using a combined metric: correlation² /
normalized RMSE, where models with higher correlation and lower
prediction error receive greater weight. This approach has been shown to
improve projection accuracy compared to simple model averaging.

**Key Processing Steps:**

1.  **Crop to Colombia** - Reduce global rasters to national extent to
    improve computational efficiency
2.  **Calculate Performance Metrics** - Assess each GCM against
    historical baseline (1970-2000)
3.  **Weight Assignment** - Combine correlation and RMSE into unified
    weighting scheme
4.  **Ensemble Generation** - Produce weighted mean projections for each
    bioclimatic variable
5.  **Anomaly Calculation** - Quantify future deviations from historical
    baseline

------------------------------------------------------------------------

## Load Required Libraries

Import essential packages for raster processing and parallel
computation.

    library(terra)
    library(future)
    library(furrr)
    library(parallel)

## Configure High-Performance Computing Environment

Optimize memory allocation and parallel processing for intensive raster
operations. Since these are global climate layers, processing can take
days without proper optimization.

    # Clear memory and enable terra optimization flags
    gc()

    ##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells  6029705 322.1   14142814  755.4  14142814  755.4
    ## Vcells 59379597 453.1  165414027 1262.1 165414027 1262.1

    options(terra.keys = TRUE)

    # Configure terra for intensive RAM and CPU usage
    terraOptions(
      memfrac  = 0.9,                          # Use up to 90% of available RAM
      progress = FALSE,                        # Suppress progress bars for cleaner output
      threads  = detectCores(logical = FALSE) - 1  # Reserve 1 core for system operations
    )

    # Configure parallel processing with future/furrr
    num_cores <- detectCores(logical = FALSE) - 1  # Leave 1 core free for OS
    plan(multisession, workers = num_cores)
    options(future.globals.maxSize = 15 * 1024^3)   # Allow up to ~15 GB for large objects

    # Set furrr options for heavy computational tasks
    furrr_options(
      seed       = TRUE,
      scheduling = Inf     # Assign heavy tasks immediately to available cores
    )

    ## <furrr_options>

## Load Historical Baseline Climate Data

Import historical bioclimatic variables (1970-2000) from WorldClim v2.1
at 2.5 arc-minute resolution. This dataset contains 19 annual mean
layers serving as the calibration baseline for GCM performance
evaluation.

    # Load historical bioclimatic baseline (19 layers)
    bio_data <- rast("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/data/2.5 minutos/bioclim_data.tif")

## Load Future Climate Projections (SSP2-4.5)

Import 10 CMIP6 GCMs for SSP2-4.5 scenario (moderate mitigation). Each
model contains 19 bioclimatic layers representing annual means for
2021-2040.

    # Load 10 GCM projections for SSP2-4.5 (190 layers total: 10 models × 19 variables)
    models_bio <- paste0("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/BIOmodelos 2021 - 2040/245/biomodel", 1:10, ".tif")
    models <- rast(models_bio)
    names(models)

    ##   [1] "wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_1"    
    ##   [2] "wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_2"    
    ##   [3] "wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_3"    
    ##   [4] "wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_4"    
    ##   [5] "wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_5"    
    ##   [6] "wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_6"    
    ##   [7] "wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_7"    
    ##   [8] "wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_8"    
    ##   [9] "wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_9"    
    ##  [10] "wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_10"   
    ##  [11] "wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_11"   
    ##  [12] "wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_12"   
    ##  [13] "wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_13"   
    ##  [14] "wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_14"   
    ##  [15] "wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_15"   
    ##  [16] "wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_16"   
    ##  [17] "wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_17"   
    ##  [18] "wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_18"   
    ##  [19] "wc2.1_2.5m_bioc_ACCESS-CM2_ssp245_2021-2040_19"   
    ##  [20] "wc2.1_2.5m_bioc_CMCC-ESM2_ssp245_2021-2040_1"     
    ##  [21] "wc2.1_2.5m_bioc_CMCC-ESM2_ssp245_2021-2040_2"     
    ##  [22] "wc2.1_2.5m_bioc_CMCC-ESM2_ssp245_2021-2040_3"     
    ##  [23] "wc2.1_2.5m_bioc_CMCC-ESM2_ssp245_2021-2040_4"     
    ##  [24] "wc2.1_2.5m_bioc_CMCC-ESM2_ssp245_2021-2040_5"     
    ##  [25] "wc2.1_2.5m_bioc_CMCC-ESM2_ssp245_2021-2040_6"     
    ##  [26] "wc2.1_2.5m_bioc_CMCC-ESM2_ssp245_2021-2040_7"     
    ##  [27] "wc2.1_2.5m_bioc_CMCC-ESM2_ssp245_2021-2040_8"     
    ##  [28] "wc2.1_2.5m_bioc_CMCC-ESM2_ssp245_2021-2040_9"     
    ##  [29] "wc2.1_2.5m_bioc_CMCC-ESM2_ssp245_2021-2040_10"    
    ##  [30] "wc2.1_2.5m_bioc_CMCC-ESM2_ssp245_2021-2040_11"    
    ##  [31] "wc2.1_2.5m_bioc_CMCC-ESM2_ssp245_2021-2040_12"    
    ##  [32] "wc2.1_2.5m_bioc_CMCC-ESM2_ssp245_2021-2040_13"    
    ##  [33] "wc2.1_2.5m_bioc_CMCC-ESM2_ssp245_2021-2040_14"    
    ##  [34] "wc2.1_2.5m_bioc_CMCC-ESM2_ssp245_2021-2040_15"    
    ##  [35] "wc2.1_2.5m_bioc_CMCC-ESM2_ssp245_2021-2040_16"    
    ##  [36] "wc2.1_2.5m_bioc_CMCC-ESM2_ssp245_2021-2040_17"    
    ##  [37] "wc2.1_2.5m_bioc_CMCC-ESM2_ssp245_2021-2040_18"    
    ##  [38] "wc2.1_2.5m_bioc_CMCC-ESM2_ssp245_2021-2040_19"    
    ##  [39] "wc2.1_2.5m_bioc_EC-Earth3-Veg_ssp245_2021-2040_1" 
    ##  [40] "wc2.1_2.5m_bioc_EC-Earth3-Veg_ssp245_2021-2040_2" 
    ##  [41] "wc2.1_2.5m_bioc_EC-Earth3-Veg_ssp245_2021-2040_3" 
    ##  [42] "wc2.1_2.5m_bioc_EC-Earth3-Veg_ssp245_2021-2040_4" 
    ##  [43] "wc2.1_2.5m_bioc_EC-Earth3-Veg_ssp245_2021-2040_5" 
    ##  [44] "wc2.1_2.5m_bioc_EC-Earth3-Veg_ssp245_2021-2040_6" 
    ##  [45] "wc2.1_2.5m_bioc_EC-Earth3-Veg_ssp245_2021-2040_7" 
    ##  [46] "wc2.1_2.5m_bioc_EC-Earth3-Veg_ssp245_2021-2040_8" 
    ##  [47] "wc2.1_2.5m_bioc_EC-Earth3-Veg_ssp245_2021-2040_9" 
    ##  [48] "wc2.1_2.5m_bioc_EC-Earth3-Veg_ssp245_2021-2040_10"
    ##  [49] "wc2.1_2.5m_bioc_EC-Earth3-Veg_ssp245_2021-2040_11"
    ##  [50] "wc2.1_2.5m_bioc_EC-Earth3-Veg_ssp245_2021-2040_12"
    ##  [51] "wc2.1_2.5m_bioc_EC-Earth3-Veg_ssp245_2021-2040_13"
    ##  [52] "wc2.1_2.5m_bioc_EC-Earth3-Veg_ssp245_2021-2040_14"
    ##  [53] "wc2.1_2.5m_bioc_EC-Earth3-Veg_ssp245_2021-2040_15"
    ##  [54] "wc2.1_2.5m_bioc_EC-Earth3-Veg_ssp245_2021-2040_16"
    ##  [55] "wc2.1_2.5m_bioc_EC-Earth3-Veg_ssp245_2021-2040_17"
    ##  [56] "wc2.1_2.5m_bioc_EC-Earth3-Veg_ssp245_2021-2040_18"
    ##  [57] "wc2.1_2.5m_bioc_EC-Earth3-Veg_ssp245_2021-2040_19"
    ##  [58] "wc2.1_2.5m_bioc_GISS-E2-1-G_ssp245_2021-2040_1"   
    ##  [59] "wc2.1_2.5m_bioc_GISS-E2-1-G_ssp245_2021-2040_2"   
    ##  [60] "wc2.1_2.5m_bioc_GISS-E2-1-G_ssp245_2021-2040_3"   
    ##  [61] "wc2.1_2.5m_bioc_GISS-E2-1-G_ssp245_2021-2040_4"   
    ##  [62] "wc2.1_2.5m_bioc_GISS-E2-1-G_ssp245_2021-2040_5"   
    ##  [63] "wc2.1_2.5m_bioc_GISS-E2-1-G_ssp245_2021-2040_6"   
    ##  [64] "wc2.1_2.5m_bioc_GISS-E2-1-G_ssp245_2021-2040_7"   
    ##  [65] "wc2.1_2.5m_bioc_GISS-E2-1-G_ssp245_2021-2040_8"   
    ##  [66] "wc2.1_2.5m_bioc_GISS-E2-1-G_ssp245_2021-2040_9"   
    ##  [67] "wc2.1_2.5m_bioc_GISS-E2-1-G_ssp245_2021-2040_10"  
    ##  [68] "wc2.1_2.5m_bioc_GISS-E2-1-G_ssp245_2021-2040_11"  
    ##  [69] "wc2.1_2.5m_bioc_GISS-E2-1-G_ssp245_2021-2040_12"  
    ##  [70] "wc2.1_2.5m_bioc_GISS-E2-1-G_ssp245_2021-2040_13"  
    ##  [71] "wc2.1_2.5m_bioc_GISS-E2-1-G_ssp245_2021-2040_14"  
    ##  [72] "wc2.1_2.5m_bioc_GISS-E2-1-G_ssp245_2021-2040_15"  
    ##  [73] "wc2.1_2.5m_bioc_GISS-E2-1-G_ssp245_2021-2040_16"  
    ##  [74] "wc2.1_2.5m_bioc_GISS-E2-1-G_ssp245_2021-2040_17"  
    ##  [75] "wc2.1_2.5m_bioc_GISS-E2-1-G_ssp245_2021-2040_18"  
    ##  [76] "wc2.1_2.5m_bioc_GISS-E2-1-G_ssp245_2021-2040_19"  
    ##  [77] "wc2.1_2.5m_bioc_INM-CM5-0_ssp245_2021-2040_1"     
    ##  [78] "wc2.1_2.5m_bioc_INM-CM5-0_ssp245_2021-2040_2"     
    ##  [79] "wc2.1_2.5m_bioc_INM-CM5-0_ssp245_2021-2040_3"     
    ##  [80] "wc2.1_2.5m_bioc_INM-CM5-0_ssp245_2021-2040_4"     
    ##  [81] "wc2.1_2.5m_bioc_INM-CM5-0_ssp245_2021-2040_5"     
    ##  [82] "wc2.1_2.5m_bioc_INM-CM5-0_ssp245_2021-2040_6"     
    ##  [83] "wc2.1_2.5m_bioc_INM-CM5-0_ssp245_2021-2040_7"     
    ##  [84] "wc2.1_2.5m_bioc_INM-CM5-0_ssp245_2021-2040_8"     
    ##  [85] "wc2.1_2.5m_bioc_INM-CM5-0_ssp245_2021-2040_9"     
    ##  [86] "wc2.1_2.5m_bioc_INM-CM5-0_ssp245_2021-2040_10"    
    ##  [87] "wc2.1_2.5m_bioc_INM-CM5-0_ssp245_2021-2040_11"    
    ##  [88] "wc2.1_2.5m_bioc_INM-CM5-0_ssp245_2021-2040_12"    
    ##  [89] "wc2.1_2.5m_bioc_INM-CM5-0_ssp245_2021-2040_13"    
    ##  [90] "wc2.1_2.5m_bioc_INM-CM5-0_ssp245_2021-2040_14"    
    ##  [91] "wc2.1_2.5m_bioc_INM-CM5-0_ssp245_2021-2040_15"    
    ##  [92] "wc2.1_2.5m_bioc_INM-CM5-0_ssp245_2021-2040_16"    
    ##  [93] "wc2.1_2.5m_bioc_INM-CM5-0_ssp245_2021-2040_17"    
    ##  [94] "wc2.1_2.5m_bioc_INM-CM5-0_ssp245_2021-2040_18"    
    ##  [95] "wc2.1_2.5m_bioc_INM-CM5-0_ssp245_2021-2040_19"    
    ##  [96] "wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2021-2040_1"  
    ##  [97] "wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2021-2040_2"  
    ##  [98] "wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2021-2040_3"  
    ##  [99] "wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2021-2040_4"  
    ## [100] "wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2021-2040_5"  
    ## [101] "wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2021-2040_6"  
    ## [102] "wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2021-2040_7"  
    ## [103] "wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2021-2040_8"  
    ## [104] "wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2021-2040_9"  
    ## [105] "wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2021-2040_10" 
    ## [106] "wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2021-2040_11" 
    ## [107] "wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2021-2040_12" 
    ## [108] "wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2021-2040_13" 
    ## [109] "wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2021-2040_14" 
    ## [110] "wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2021-2040_15" 
    ## [111] "wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2021-2040_16" 
    ## [112] "wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2021-2040_17" 
    ## [113] "wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2021-2040_18" 
    ## [114] "wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2021-2040_19" 
    ## [115] "wc2.1_2.5m_bioc_MIROC6_ssp245_2041-2060_1"        
    ## [116] "wc2.1_2.5m_bioc_MIROC6_ssp245_2041-2060_2"        
    ## [117] "wc2.1_2.5m_bioc_MIROC6_ssp245_2041-2060_3"        
    ## [118] "wc2.1_2.5m_bioc_MIROC6_ssp245_2041-2060_4"        
    ## [119] "wc2.1_2.5m_bioc_MIROC6_ssp245_2041-2060_5"        
    ## [120] "wc2.1_2.5m_bioc_MIROC6_ssp245_2041-2060_6"        
    ## [121] "wc2.1_2.5m_bioc_MIROC6_ssp245_2041-2060_7"        
    ## [122] "wc2.1_2.5m_bioc_MIROC6_ssp245_2041-2060_8"        
    ## [123] "wc2.1_2.5m_bioc_MIROC6_ssp245_2041-2060_9"        
    ## [124] "wc2.1_2.5m_bioc_MIROC6_ssp245_2041-2060_10"       
    ## [125] "wc2.1_2.5m_bioc_MIROC6_ssp245_2041-2060_11"       
    ## [126] "wc2.1_2.5m_bioc_MIROC6_ssp245_2041-2060_12"       
    ## [127] "wc2.1_2.5m_bioc_MIROC6_ssp245_2041-2060_13"       
    ## [128] "wc2.1_2.5m_bioc_MIROC6_ssp245_2041-2060_14"       
    ## [129] "wc2.1_2.5m_bioc_MIROC6_ssp245_2041-2060_15"       
    ## [130] "wc2.1_2.5m_bioc_MIROC6_ssp245_2041-2060_16"       
    ## [131] "wc2.1_2.5m_bioc_MIROC6_ssp245_2041-2060_17"       
    ## [132] "wc2.1_2.5m_bioc_MIROC6_ssp245_2041-2060_18"       
    ## [133] "wc2.1_2.5m_bioc_MIROC6_ssp245_2041-2060_19"       
    ## [134] "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060_1" 
    ## [135] "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060_2" 
    ## [136] "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060_3" 
    ## [137] "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060_4" 
    ## [138] "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060_5" 
    ## [139] "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060_6" 
    ## [140] "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060_7" 
    ## [141] "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060_8" 
    ## [142] "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060_9" 
    ## [143] "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060_10"
    ## [144] "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060_11"
    ## [145] "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060_12"
    ## [146] "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060_13"
    ## [147] "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060_14"
    ## [148] "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060_15"
    ## [149] "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060_16"
    ## [150] "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060_17"
    ## [151] "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060_18"
    ## [152] "wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2041-2060_19"
    ## [153] "wc2.1_2.5m_bioc_MRI-ESM2-0_ssp245_2021-2040_1"    
    ## [154] "wc2.1_2.5m_bioc_MRI-ESM2-0_ssp245_2021-2040_2"    
    ## [155] "wc2.1_2.5m_bioc_MRI-ESM2-0_ssp245_2021-2040_3"    
    ## [156] "wc2.1_2.5m_bioc_MRI-ESM2-0_ssp245_2021-2040_4"    
    ## [157] "wc2.1_2.5m_bioc_MRI-ESM2-0_ssp245_2021-2040_5"    
    ## [158] "wc2.1_2.5m_bioc_MRI-ESM2-0_ssp245_2021-2040_6"    
    ## [159] "wc2.1_2.5m_bioc_MRI-ESM2-0_ssp245_2021-2040_7"    
    ## [160] "wc2.1_2.5m_bioc_MRI-ESM2-0_ssp245_2021-2040_8"    
    ## [161] "wc2.1_2.5m_bioc_MRI-ESM2-0_ssp245_2021-2040_9"    
    ## [162] "wc2.1_2.5m_bioc_MRI-ESM2-0_ssp245_2021-2040_10"   
    ## [163] "wc2.1_2.5m_bioc_MRI-ESM2-0_ssp245_2021-2040_11"   
    ## [164] "wc2.1_2.5m_bioc_MRI-ESM2-0_ssp245_2021-2040_12"   
    ## [165] "wc2.1_2.5m_bioc_MRI-ESM2-0_ssp245_2021-2040_13"   
    ## [166] "wc2.1_2.5m_bioc_MRI-ESM2-0_ssp245_2021-2040_14"   
    ## [167] "wc2.1_2.5m_bioc_MRI-ESM2-0_ssp245_2021-2040_15"   
    ## [168] "wc2.1_2.5m_bioc_MRI-ESM2-0_ssp245_2021-2040_16"   
    ## [169] "wc2.1_2.5m_bioc_MRI-ESM2-0_ssp245_2021-2040_17"   
    ## [170] "wc2.1_2.5m_bioc_MRI-ESM2-0_ssp245_2021-2040_18"   
    ## [171] "wc2.1_2.5m_bioc_MRI-ESM2-0_ssp245_2021-2040_19"   
    ## [172] "wc2.1_2.5m_bioc_UKESM1-0-LL_ssp245_2021-2040_1"   
    ## [173] "wc2.1_2.5m_bioc_UKESM1-0-LL_ssp245_2021-2040_2"   
    ## [174] "wc2.1_2.5m_bioc_UKESM1-0-LL_ssp245_2021-2040_3"   
    ## [175] "wc2.1_2.5m_bioc_UKESM1-0-LL_ssp245_2021-2040_4"   
    ## [176] "wc2.1_2.5m_bioc_UKESM1-0-LL_ssp245_2021-2040_5"   
    ## [177] "wc2.1_2.5m_bioc_UKESM1-0-LL_ssp245_2021-2040_6"   
    ## [178] "wc2.1_2.5m_bioc_UKESM1-0-LL_ssp245_2021-2040_7"   
    ## [179] "wc2.1_2.5m_bioc_UKESM1-0-LL_ssp245_2021-2040_8"   
    ## [180] "wc2.1_2.5m_bioc_UKESM1-0-LL_ssp245_2021-2040_9"   
    ## [181] "wc2.1_2.5m_bioc_UKESM1-0-LL_ssp245_2021-2040_10"  
    ## [182] "wc2.1_2.5m_bioc_UKESM1-0-LL_ssp245_2021-2040_11"  
    ## [183] "wc2.1_2.5m_bioc_UKESM1-0-LL_ssp245_2021-2040_12"  
    ## [184] "wc2.1_2.5m_bioc_UKESM1-0-LL_ssp245_2021-2040_13"  
    ## [185] "wc2.1_2.5m_bioc_UKESM1-0-LL_ssp245_2021-2040_14"  
    ## [186] "wc2.1_2.5m_bioc_UKESM1-0-LL_ssp245_2021-2040_15"  
    ## [187] "wc2.1_2.5m_bioc_UKESM1-0-LL_ssp245_2021-2040_16"  
    ## [188] "wc2.1_2.5m_bioc_UKESM1-0-LL_ssp245_2021-2040_17"  
    ## [189] "wc2.1_2.5m_bioc_UKESM1-0-LL_ssp245_2021-2040_18"  
    ## [190] "wc2.1_2.5m_bioc_UKESM1-0-LL_ssp245_2021-2040_19"

## Crop Global Rasters to Colombia

**CRITICAL STEP**: Reduce computational burden by cropping global
climate layers to Colombian extent before ensemble processing. Without
this step, processing would require days of computation.

    # Load Colombian administrative boundaries
    colombia_shp <- vect("C:/Users/manuel.cortes/Desktop/PMP2024/Bactericera/Colombia.shp")

    # Reproject shapefile to match climate data CRS if necessary
    if (!crs(colombia_shp) == crs(bio_data)) {
      colombia_shp <- project(colombia_shp, crs(bio_data))
    }

    if (!crs(colombia_shp) == crs(models)) {
      colombia_shp <- project(colombia_shp, crs(models))
    }

    # Crop and mask historical baseline to Colombia
    cropped <- crop(bio_data, colombia_shp)
    bio_data <- mask(cropped, colombia_shp)

    # Crop and mask GCM projections to Colombia
    cropped_models <- crop(models, colombia_shp)
    models <- mask(cropped_models, colombia_shp)

## Preprocess Climate Layers

Remove areas with missing data to ensure consistent spatial coverage
across all datasets.

    # Mask areas without valid historical data in both datasets
    models <- mask(models, !is.na(bio_data))  
    bio_data <- mask(bio_data, !is.na(bio_data))

## Initialize Storage for Results

Create empty lists to store ensemble outputs and performance metrics for
each bioclimatic variable.

    # Initialize lists for 19 bioclimatic variables
    final_ensemble_list <- vector("list", 19)
    correlation_list <- vector("list", 19)

------------------------------------------------------------------------

## Weighted Ensemble Generation

This loop processes each of the 19 bioclimatic variables independently,
calculating performance-based weights and generating ensemble
projections.

### Algorithm Overview:

For each bioclimatic variable (BIO1-BIO19):

1.  Extract historical layer and corresponding GCM projections
2.  Calculate RMSE (prediction error) and correlation (agreement) for
    each GCM
3.  Normalize RMSE by historical standard deviation for
    scale-independent comparison
4.  Compute weights: (correlation² / normalized RMSE)
5.  Generate weighted mean projection across 10 GCMs

<!-- -->

    # Iterate over 19 bioclimatic variables
    for (var_index in 1:19) {
      
      # Extract historical data and GCM projections for current variable
      hist_data <- bio_data[[var_index]]
      models_var <- models[[seq(var_index, nlyr(models), by = 19)]]  # Extract every 19th layer (10 models)
      
      # Initialize performance metric vectors
      rmse_list <- numeric(nlyr(models_var))
      corr_list <- numeric(nlyr(models_var))
      sd_hist <- sd(values(hist_data), na.rm = TRUE)  # Historical SD for RMSE normalization
      
      # Calculate RMSE and correlation for each GCM
      for (i in 1:nlyr(models_var)) {
        model <- models_var[[i]]
        
        # Handle missing models
        if (is.null(model)) {
          rmse_list[i] <- NA
          corr_list[i] <- NA
          next
        }
        
        # Extract valid pixel values (exclude NA)
        model_values <- values(model)
        hist_values <- values(hist_data)
        valid_indices <- !is.na(model_values) & !is.na(hist_values)
        
        # Calculate RMSE (root mean square error)
        rmse_list[i] <- sqrt(mean((model_values[valid_indices] - hist_values[valid_indices])^2))
        
        # Calculate Pearson correlation coefficient
        corr_list[i] <- cor(model_values[valid_indices], hist_values[valid_indices], method = "pearson")
      }
      
      # Compute ensemble weights using combined metric
      rmse_normalized <- rmse_list / sd_hist  # Scale-independent RMSE
      weights <- (corr_list^2) / rmse_normalized  # Combined performance metric
      weights[is.na(weights)] <- 0  # Handle NA weights
      weights <- weights / sum(weights, na.rm = TRUE)  # Normalize to sum = 1
      
      # Store correlation coefficients for diagnostics
      correlation_list[[var_index]] <- corr_list
      
      # Generate weighted ensemble for current variable
      final_ensemble <- app(models_var, fun = function(x) weighted.mean(x, w = weights, na.rm = TRUE))
      final_ensemble_list[[var_index]] <- final_ensemble
    }

------------------------------------------------------------------------

## Combine Ensemble Layers

Stack all 19 weighted ensemble projections into a single multi-layer
raster object.

    # Combine 19 ensemble layers into single raster stack
    combined_final_ensemble <- rast(final_ensemble_list)

    # Visualize ensemble projections
    plot(combined_final_ensemble, main = "Weighted Ensemble Projections (19 BioClim Variables)")

![](C:\Users\manuel.cortes\Desktop\PMP2024\Bactericera\Paper_comnplejo_punta_morada\GCMs%20ensemble_files/figure-markdown_strict/combine_ensemble-1.png)

## Calculate Climate Anomalies

Compute future deviations from historical baseline by subtracting
historical means from ensemble projections. These anomalies represent
the magnitude of projected climate change for each variable.

    # Calculate anomalies (Future - Historical) for each variable
    combined_anomaly_list <- lapply(1:19, function(i) final_ensemble_list[[i]] - bio_data[[i]])
    combined_anomaly <- rast(combined_anomaly_list)

    # Visualize climate anomalies
    plot(combined_anomaly, main = "Climate Anomalies: 2021-2040 vs 1970-2000 Baseline")

![](C:\Users\manuel.cortes\Desktop\PMP2024\Bactericera\Paper_comnplejo_punta_morada\GCMs%20ensemble_files/figure-markdown_strict/calculate_anomalies-1.png)

------------------------------------------------------------------------

## Summary

This weighted ensemble approach successfully integrated 10 CMIP6 GCMs
into unified climate projections for Colombia. Key advantages of this
methodology:

-   **Uncertainty Reduction**: Performance-based weighting reduces bias
    from poorly-performing GCMs
-   **Computational Efficiency**: Cropping to national extent reduced
    processing time from days to hours
-   **Scale-Independent Metrics**: RMSE normalization ensures fair
    comparison across variables with different units
-   **Transparency**: Correlation coefficients stored for post-hoc model
    evaluation

The resulting ensemble projections serve as inputs for species
distribution modeling under SSP2-4.5 and SSP5-8.5 scenarios, providing
climate-adapted forecasts for *B. cockerelli* establishment potential
through 2040.
