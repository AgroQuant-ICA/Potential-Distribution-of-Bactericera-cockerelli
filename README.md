##🐛 Potential Distribution of Bactericera cockerelli in Colombia

This repository contains the occurrence data used in the study:

**"Potential Distribution of *Bactericera cockerelli* in Colombia: Multi-Scale Epidemiological Approach Under Climate Change Scenarios"**

📂 Repository Structure
🗂️ Data

Bactericera cockerelli Colombia Oficial surveillance ICA.xlsx
Confirmed presence records of B. cockerelli (2021–2024) from ICA official surveillance in Nariño department.

Bactericera cockerelli Global.xlsx
Global occurrence data compiled from GBIF and published literature.

💻 Scripts
1. 🎯 Environmental Stratified Sampling

Generates background (pseudoabsence) points using PCA + K-means clustering to ensure balanced representation of environmental gradients and minimize spatial autocorrelation.

2. ⚙️ Model Calibration and Hyperparameter Tuning

MaxEnt_tuning.R: explores combinations of feature classes and regularization multipliers using a genetic algorithm for optimal model complexity.

RF_tuning.R: tunes number of trees and mtry parameter using out-of-bag error minimization with stratified cross-validation.

3. 🌡️ Weighted Ensemble of GCMs (2021–2040)

Implements a Weighted Ensemble Mean (WEM) combining ten CMIP6 models from WorldClim v2.1 using RMSE and correlation-based weights for Tmax, Tmin, and precipitation projections under SSP2-4.5 and SSP5-8.5 scenarios.

4. 🧭 MESS Analysis

Identifies climatically analogous zones across Colombia to constrain model projections within ecologically realistic environments.

5. 🧩 Spatially Explicit Pixel-Change Analysis

Compares binary habitat suitability maps across current and future climate scenarios to classify areas as:

🟦 Stable Areas

🟥 New Areas (↑ Risk)

🟩 Area Loss (↓ Risk)

6. 🧠 Model Evaluation and Ensemble Integration

Combines outputs from MaxEnt and RF to assess agreement and uncertainty, identifying high-confidence zones for B. cockerelli establishment.

🗺️ Data Sources

WorldClim v2.1 (https://worldclim.org
) – 19 bioclimatic predictors

ENVIREM (https://envirem.github.io/
) – additional environmental predictors

CMIP6 GCMs (2021–2040) – SSP2-4.5 and SSP5-8.5 scenarios

UPRA (https://sipra.upra.gov.co/nacional
) – potato crop suitability rasters for MESS reference zones

📊 Reproducibility

All analyses were performed in R using the terra, dismo, randomForest, ENMeval, envirem, and caret packages. The full workflow is designed for reproducibility, allowing users to replicate the entire modeling pipeline from environmental sampling to spatial projections.

## 📬 Contact
For questions or data requests, please contact:  
**Manuel Alejandro Cortes Quiceno** – manuel.cortes@ica.gov.co  
Instituto Colombiano Agropecuario (ICA), Subgerencia de Proteccion vegetal, Dirección Técnica de Epidemiología y Vigilancia Fitosanitaria, Bogotá, Colombia.
