## 🐛 Potential Distribution of Bactericera cockerelli in Colombia

This repository accompanies the manuscript:
“Potential Distribution of Bactericera cockerelli (Sulc) in Colombia: Multi-Scale Epidemiological Approach Under Climate Change Scenarios”

All scripts and datasets support transparent, reproducible, and scalable modeling for phytosanitary surveillance and climate-informed pest risk assessment.

## 📂 Repository Structure
## 🗂️ Data

Bactericera cockerelli Colombia Oficial surveillance ICA.xlsx
Confirmed presence records of B. cockerelli (2021–2024) from ICA official surveillance in Nariño department.

Bactericera cockerelli Global.xlsx
Global occurrence data compiled from GBIF and published literature.

## 💻 Scripts
## 🎯 1. Stratified Sampling in space Environmental

Generates background (pseudoabsence) points using PCA + K-means clustering to ensure balanced representation of environmental gradients and minimize spatial autocorrelation.

## ⚙️ 2. Model Calibration and Hyperparameter Tuning

MaxEnt_tuning.R – explores combinations of feature classes and regularization multipliers using a genetic algorithm to optimize model complexity.

RF_tuning.R – tunes the number of trees and mtry parameter using out-of-bag error minimization with stratified cross-validation.

## 🌡️ 3. Ensemble of GCMs (2021–2040)

Implements a Weighted Ensemble Mean (WEM) combining ten CMIP6 models from WorldClim v2.1.
Weights are assigned based on Root Mean Square Error (RMSE) and correlation coefficients, producing ensemble projections for Tmax, Tmin, and precipitation under SSP2-4.5 and SSP5-8.5.

## 🧭 4. Multivariate Environment Similary Superface Analysis

Identifies climatically analogous zones across Colombia to constrain model projections within ecologically realistic conditions, reducing extrapolation uncertainty.

## 🧩 5. Spatially Explicit Pixel-Change Analysis

Compares binary habitat suitability maps between scenarios (Current vs SSP2-4.5, SSP5-8.5) to classify:

🟦 Stable Areas

🟥 New Areas (↑ Risk)

🟩 Area Loss (↓ Risk)

This approach quantifies spatial risk transitions and highlights areas of potential vector expansion.

## 🧠 6. Model Evaluation and Ensemble Integration

Integrates MaxEnt and Random Forest model outputs to assess inter-algorithm agreement, identifying high-confidence establishment zones and quantifying uncertainty in presence-only modeling.

🗺️ Data Sources

🌍 WorldClim v2.1 – 19 bioclimatic predictors (https://worldclim.org))

🌿 ENVIREM – additional environmental predictors (https://envirem.github.io/))

🌡️ CMIP6 GCMs (2021–2040) – SSP2-4.5 and SSP5-8.5 projections

🥔 UPRA (Colombia) – potato crop suitability rasters for MESS reference zones (https://sipra.upra.gov.co/nacional))

## 📊 Reproducibility

All analyses were performed in R using the following packages:
terra, dismo, tidyverse, tidyterra, SDMtune, randomForest, envirem, and RamdonForestExplainer.
The complete workflow is reproducible from environmental sampling to final spatial projections.

## 📬 Contact
For questions or data requests, please contact:  
**Manuel Alejandro Cortes Quiceno** – manuel.cortes@ica.gov.co  
Instituto Colombiano Agropecuario (ICA), Subgerencia de Proteccion vegetal, Dirección Técnica de Epidemiología y Vigilancia Fitosanitaria, Bogotá, Colombia.
