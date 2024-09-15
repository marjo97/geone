# Geostatistical Analysis and Kriging Interpolation Script
## Overview
This Octave script is designed to perform various geostatistical analyses and kriging interpolation, making use of the open-source capabilities of Octave (version `8.3.0`). It includes functionality for data visualization, spatial analysis, variogram modeling, parameter estimation via the maximum likelihood method, cross-validation, and error analysis.

The script provides flexibility through adjustable parameters, allowing users to tailor the code to their specific geostatistical problems.

### **Authors**

* Maria-Despoina Koltsidopoulou
* Andreas Pavlidis
* Emmanouil Varouchakis
* Dionisis Hristopoulos

### **Requirements**

- Octave Version: >= 4.0
- Toolboxes: Statistics and Optimization Toolboxes are required for full functionality.

### **Key Features**

1. **Data Loading and Preprocessing**: Loads data from .mat files and organizes it for geostatistical analysis.
2. **Data Visualization**: Includes scatter plots, histograms, and normal probability plots to visualize spatial distributions and statistical properties of the data.
3. **Statistical Analysis**: Computes basic statistics like mean, median, standard deviation, skewness, kurtosis, and empirical cumulative distribution function (CDF) comparisons.
4. **Variogram Modeling**: Performs variogram fitting using various models (e.g., exponential, SSRF) and includes a nugget effect.
5. **Maximum Likelihood Parameter Estimation**: Implements parameter estimation for variogram models using the Maximum Likelihood method.
6. **Cross-validation and Error Analysis**: Performs cross-validation, computes errors such as Mean Absolute Error (MAE) and Root Mean Squared Error (RMSE), and identifies outliers.
7. **Kriging Interpolation**: Uses ordinary kriging to predict values at unsampled locations and visualizes the results, including uncertainty measures.

### **Usage Instructions**

1. **Setting up Models and Parameters**: The script includes adjustable parameters for different models. You can change the `MODEL` variable and other parameters to customize the analysis for your specific needs.

- Default Model: Exponential (`MODEL = 1`)
- Distance Model: Minkowski distance (`dist_model = 'minkowski'`)
- Step Size: Can be adjusted in the meshgrid setup (`step_size = 0.1`)

2. **Loading Data:** Make sure your `.mat` files are in the correct path specified in the script (`base_path = 'NormalizedData2'`).

3. **Running the Script**: Run the script directly in Octave. The workspace will be cleared at the start (`clear all; close all; clc;`), and various sections will execute in sequence.

4. **Customization**: You can modify several indicative values to suit the dataset or problem being addressed, such as:

- Number of neighbors in kriging.
- Maximum distance for neighbors (`maxdis`).
- Model parameters (variance, range, nugget).

### Output
- **Visualizations**: Scatter plots, histograms, variograms, and kriging results are visualized in figures.
- **Error Metrics**: MAE, RMSE, correlation coefficients, and Mean Error (ME) are printed to the console for model evaluation.

### Contact
For questions or issues, please reach out to the authors listed above.
