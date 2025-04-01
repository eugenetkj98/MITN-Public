# Welcome Page

This is the documentation for the ITN Stock and Flow (SNF) model for the Malaria Atlas Project.

Author: Eugene Tan

Last Updated: 6/8/2024 

## Overview

The MAP ITN model attempts to estimate the quantity, coverage and effects of insecticide treated nets within African countries. To do this, estimates regarding the actual quantities of net crop in the country are required. This model attempts to infer this from various sources of net delivery and distribution data, IMH population surveys, DHS household surveys among several others. The first iteration of this model was first proposed by Bhatt et al. (2015) and used to estimate net coverage dynamics over the years 2000-2015. This was extended upon by Bertozzi-Villa et al. (2021) to model net access. They also proposed additional measures of access gap and use gap, which are subsequently used to infer local variations in net coverage and usage, attrition times and the construction of a geospatial map of net coverage. 

This current iteration of the MAP ITN aims to clean up and refine on the Bertozzi-Villa (BV) model to be more interpretable for third-party use. Additionally, it also aims to maximise the amount of information that can be extracted from existing datasets, while extending for newer and richer incoming datasets that are of a higher quality and resolution. See below for a non-exhaustive list of aims.

**Implementation Aims and Roadmap**
- Interpretability
    - Creating full documentation for the model to allow better accessibility by third-party and successive inheritors
    - Detailed write-ups on mathematics and theory of the models for future potential improvements
    - Defining required structures of input datasets for future users
- Extensibility
    - Implementing the model that utilises optimum languages for each component task (e.g. Julia - algorithm design, data wrangling and machine learning; R - regression, complex Bayesian models)
    - Outputting intermediate and final output datasets with standardised structures and components that may be reused for other tasks and workflows.
- Mathematical Extensions
    - Maximising information extraction from existing datasets (movement to monthly resolution)
    - Integrating transition to analysis on new incoming datasets: higher resolution, no need for delivery data
    - Adding support for multiple net types with different efficacies
    - Extending analysis from national (ADMIN0) to subnational (ADMIN1)
