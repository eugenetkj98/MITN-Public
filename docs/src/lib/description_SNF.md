# Net Crop Stock and Flow (SNF)

The net crop stock and flow (SNF) component of the model infers estimated values of national net crop based on reported net delivery data (organisations to governments), distribution data (governments to communities) and household survey data (actual reported net access and usage). This model is also used to learn the failure behaviour of insecticide treated nets (ITNs) and their variances across different countries.

## Compartmental Stock and Flow (SNF)


## Upsampling and Inferring Monthly Distribution


## Model Fitting


## Net Access Interpolation


Structure
- Aim of section -> To estimate net crop at any given year
- Compartmental Stock and Flow (SNF)
- Net Crop Model and Matching Distribution
    - Upsampling to monthly resolution -> Inferring distribution times
- Model Fitting and Bayesian Inference
    - RWMH MCMC for parameters in stock and flow
    - Gradient Descent with custom loss function for monthly distribution weights
    - Need to justify reason for DIC
    - Iterative scheme until convergence
- ITN Net Access Interpolation
    - How to generate forward -> Go through with code