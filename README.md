# Code associated with "Resistance-minimizing strategies for introducing a novel antibiotic for gonorrhea treatment: a mathematical modeling study"

Authors: Reichert E, Yaesoubi R, Rönn MM, Gift TL, Salomon JA, Grad YH

Correspondence: ereichert@hsph.harvard.edu

Preprint: https://www.medrxiv.org/content/10.1101/2023.02.14.23285710v1

**Use Instructions:**

The three core files used to define the mathematical models and generate output are GCfunctions.R, calibration_ER.R, and GCtransmission.Rmd. 

- Custom functions that are used throughout are defined in GCfunctions.R. 
- The calibration code, calibration_ER, uses maximum likelihood estimation (MLE) to estimate a set of model parameters. Since the parameters that maximize the likelihood of a 3.0% mean gonorrhea prevalence at baseline are already defined in the GCtransmission notebook, this code only needs to be rerun if you would like to recalibrate the model to a different baseline equilibrium gonorrhea prevalence.
- The notebook "GCtransmission.Rmd" generates output for each of the antibiotic introduction strategies of interest under baseline model parameterization.

The folder "Baseline Visualization" contains code to further analyze and display the GCtransmission.Rmd output.

The folder "Sensitivity Analysis" contains code to regenerate model output under alternative parameterizations for Drug B and, to a lesser extent Drug A, within the SensAnalysis_Run.R script. The SensAnalysis_Viz.R script can then be used to analyze and visualize results of sensitivity analysis.
