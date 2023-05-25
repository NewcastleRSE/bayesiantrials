# bayesiantrials
This R package has been written to accompany the [Leveraging External Information in Clinical Trials](https://www.newcastle-biostatistics.com/courses/external_information/) short course provided by the Population Health Sciences Institute [Biostatistics Research Group (PHSI-BRG) at Newcastle University](https://www.newcastle-biostatistics.com). The functions in this package follow the approach detailed in Wilson (2022) for the Bayesian design and analysis of two-arm cluster randomised trials using assurance.

# Installation

If not already installed you will first need to manually install the Bayesian sampling program JAGS. Please see the [official JAGS documentation](https://mcmc-jags.sourceforge.io) for instructions on how to do this.

Next you can install the latest version of the package from Github using the following R commands:
```R
install.packages("devtools")
devtools::install_github("NewcastleRSE/bayesiantrials")
```

# References
Wilson, Kevin J. 2022. "Bayesian Design and Analysis of Two-Arm Cluster Randomised Trials Using Assurance." arXiv Preprint arXiv:2208.12509. [https://arxiv.org/abs/2208.12509](https://arxiv.org/abs/2208.12509).
