This repository provides the R routines for the articles derived from my doctoral thesis and its main spin-offs. All the following projects are joint work with my supervisor [Kelly C. M. Gonçalves](https://sites.google.com/dme.ufrj.br/kelly/). 

### 1. Bayesian quantile regression models for complex survey data under informative sampling[^a]

The interest in considering the relation among random variables in quantiles instead of the mean has emerged in various fields, and data collected from complex survey designs are of fundamental importance to different areas. Despite the extensive literature on survey data analysis and quantile regression models, research papers exploring quantile regression estimation accounting for an informative design have primarily been restricted to a frequentist framework. The paper introduces different Bayesian methods relying on the survey-weighted estimator and the estimating equations. These tools are illustrated in a design-based simulation study that uses data from Prova Brasil 2011. 

The repo includes:

- **MCMC_BWQR_AL.R**: MCMC routine for the weighted quantile regression based on the Asymmetric Laplace distribution;
- **MCMC_BWQR_SL.R**: MCMC routine for the weighted quantile regression based on the score likelihood;
- **MCMC_BWQR_AP.R**: MCMC routine for the weighted quantile regression based on the approximate method;
- **data_provabrasil.txt**: dataset utilized in the real-data-based simulation study.

<!--
2. **Bayesian quantile regression models for bounded count data under informative sampling**. 

- **MCMC_BWQR_AL_count.R**: MCMC routine for the weighted quantile regression based on the Asymmetric Laplace distribution for count data;
- **MCMC_BWQR_AL_bounded_count.R**: MCMC routine for the weighted quantile regression based on the Asymmetric Laplace distribution for bounded count data.
- **MCMC_BWQR_PL_count.R**: MCMC routine for the weighted quantile regression based on the pseudo posterior for count data;
- **MCMC_BWQR_PL_bounded_count.R**: MCMC routine for the weighted quantile regression based on the pseudo posterior for bounded count data.
-->

### 2. A Bayesian approach to multiple-output quantile regression analysis under informative sampling

The paper presents a Bayesian multiple-output quantile regression for complex survey data under informative sampling. Our approach relies on the asymmetric Laplace distributional assumption. From the location-scale mixture representation of this distribution, we introduce an Expectation–Maximization algorithm that provides a less computationally intensive alternative to the commonly used Markov Chain Monte Carlo algorithm for posterior inference. Our developments are mainly motivated by the joint analysis of growth indexes from Brazilian children under five.

The repo includes:

- **EM_BWQR_AL_MO.R**: EM routine for the multiple-output weighted quantile regression based on the Asymmetric Laplace distribution.
- **data_nhds.txt**: dataset utilized as a motivating example.

### 3. An Expectation-Maximization algorithm for noncrossing Bayesian quantile regression analysis under informative sampling

When quantiles are fitted separately, the resultant regression lines may cross, violating the basic probabilistic rule that quantiles are monotonic functions and possibly causing problems for inference and interpretation in practice. The paper introduces a method for handling crossing issues regarding the analysis of complex survey data under informative sampling. Using the location-scale mixture representation of the asymmetric Laplace distribution, we write a joint posterior density function for the quantile levels of interest and develop a constrained Expectation-Maximization algorithm. Data from the Brazilian National Demographic Health Survey of Women and Children is analyzed to verify and illustrate the algorithm's effectiveness.

The repo includes:

- **NonCrossingBWQR_AL.R**: routine for noncrossing Bayesian quantile regression analysis under informative sampling. File **MCMC_BWQR_AL.R** is required to run the function.
- **data_nhds_rural_northeast.txt**: dataset utilized as a real-data illustration.

[^a]: Nascimento, M. L., Gonçalves, K. C. M. [Bayesian quantile regression models for complex survey data under informative sampling](https://doi.org/10.1093/jssam/smae015). Journal of Survey Statistics and Methodology, 12(4), 1105-1130. 
