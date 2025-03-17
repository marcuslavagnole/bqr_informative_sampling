This repository provides the R routines from my doctoral thesis and its main spin-offs.

### 1. Bayesian quantile regression models for complex survey data under informative sampling[^a]

The repo includes:

- **MCMC_BWQR_AL.R**: MCMC routine for the weighted quantile regression based on the Asymmetric Laplace distribution;
- **MCMC_BWQR_SL.R**: MCMC routine for the weighted quantile regression based on the score likelihood;
- **MCMC_BWQR_AP.R**: MCMC routine for the weighted quantile regression based on the approximate method;
- **data_provabrasil.txt**: dataset utilized in the real-data-based simulation study.

<!--
2. **Bayesian quantile regression models for bounded count data under informative sampling**. 

- **MCMC_BWQR_AL_count.R** : MCMC routine for the weighted quantile regression based on the Asymmetric Laplace distribution for count data;
- **MCMC_BWQR_AL_bounded_count.R** : MCMC routine for the weighted quantile regression based on the Asymmetric Laplace distribution for bounded count data.
- **MCMC_BWQR_PL_count.R** : MCMC routine for the weighted quantile regression based on the pseudo posterior for count data;
- **MCMC_BWQR_PL_bounded_count.R** : MCMC routine for the weighted quantile regression based on the pseudo posterior for bounded count data.
-->

### 2. A Bayesian approach to multiple-output quantile regression analysis under informative sampling

The repo includes:

- **EM_BWQR_AL_MO.R**: EM routine for the multiple-output weighted quantile regression based on the Asymmetric Laplace distribution.
- **data_nhds.txt**: dataset utilized as a motivating example.

### 3. An Expectation-Maximization algorithm for noncrossing Bayesian quantile regression analysis under informative sampling

The repo includes:

- **NonCrossingBWQR_AL.R**: routine for noncrossing Bayesian quantile regression analysis under informative sampling. File **MCMC_BWQR_AL.R** is required to run the function.
- **data_nhds_rural_northeast.txt**: dataset utilized as a real-data illustration.

[^a]: Nascimento, M. L., Gonçalves, K. C. M. [Bayesian quantile regression models for complex survey data under informative sampling](https://doi.org/10.1093/jssam/smae015). Journal of Survey Statistics and Methodology, 12(4), 1105-1130. 
