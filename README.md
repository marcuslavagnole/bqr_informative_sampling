Directory with R routines used in three main chapters of my doctoral dissertation:

1. **Bayesian quantile regression models for complex survey data under informative sampling**, Journal of Survey Statistics and Methodology, 12(4), 1105–1130, 2024. 

- **MCMC_BWQR_AL.R** : MCMC routine for the weighted quantile regression based on the Asymmetric Laplace distribution;
- **MCMC_BWQR_SL.R** : MCMC routine for the weighted quantile regression based on the score likelihood;
- **MCMC_BWQR_AP.R** : MCMC routine for the weighted quantile regression based on the approximate method.

2. **Bayesian quantile regression models for bounded count data under informative sampling**. 

- **MCMC_BWQR_AL_count.R** : MCMC routine for the weighted quantile regression based on the Asymmetric Laplace distribution for count data;
- **MCMC_BWQR_AL_bounded_count.R** : MCMC routine for the weighted quantile regression based on the Asymmetric Laplace distribution for bounded count data.
- **MCMC_BWQR_PL_count.R** : MCMC routine for the weighted quantile regression based on the pseudo posterior for count data;
- **MCMC_BWQR_PL_bounded_count.R** : MCMC routine for the weighted quantile regression based on the pseudo posterior for bounded count data.

3. **Bayesian multiple-output quantile regression for complex survey data under informative sampling**.

- **EM_BWQR_AL_Mult** : EM algorithm for the Bayesian multiple-output quantile regression.
