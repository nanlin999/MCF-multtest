# MCF-multtest
Implementation of the MCF-based large-scale multiple discrete testing method, and codes for the paper:

Xiaoyu Dai, Nan Lin, Daofeng Li, Ting Wang. A non-randomized procedure for heterogeneous discrete multiple testing based on randomized tests.

* BT.R: Simulations of Binomial tests.
* Condition_check.R: plots to check Condition 3.1 in main paper, based on simulated datasets.
* FET.R: Simulations of Fisher's Exact Tests.
* FNR_bound.R: plots of the FNR upper bound as in Theorem 2 in main paper, based on simulated datastes.
* MT.R: Impementation of the MCF-based multiple testing method.
* Methy.R: Main code for WGBS study.
* Methy_Z.R: Plots of the standardized methylation level difference.
* Methy_addition_dif.R: Plots of the differencen between p-value and p-value-minus.
* Methy_power.R: Plots of the power of DMC detection.
* dmr_merge.R: Merge DMCs into DMRs.
* mcf.cpp: Cpp implemention for calculating MCF values.
* pi0.R: Implemention for pi0 estimators.
* pi0_rep.R: Comparison of differnt threshold in Storey's pi0 estimator.
* validation.R: Validating our theoretical results.
