# MCF-multtest
Implementation of the MCF-based large-scale multiple discrete testing method, and codes for the paper:

Xiaoyu Dai, Nan Lin, Daofeng Li, Ting Wang. A non-randomized procedure for heterogeneous discrete multiple testing based on randomized tests.

* BT.R: simulations of Binomial tests.
* Condition_check.R: plots to check Condition 3.1 in main paper, based on simulated datasets.
* FET.R: simulations of Fisher's Exact Tests.
* FNR_bound.R: plots of the FNR upper bound as in Theorem 2 in main paper, based on simulated datastes.
* MT.R: impementation of the MCF-based multiple testing method.
* Methy.R: main code for WGBS study.
* Methy_Z.R: plots of the standardized methylation level difference.
* Methy_addition_dif.R: plots of the differencen between p-value and p-value-minus.
* Methy_power.R: plots of the power of DMC detection.
* dmr_merge.R: merge DMCs into DMRs.
* mcf.cpp: cpp implemention for calculating MCF values.
* pi0.R: implemention for pi0 estimators.
* pi0_rep.R: comparison of differnt threshold in Storey's pi0 estimator.
* validation.R: validating the theoretical results in main paper.
