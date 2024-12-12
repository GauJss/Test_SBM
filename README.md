# Test_SBM
This file provides the codes of our proposed testing methods in "Testing the Stochastic Block Models Based on Maximum Sampling Entry-Wise Deviations".

Test_I contains the codes for testing hypothesis (I)：
"Simu_Main.m": The central script orchestrating the simulation process, calling other functions as needed.
"Generating.m": Responsible for generating network data based on predefined parameters and models.
"MLE.m": Computes the maximum likelihood estimators for the edge probability matrix $Q$ and the entry-wise deviations $\rho_{iv}$.
"SP_kmeans": Implements spectral clustering to estimate the membership vector $g$.
"Zhang_Aug.m": Calculates the augmented test statistic proposed by Hu et al. (2021).
"Zhang_AugBoot.m": Computes the bootstrap-corrected augmented test statistic introduced by Hu et al. (2021).
"Proposed.m": Implements the proposed test statistic $T_n$.
"Proposed_Boot.m": Computes the proposed bootstrap-corrected test statistic $T_{n,boot}$
"Proposed_AugBoot.m": Computes the proposed bootstrap-corrected augmented test statistic $T_{n,boot}^+$.


Test_II contains the codes for testing hypothesis (II):
"Simu_Main.m": Serves as the primary script coordinating the simulation process and invoking other supporting functions as required.
"Generating.m": Generates network data using predefined parameters and models.
"Corrupted_Zhang.m": Calculates the test statistic proposed by Hu et al. (2021), considering z% corruption in the membership vector 
$g$ (z=0,5,10).
"Corrupted_Zhang_AugBoot.m": Computes the bootstrap-corrected augmented test statistic by Hu et al. (2021), incorporating z% corruption in $g$.
"Corrupted_Proposed.m": Implements the proposed test statistic $T_n$ while accounting for z% corruption in $g$ corrupted.
"Corrupted_Proposed_AugBoot.m": Computes the proposed bootstrap-corrected augmented test statistic $T_{n,boot}^+$, considering z% corruption in $g$.
