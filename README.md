# Randomised-bcrm
The codes for the paper P. Mozgunov, T. Jaki and X. Paoletti (2018), Randomized dose-escalation designs for drug combination cancer trials with immunotherapy

This folder includes files for the paper 

P. Mozgunov, T.Jaki and X. Paoletti (2018) "Radnomized dose-escalation designs for drug combination cancer trials with immunotherapy", Journal of Biopharmaceutical Statistics.

The codes are the modification of the functions of the bcrm package by Michael Sweeting available on CRAN https://cran.r-project.org/web/packages/bcrm/. The functions of the original package are described in the paper

Sweeting, Michael, Adrian Mander, and Tony Sabin. "Bcrm: Bayesian continual reassessment method designs for phase I dose-finding trials." Journal of Statistical Software 54.13 (2013): 1-26.

which can be found here https://www.jstatsoft.org/article/view/v054i13/bcrm_Bayesian_Continual_Reassessment_Method_Designs_for_Phase_I_Dose-Finding_Trials.pdf

This folder includes the following files:

"emax_code" - the main bcrm function with the 4-parameter Emax model incorporated as one of the choice of the model with prior parameters specified as in the paper and the randomisation to the conrol arm. The file also includes the code for computing the benchmark estimates. This file should be run before proceeding with the next one.

"emax performance example" - the file demonstrates the bcrm functions with the Emax model. The randomisation ratio can be customised using parameters "cohort" (for the investigational arm) and "cohort.control" (for the control arm). The rest of the functions are as in the original package.

"BCRM package with the control group" includes the original bcrm function with the extenstion allowing randomisation to the control group.  The randomisation ratio can be customised using parameters "cohort" (for the investigational arm) and "cohort.control" (for the control arm). The rest of the functions are as in the original package. This file should be run before proceeding with th next one.

"L2R example" - the file demonstrates the performance of L2R model as specified in the original paper.




