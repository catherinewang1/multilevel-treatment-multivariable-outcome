# multilevel-treatment-multivariable-outcome
analysis for scCRISPR screens with multilevel treatments and multivariable outcomes. focusing on shrinkage



# `notebooks/shrinkEstimates.qmd` Notebook


Notebook to analyze the cleary dataset with a combination of methods (SCEPTRE, Poisson, Negative Binomial) with clustering (kmeans and mixture modeling) and robust empirical bayes shrinkage.

The structure of the analysis and notebook is:

-   [Perform SCEPTRE w/ cell sample splitting](#sceptre)
-   [Perform GLM (Poisson and Negative Binomial) w/ cell sample splitting](#glm)
-   [Perform robust shrinkage (on estimates or tstat, not decided yet) after clustering (kmeans or mixture models)](#ebci)
-   [Plot Results](#plots)

This is a clean version of previous `robustempiricalbayestest.qmd`, which has more unnecessary code which is now stored in utils files. See the utils and/or previous version for more details.

## some details

-   uses the util code files:

    -   `../utils/perform_sceptre_cleary.r` (get the original estimates)
    -   `../utils/cluster_and_ebci_shrinkage.r` (perform shrinkage)

-   saves plots in `../plots/`

-   saves objects (rds objects, csv, sceptre object, etc...) in `../saves/`

-   Each section should save the important objects, and then the following sections can just load the relevant objects. This is done so that each section can be run at different times because each section might take a while to run, and then it is easy to restart at any section.

## SCEPTRE {#sceptre}

SCEPTRE (w/ cell sample splitting)

Main created objects to be saved are

-   `../saves/sceptre/sceptre_obj_train.rds`
-   `../saves/sceptre/sceptre_obj_test.rds`
-   `../saves/sceptre/sceptre_obj_all.rds`



## GLM {#glm}

GLM (w/ sample splitting) on same tests as SCEPTRE

Main created objects to be saved are

-   `../saves/poisson_effects_train.csv`
-   `../saves/poisson_effects_test.csv`
-   `../saves/poisson_effects_all.csv`
-   `../saves/negativebinomial_effects_train.csv`
-   `../saves/negativebinomial_effects_test.csv`
-   `../saves/negativebinomial_effects_all.csv`


## EBCI Shrinkage {#ebci}

EBCI Shrinkage on SCEPTRE + GLM (w/ kmeans and mixture model clustering)



## Plot Results {#plots}




