# Jacobs_et_al_CurrentBiology

This repository contains the analysis code for the Current Biology paper "Cortical state fluctuations during sensory decision making".

#### How the code is organised

The code is organised according to which figure the analyses are shown in, with the core code inside folders named after the figures, and any code that is required to generate the data to be used in the core code inside the "helpers" folder.

* Within the helpers folder is a subfolder called "experimentLists", which contains the structures with the details of the experiments that were used in the analyses. There is also a zipped subfolder called "ROIPixelSelections", which contains the pixel coordinates per dataset for the ROIs in the widefield imaging datasets that are called throughout the analyses.

Many scripts have been used in several figures: when this is the case, the script is in the folder for the figure where the analysis occurs for the first time. (This is why there are no Fig4 and Fig6 folders, as the code for the analysis shown in those figures is the same as in Figures 2 and 3; the only differences are the behavioural comparisons that are specified.) The list below explicitly explains which scrips this applies to.

* behaviourPerPower_imaging\
This script computes the probability of a certain response (Miss, Incorrect Choice, etc) as a function of low frequency power for the widefield imaging datasets. This type of plot occurs for the first time in Figure 2, and the script is therefore in that folder, but the same script was re-used for the analyses in Figures 3 and 4, and Supplementary Figures 2 and 6.
* examplePowerDifferenceMaps\
This script plots the example power difference maps in Figures 3-7, and Supplementary Figures 5 and 7. It is inside the folder for Figure 3.
* compute_powerDifferences\
This script computes the power difference between two behavioural conditions (for example Choice Miss, Correct Incorrect) per ROI. This analysis occurs for the first time in Figure 3, but is also applied in Figures 4-7, and Supplementary Figures 6-7.
* plot_powerDifferences\
This script plots the results from compute_powerDifferences.
* powerDifference_statistics\
Applies nested mixed effects models to the results from compute_powerDifferences.
* plot_powerDifference_comparisons\
This script plots the average Correct - Miss vs Incorrect - Miss, and Choice - Miss vs Correct - Incorrect power differences per ROI. This analysis has been repeated in Figure 7.

Please note that this code has been adapted from scripts that were structured slightly differently; I tried to make the ones published here more generalisable. If you encounter any bugs, please raise an issue and I’ll look into it.


#### Requirements
All the code is written in [Matlab](https://www.mathworks.com/products/matlab.html), using versions 2016b and 2018b. \
Some scripts call functions from other repositories on GitHub, mostly from the [Cortex-Lab](https://github.com/cortex-lab); but this is specified for each script.
