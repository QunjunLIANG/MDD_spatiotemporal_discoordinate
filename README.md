# Spatiotemporal topology abnomality in MDD

#### Description
This repository provides the code and original data used for analyzing the spatiotemporal topology of the MDD brain.

#### File Architecture

**/Analyze_script** contains the R code to run the main analysis, including spatiotemporal topology estimation,
machine learning classification, inter-subject representational similarity analysis, etc.

**/Analyze_script/Analyz+number...** are the codes for cohort1 HC-MDD dataset (91 matched pairs of subject). Number indicates the analysis sequence.

**/Analyze_script/AnalyzValidation+number...** are the codes for testing the relationship between spatiotemporal topology and disease severity with stand-alone MDD data in cohort1 (114 MDD subject).

**/Analyze_script/cohort2_Analyz+number...**  are the codes for reproduciable analysis with cohort2 (the REST-meta-MDD Project from DIRECT Consortium on [ScienceDB](https://www.scidb.cn/en/detail?dataSetId=cbeb3c7124bf47a6af7b3236a3aaf3a8))

**/Analyze_script/function_...**  are the custome R functions for this project.

**/inputs** contains the time delay matrix or raw functional gradient data utilized for estimating spatiotemporal topology.

**/Python_MATLAB_scripts** includes the Python script utilized for functional gradient estimation, as well as the MATLAB script employed for time delay estimation.

#### contact information

Please contact me with the email address [here](2015021358@m.scnu.edu.cn).

It should be noted that my expertise lies in neuroscience rather than software engineering, meaning that my capacity for maintaining this repository is limited.

For programming-related issues, please try searching on Google or visiting ChatGPT for assistance.


