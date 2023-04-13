# Likelihood-Free-Population-Calibration

Matlab and Julia code to support the article Population Calibration using Likelihood-Free Bayesian Inference (https://arxiv.org/pdf/2202.01962.pdf).

Each folder contains code for an example of the paper:

- mixture example: start with run_sl.m
- growth example: start with run.m
- internalisation example: This example is taken from a previous paper (Browning et al 2021) where we originally applied our method. 
                           Calls code available at https://github.com/ap-browning/internalisation (run from the root directory).
                           The file InternalisationPredictive.jl then produces more results.

Within each folder there is a produce_results.m file that can be run to reproduce the figures in the paper.

Some files to facilitate the visualisation of the results have been obtained from freely available sources:

Aslak Grinsted (2022). Subaxis - Subplot (https://www.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot), 
MATLAB Central File Exchange. 

Nico Schlömer (2022). matlab2tikz/matlab2tikz (https://github.com/matlab2tikz/matlab2tikz), GitHub.

Other References

Browning, A. P., Ansari, N., Drovandi, C., Johnston, A., Simpson, M. J., and Jenner, A. L. (2021).
Identifying cell-to-cell variability in internalisation using flow cytometry. bioRxiv.

