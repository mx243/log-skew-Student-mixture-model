# log-skew-Student-mixture-model
This repository contains a Matlab implementation of a Bayesian network model. The model is capable of fiting given data (1D or higher-dimensional) using a linear combination of multiple skew-Student distributions. It has been used to effectively model the joint distribution of prostate cancer patients' log time of survival and their various biomarker measurements, and therefore predict a patient's survival based on his biomarker levels.

A full description of the model can be found in section B of [[1]](#1).

The current version of the model is tuned to model the distributions of log cytokine levels of tuberculous meningitis patients. Instructions on how to use the software for cytokine modelling can be found in Instructions.pdf.

The software is joint work by Dr Roger Sewell (rfs34), Mingtong Xu (mx243), Jacob Coxon (jc2062) and Tommy Walker Mackay (tow24).

## References
<a id="1">[1]</a> 
Roger Sewell. 
Assessment of the quality of a prediction. 2025. 
arXiv: 2404.15764 [math.ST].
URL: https://arxiv.org/abs/2404.15764.
