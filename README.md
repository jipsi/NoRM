# NoRM
Activation induced heterogeneity: Single cells in a macrophage population respond differently to stimulus. Upon secondary stimuli some of these cells can become non-responsive.
This non-responsiveness can be permanent (at least over an observed period of time) or even reversible as seen by delayed kinetics. This code models to appearance and disappearance of these proposed states (Dey et al, unpublished) by using a Doob-Gillespie algorithm. This code traces the time-evolution of these states and predicts transition rates that maybe closest to predicting experimental data. These parameters are referred to as NoRM parameters.

Priors
*You will require all csv files to be in the same folder as your working directory
*Select a relative/absolute path in fit_2_challenge_time_points_param_pies_test.m
*Select a relative/absolute path for saving results in 

Start at execute.m
*Executes fit_2_challenge_time_points_param_pies_test.m with intialising parameters

minimum_deviation=0.2; #minimum deviation from experimental results
tpts_bmdm=[4 8 12]; #timepoints for BMDMs
tpts_raw=[8 12 16]; #timepoints for RAW cells
search=0; #decides whether to search for parameters or select parameters saved in a csv

