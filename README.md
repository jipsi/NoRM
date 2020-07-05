# NoRM
#Activation induced heterogeneity: Single cells in a macrophage population respond differently to stimulus. Upon secondary stimuli some of these cells can become non-responsive.
#This non-responsiveness can be permanent (at least over an observed period of time) or even reversible as seen by delayed kinetics. This code models to appearance and disappearance 
#of these proposed states (Dey et al, unpublished) by using a Doob-Gillespie algorithm. This code traces the time-evolution of these states and predicts transition rates that may
#be closest to predicting experimental data. These parameters are referred to as NoRM parameters.

#Time-evolution of these cell-states can be tracked using the following command (example):
x=fit_2_challenge_time_points_param_pies(param1, param2, param3, param4, param5)
#param1 and param2 are the name and formatted name of the protein under study
#param3 is can be any string (currently used to distinguish between different files) but is not essential
#param4 specifies if the simulations are to be run based on values obtained from a file or generated random numbers from 
#normal, uniform and negative binomial distributions to serve as test NoRM parameters. param4 value should be 1 to enable this.
#param5 should be numeric and specifies the minimum rmsd value that is acceptable to the model and empirical output.

