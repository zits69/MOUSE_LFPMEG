# MOUSE_LFPMEG
Data and scripts from Patai, E. Z., Foltynie, T., Limousin, P., Akram, H., Zrinzo, L., Bogacz, R., & Litvak, V. (2022). Conflict Detection in a Sequential Decision Task Is Associated with Increased Cortico-Subthalamic Coherence and Prolonged Subthalamic Oscillatory Response in the β Band. Journal of Neuroscience, 42(23), 4681–4692. https://doi.org/10.1523/JNEUROSCI.0572-21.2022

https://www.jneurosci.org/content/42/23/4681

Scripts for modelling written by RB and EZP, LFP/MEG analysis by RB, VL and EZP.

Requirements: Matlab, SPM12 & Fieldtrip (for LFP/MEG data - note: SPM12 includes FT)

Scripts for behavioural and LFP analysis. 
LFP (patients) data is available on: https://data.mrc.ox.ac.uk/data-set/human-lfp-recordings-stn-during-sequential-conflict-task
MEG (controls) data is available on: https://openneuro.org/datasets/ds002908/versions/1.0.0

STN_behaviour.m:  behavioural analysis

STN_behaviour_model_comparisons.m: results of modelling behaviour with explicit mathematical models (uses scripts M1-M4 parametersearch) - supplemental data

compare_strategies.m: separating behaviour into those who integrate evidence or just respond to 'same' stimuli in a row


lfp_regression_allcontacts_figures: recreated Figure 2

plot_same_same.m: recreates Figure 3

MEGLFP_coherence_final.m: recreates Figure 5


extra support scripts:
barwitherr,
trial_predictors_stim_by_choice,
plotpatch
