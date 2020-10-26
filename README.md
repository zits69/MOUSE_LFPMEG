# MOUSE_LFPMEG
Data and scripts from Patai, Z. E., Foltynie, T., Limousin, P., Hariz, M. I., Zrinzo, L., Bogacz, R., Litvak, V. (2020). Conflict detection in a sequential decision task is associated with increased cortico-subthalamic coherence and prolonged subthalamic oscillatory response in the beta band. (In Submission).

Scripts for modelling written by RB and EZP, LFP/MEG analysis by RB, VL and EZP.

Requirements: Matlab, SPM12 & Fieldtrip (for LFP/MEG data - note: SPM12 includes FT)

Scripts for behavioural and LFP analysis. 
LFP (patients) data is available on: https://data.mrc.ox.ac.uk/data-set/human-lfp-recordings-stn-during-sequential-conflict-task
MEG (controls) data is available on: https://openneuro.org/datasets/ds002908/versions/1.0.0

STN_behaviour.m:  behavioural analysis

STN_behaviour_model_comparisons.m: results of modelling behaviour with explicit mathematical models (uses scripts M1-M4 parametersearch) - supplemental data

compare_strategies.m: separating behaviour into those who integrate evidence or just respond to 'same' stimuli in a row

plot_same_same.m: recreates Figure 2

MEGLFP_coherence.m: recreates Figure 4


extra support scripts:
barwitherr,
trial_predictors_stim_by_choice,
plotpatch
