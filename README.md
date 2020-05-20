# MOUSE_LFPMEG
Data and scripts from Patai, Z. E., Foltynie, T., Limousin, P., Hariz, M. I., Zrinzo, L., Bogacz, R., Litvak, V. (2020). Conflict detection in a sequential decision task is associated with increased cortico-subthalamic coherence and prolonged subthalamic oscillatory response in the beta band. (In Submission).

Scripts for modelling written by RB and EZP, LFP/MEG analysis by RB, VL and EZP.

Requirements: Matlab, SPM12 & Fieldtrip (for LFP/MEG data - note: SPM12 includes FT)

Scripts for behavioural and LFP analysis, and LFP data. MEG data is available on: MRC BDNU 

STN_behaviour.m:  behavioural analysis

STN_behaviour_model_comparisons.m: results of modelling behaviour with explicit mathematical models (uses scripts M1-M6 parametersearch) - supplemental data

compare_strategies.m: separating behaviour into those who integrate evidence or just respond to 'same' stimuli in a row

plot_same_same.m: recreates Figure 2

MEGLFP_coherence.m: recreates Figure 4


extra support scripts:
barwitherr,
trial_predictors_stim_by_choice,
plotpatch
