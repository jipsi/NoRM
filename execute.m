%config
%%%%section bmdm - sd
minimum_deviation=0.2;
tpts_bmdm=[4 8 12];
tpts_raw=[8 12 16];
search=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BMDM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%CHANGE TIMEPOINTS TO 4,8,12 before running%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BMDM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%3 state
fit_2_challenge_time_points_param_pies_test('tnf', 'TNF', 'bmdm', 3, search, 0.05, tpts_bmdm)
% fit_2_challenge_time_points_param_pies_test('il1b', 'pro-IL1b', 'bmdm', 3, search, 0.2, tpts_bmdm)
% fit_2_challenge_time_points_param_pies_test('il6', 'IL-6', 'bmdm', 3, search, 0.08, tpts_bmdm)
% fit_2_challenge_time_points_param_pies_test('nos2', 'NOS2', 'bmdm', 3, search, 0.2, tpts_bmdm)
% %%%4 state
fit_2_challenge_time_points_param_pies_test('tnf', 'TNF', 'bmdm', 4, search, 0.05, tpts_bmdm)
% fit_2_challenge_time_points_param_pies_test('il1b', 'pro-IL1b', 'bmdm', 4, search, 0.2, tpts_bmdm)
% fit_2_challenge_time_points_param_pies_test('il6', 'IL-6', 'bmdm', 4, search, 0.08, tpts_bmdm)
% fit_2_challenge_time_points_param_pies_test('nos2', 'NOS2', 'bmdm', 4, search, 0.2, tpts_bmdm)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RAW%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%CHANGE TIMEPOINTS TO 8,12,16 before running%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RAW%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%3 state
%fit_2_challenge_time_points_param_pies_test('tnf', 'TNF', 'raw_rep_avg', 3, search, 0.15, tpts_raw)
%fit_2_challenge_time_points_param_pies_test('il1b', 'pro-IL1b', 'raw_rep_avg', 3, search, 0.15, tpts_raw)
%fit_2_challenge_time_points_param_pies_test('il6', 'IL-6', 'raw_rep_avg', 3, search, 0.2, tpts_raw)
%fit_2_challenge_time_points_param_pies_test('nos2', 'NOS2', 'raw_rep_avg', 3, search, 0.17, tpts_raw)
% %%%%%4 state
% %%0.05 not goodnano
%fit_2_challenge_time_points_param_pies_test('tnf', 'TNF', 'raw_rep_avg', 4, search, 0.2, tpts_raw)
%fit_2_challenge_time_points_param_pies_test('il1b', 'pro-IL1b', 'raw_rep_avg', 4, search, 0.15, tpts_raw)
%fit_2_challenge_time_points_param_pies_test('il6', 'IL-6', 'raw_rep_avg', 4, search, 0.2, tpts_raw)
%fit_2_challenge_time_points_param_pies_test('nos2', 'NOS2', 'raw_rep_avg', 4, search, 0.17, tpts_raw)
