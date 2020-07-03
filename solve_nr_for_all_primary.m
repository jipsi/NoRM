function fun_solve_nr_for_all_primary(experiment, cytokine)
%%drive path
str_path=string_path;
%'C:\Users\shoum\Google Drive\silence\';
%'D:\GoogleDrive\silence\';
%%The common parameters
cell_in_soup=[0 0 100 0];
plusminus=15;
simulations=3e6;
%the not so common parameters
startTime=0;
stopTime=24;
%pick experiment
exp_str='r75';
%pick cytokine
cyto_str='tnf';
%make file name
exp_cyto=append(exp_str,'_',cyto_str);
%%load_het - experimental data loads here
expTab=readtable(append(str_path,'PhD\MATLAB\8_state_models\May_2019\Individual_model_2_non_responsive_states\challenge_fitter\',exp_cyto,'.csv'),'ReadVariableNames',true);
%NoT 
expData=[expTab.LPS1_NoT expTab.LPS10_NoT expTab.LPS100_NoT expTab.LPS1000_NoT]/100;
x=fun_sim_explore_nr_primary(cell_in_soup, startTime, stopTime, expData, plusminus,simulations);
csvwrite(append(str_path,'PhD\MATLAB\8_state_models\May_2019\Individual_model_2_non_responsive_states\challenge_fitter\results\',cyto_str,'\',exp_cyto,'_estimates_NoT.csv'),x)
%csvwrite(append(str_path,'PhD\MATLAB\8_state_models\May_2019\Individual_model_2_non_responsive_states\challenge_fitter\results\',cyto_str,'\',exp_cyto,'NoTTESTING.csv'),x)

%NTC
expData=[expTab.LPS1_NTC expTab.LPS10_NTC expTab.LPS100_NTC expTab.LPS1000_NTC]/100;
x=fun_sim_explore_nr_primary(cell_in_soup, startTime, stopTime, expData, plusminus,simulations);
csvwrite(append(str_path,'PhD\MATLAB\8_state_models\May_2019\Individual_model_2_non_responsive_states\challenge_fitter\results\',cyto_str,'\',exp_cyto,'_estimates_NTC.csv'),x)
%csvwrite(append(str_path,'PhD\MATLAB\8_state_models\May_2019\Individual_model_2_non_responsive_states\challenge_fitter\results\',cyto_str,'\',exp_cyto,'NTCTESTING.csv'),x)
% 
% %siDCR
expData=[expTab.LPS1_siDCR expTab.LPS10_siDCR expTab.LPS100_siDCR expTab.LPS1000_siDCR]/100;
x=fun_sim_explore_nr_primary(cell_in_soup, startTime, stopTime, expData, plusminus,simulations);
csvwrite(append(str_path,'PhD\MATLAB\8_state_models\May_2019\Individual_model_2_non_responsive_states\challenge_fitter\results\',cyto_str,'\',exp_cyto,'_estimates_siDCR.csv'),x)
%csvwrite(append(str_path,'PhD\MATLAB\8_state_models\May_2019\Individual_model_2_non_responsive_states\challenge_fitter\results\',cyto_str,'\',exp_cyto,'siDCRTESTING.csv'),x)
