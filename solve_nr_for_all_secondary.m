%%drive path
str_home='C:\Users\shoum\Google Drive\silence\';
%str_work='D:\GoogleDrive\silence\';
%%The common parameters
cell_in_soup=[0 0 100 0];
plusminus=45;
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
expTab=readtable(append(str_home,'PhD\MATLAB\8_state_models\May_2019\Individual_model_2_non_responsive_states\challenge_fitter\',exp_cyto,'_2ndary.csv'),'ReadVariableNames',true);

%NoT 
expData=[expTab.M_1000_NoT expTab.x10_1000_NoT expTab.x10_1000_NoT]/100;
% x=fun_sim_explore_nr_primary(cell_in_soup, expData, plusminus,simulations);
% csvwrite(append(str_home,'PhD\MATLAB\8_state_models\May_2019\Individual_model_2_non_responsive_states\challenge_fitter\results\',cyto_str,'\',exp_cyto,'_2ndary_estimates_NoT.csv'),x)
%NTC
expData=[expTab.M_1000_NTC expTab.x10_1000_NTC expTab.x10_1000_NTC]/100; 
% x=fun_sim_explore_nr_primary(cell_in_soup, expData, plusminus,simulations);
% csvwrite(append(str_home,'PhD\MATLAB\8_state_models\May_2019\Individual_model_2_non_responsive_states\challenge_fitter\results\',cyto_str,'\',exp_cyto,'_2ndary_estimates_NTC.csv'),x)
% %siDCR
expData=[expTab.M_1000_siDCR expTab.x10_1000_siDCR expTab.x10_1000_siDCR]/100;
% x=fun_sim_explore_nr_primary(cell_in_soup, expData, plusminus,simulations);
% csvwrite(append(str_home,'PhD\MATLAB\8_state_models\May_2019\Individual_model_2_non_responsive_states\challenge_fitter\results\',cyto_str,'\',exp_cyto,'_2ndary_estimates_siDCR.csv'),x)
