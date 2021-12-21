function [filtered_params] = fit_2_challenge_time_points_param_pies(cytokine, cytokine_formatted, type, search, minimum_deviation)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%Initialise%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    startTime=0;
    stopTime=24;
    %number of cells in various states cells[2]=negative cells
    cells=[0 500 0 0];
    
    %
    LPS_small_first_dose=10;
    LPS_first_dose=1000;
    LPS_second_dose=1000;
    %load experimental values
    %add your csv file here with percentage positives with different time
    %points
    E2 = readtable(strcat('bmdm1_',cytokine,'_',type,'.csv'),'ReadVariableNames',true);
    
    
    if (search==1)
        hyper_combinations = [0.0 0.0005 0.0005 0.001  0.005 0.05];
        param_search_mat_latin_hypercube = fun_latin_hypercube(hyper_combinations, 1);
        param_search_mat = param_search_mat_latin_hypercube;
        param_search_mat_norm = abs(normrnd(3,10,[500000,5]))./100;
        param_search_mat_negbinom = nbinrnd(0.7,0.05,55000,5)./100;
        param_search_mat =  [param_search_mat; param_search_mat_norm; param_search_mat_negbinom];
    %use search=0 if you have a pre-saved list of parameters to test against
    else
        path2read='D:\GoogleDrive\silence\PAPER\Macrophage AIH\FRONTIERS_review\modelling\NoRM\results_with_mu\';                                        
        param_search_mat = readtable(strcat(path2read,'\',cytokine_formatted,'\','exp_',cytokine,'_',type,'.csv'),'ReadVariableNames',true);
        param_search_mat = param_search_mat(:,1:5);
        param_search_mat = table2array(param_search_mat);
    end

    %param_search_mat =  [nbinrnd(0.7,0.05,3000,5)./100];
    % end    
    operatingMatrix = param_search_mat;
    %QC step. beta1 and gamma1 will be zero if gamma=0.
    ind1 = operatingMatrix(:,3)==0; %grab indices along third column where value is zero
    operatingMatrix(ind1,4) = 0; %set gamma2 for those indices to 0
    operatingMatrix(ind1,5) = 0; %set gamma2 for those indices to 0
    %sorted_param_search=sortrows(param_search_mat,3);
    [rows, columns]=size(operatingMatrix);
    final_params=[];
    primary_community=[];
    primary_16_community=[];
    secondary_community=[];
    plotAllEvolutionMat=[];
    plotLowAllEvolutionMat=[];
    bigChallengeMat=[];
    bigReChallengeMat=[];
    bigLowChallengeMat=[];
    bigReLowChallengeMat=[];
    %The following is dependedent on currently hard-coded column names.
    %Change to generic in next release
    %%challenge
    %8hr
    e2chall4 =  E2.M_1000_8/100;
    %12hr
    e2chall8 =  E2.M_1000_12/100;
    %16hr
    e2chall12 =  E2.M_1000_16/100;

    %%twice-challenged - 1 0/1000
    %8hr
    e2loTol4 =  E2.x10_1000_8/100;
    %12hr
    e2loTol8 =  E2.x10_1000_12/100;
    %16hr
    e2loTol12 =  E2.x10_1000_16/100;

    %%twice-challenged - 1000/1000
    %8hr
    e2Tol4 =  E2.x1000_1000_8/100;
    %12hr
    e2Tol8 =  E2.x1000_1000_12/100;
    %16hr
    e2Tol12 =  E2.x1000_1000_16/100;

    sendToPlottingFunction = [e2chall4, e2chall8, e2chall12; e2loTol4, e2loTol8, e2loTol12; e2Tol4, e2Tol8, e2Tol12];

    total_time_in_hours = 0:47;

    for row=1:rows

        %[1. p_minus_TO_p_plus 2. p_plus_TO_p_minus 3. p_plus_TO_p_nr 4. p_plus_TO_p_nr_minus 5. p_nr_TO_p_minus]
        prot_rate_coeff=operatingMatrix(row,:);
        low_LPS=LPS_small_first_dose;
        LPS=LPS_first_dose;
        %Gillespie 1: For low dose and high dose
        [prot_lps_timeCommunity, challenge_lps_vector_time] = Gillespie_4_state_5_rate_memory(cells, prot_rate_coeff, LPS, startTime, stopTime);
        [prot_low_lps_timeCommunity, challenge_low_lps_vector_time] = Gillespie_4_state_5_rate_memory(cells, prot_rate_coeff, low_LPS, startTime, stopTime);

        challenged_species=prot_lps_timeCommunity(1,:);
        low_challenged_species=prot_low_lps_timeCommunity(1,:);

        %Take first 24 items
        challenged_species_proportion=challenged_species(1:stopTime)./sum(cells);
        low_challenged_species_proportion=low_challenged_species(1:stopTime)./sum(cells);
        %%%%%%%%%%%%%%%%BLOCK2 CHALLENGE 2%%%%%%%%%%%%%%%%%%%
        %take the snapshot of the challenge community at 24hours%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        LPS=LPS_second_dose;
        %pick current state of community after first dose
        primary_community_at_16_hours=transpose(prot_lps_timeCommunity(:,16));
        new_LPS_community=transpose(prot_lps_timeCommunity(:,stopTime));
        new_low_LPS_community=transpose(prot_low_lps_timeCommunity(:,stopTime));
        %Gillespie 2: High dose only
        [prot_re_lps_timeCommunity, rechallenge_lps_vector_time] = Gillespie_4_state_5_rate_memory(new_LPS_community, prot_rate_coeff, LPS, startTime, stopTime);
        [prot_low_re_lps_timeCommunity, rechallenge_low_lps_vector_time] = Gillespie_4_state_5_rate_memory(new_low_LPS_community, prot_rate_coeff, LPS, startTime, stopTime);

        rechallenged_species=prot_re_lps_timeCommunity(1,:);
        low_rechallenged_species=prot_low_re_lps_timeCommunity(1,:);
        %Take 24 items - discard first item beause that is the memory
        rechallenged_species_proportion=rechallenged_species(2:stopTime+1)./sum(cells);
        low_rechallenged_species_proportion=low_rechallenged_species(2:stopTime+1)./sum(cells);
        %
        plotAllEvolution = [challenged_species_proportion rechallenged_species_proportion];
        plotLowAllEvolution = [low_challenged_species_proportion low_rechallenged_species_proportion];


        %calculating difference between experimental and simulated value     
        mean_squared_sum = (challenged_species_proportion(:,4) - e2chall4)^2;
        mean_squared_sum = mean_squared_sum + (challenged_species_proportion(:,8) - e2chall8)^2;
        mean_squared_sum = mean_squared_sum + (challenged_species_proportion(:,12) - e2chall12)^2;

        mean_squared_sum = mean_squared_sum + (low_rechallenged_species_proportion(:,4) - e2loTol4)^2;
        mean_squared_sum = mean_squared_sum + (low_rechallenged_species_proportion(:,8) - e2loTol8)^2;
        mean_squared_sum = mean_squared_sum + (low_rechallenged_species_proportion(:,12) - e2loTol12)^2;

        mean_squared_sum = mean_squared_sum + (rechallenged_species_proportion(:,4) - e2Tol4)^2;
        mean_squared_sum = mean_squared_sum + (rechallenged_species_proportion(:,8) - e2Tol8)^2;
        mean_squared_sum = mean_squared_sum + (rechallenged_species_proportion(:,12) - e2Tol12)^2;

        root_mean_sq_deviation = (mean_squared_sum/9)^0.5;
        
        %selecting solutions that have less than 0.5 rmsd
        if (root_mean_sq_deviation <0.5) 
            %plot(total_time_in_hours, plotAllEvolution, colorMe)
            plotAllEvolutionMat=cat(1,plotAllEvolutionMat,plotAllEvolution);
            plotLowAllEvolutionMat=cat(1,plotLowAllEvolutionMat,plotLowAllEvolution);
            %add deviation score to final params
            prot_rate_coeff_score = [prot_rate_coeff, root_mean_sq_deviation];
            final_params=cat(1,final_params,prot_rate_coeff_score);
            %I'm making pies
            primary_16_community=cat(1,primary_16_community,primary_community_at_16_hours);
            primary_community=cat(1,primary_community,new_LPS_community);
            secondary_community=cat(1,secondary_community,transpose(prot_re_lps_timeCommunity(:,stopTime)));
            %For descriptive stats
            bigChallengeMat=cat(1,bigChallengeMat,challenged_species_proportion);
            bigReChallengeMat=cat(1,bigReChallengeMat, rechallenged_species_proportion);

            bigLowChallengeMat=cat(1,bigLowChallengeMat, low_challenged_species_proportion);
            bigReLowChallengeMat=cat(1,bigReLowChallengeMat, low_rechallenged_species_proportion);
        end

    end
    meta_data_for_file_name = strcat('exp_',cytokine,'_',type,'_gating');
    filterValue=minimum_deviation;
    filtered_params = fun_plot_fit_data(final_params, plotAllEvolutionMat, bigChallengeMat, bigReChallengeMat,...
                                plotLowAllEvolutionMat, bigLowChallengeMat, bigReLowChallengeMat,...
                                primary_16_community, secondary_community,...
                                cytokine_formatted, total_time_in_hours, sendToPlottingFunction,...
                                meta_data_for_file_name, search, filterValue);
end


