function [filtered_params] = fit_2_challenge_time_points_param_pies_test(cytokine, cytokine_formatted, type, states,  search, minimum_deviation,tpts)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%Initialise%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    startTime=0;
    stopTime=24;
    sample_pri_sec_challenge_at=tpts(3);
    %number of cells in various states cells[2]=negative cells
    cells=[0 5000 0 0];
    
    %
    LPS_small_first_dose=10;
    LPS_first_dose=1000;
    LPS_second_dose=1000;
    %load experimental values
    %add your csv file here with percentage positives with different time
    %points
    E2 = readtable(strcat(type,'_',cytokine,'.csv'),'ReadVariableNames',true);

    
    if (search==1)
        hyper_combinations = [0.0 0.0005 0.001 0.01 0.1 ];
        param_search_mat_latin_hypercube = fun_latin_hypercube(hyper_combinations, 10);
        %param_search_mat = param_search_mat_latin_hypercube;
        param_search_mat_norm = abs(normrnd(3,10,[10000,5]))./100;
        param_search_mat_negbinom = nbinrnd(0.7,0.05,10000,5)./100;
        param_search_mat =  [param_search_mat_norm; param_search_mat_negbinom];
        %search LPS co-efficients mu
        mu_lps=[0.1 1 10 100];
    %use search=0 if you have a pre-saved list of parameters to test against
    else
        path2read='D:\GoogleDrive\silence\PAPER\Macrophage AIH\FRONTIERS_review\modelling\NoRM\results_june21';                                        
        param_search_mat = readtable(strcat(path2read,'\',cytokine_formatted,'\','exp_',cytokine,'_',type,'_',num2str(states),'_states.csv'),'ReadVariableNames',true);
        %keep LPS co-efficients mu as a separate matrix
        mu_lps=table2array(param_search_mat(:,6));
        param_search_mat = param_search_mat(:,1:5);
        param_search_mat = table2array(param_search_mat);
    end

    
    [num_param,~]=size(param_search_mat);
    
    if (search==1)
        num_param*length(mu_lps)
    else
        size(mu_lps)
    end
    
    %param_search_mat =  [nbinrnd(0.7,0.05,3000,5)./100];
    % end    
    operatingMatrix = param_search_mat;
    %QC step. beta1 and gamma1 will be zero if gamma=0.
    if (states>3)
      ind1 = operatingMatrix(:,3)==0; %grab indices along third column where value is zero
      operatingMatrix(ind1,4) = 0; %set gamma2 for those indices to 0
      operatingMatrix(ind1,5) = 0; %set beta2 for those indices to 0 
    elseif (states==3)  %
      operatingMatrix(:,4) = 0;   
    end
   
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
    
    [~, mu_columns]=size(mu_lps);
    
    counter=0;
    %
        
    %looping through parameters
    for mu = 1:mu_columns
        %looping through mu_lps values
        for row=1:rows

            %[1. p_minus_TO_p_plus 2. p_plus_TO_p_minus 3. p_plus_TO_p_nr 4. p_plus_TO_p_nr_minus 5. p_nr_TO_p_minus]
            prot_rate_coeff=operatingMatrix(row,1:5);
            low_LPS=LPS_small_first_dose;
            LPS=LPS_first_dose;
            if (search==1)
                mu_val=mu_lps(1,mu);
            else
                mu_val=mu_lps(row,1);
            end
            %Gillespie 1: For low dose and high dose
            [prot_lps_timeCommunity, challenge_lps_vector_time] = Gillespie_4_state_5_rate_memory(cells, prot_rate_coeff, LPS, startTime, stopTime, mu_val);
            [prot_low_lps_timeCommunity, challenge_low_lps_vector_time] = Gillespie_4_state_5_rate_memory(cells, prot_rate_coeff, low_LPS, startTime, stopTime, mu_val);

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
            primary_community_at_16_hours=transpose(prot_lps_timeCommunity(:,sample_pri_sec_challenge_at));
            new_LPS_community=transpose(prot_lps_timeCommunity(:,stopTime));
            new_low_LPS_community=transpose(prot_low_lps_timeCommunity(:,stopTime));
            %Gillespie 2: High dose only
            [prot_re_lps_timeCommunity, rechallenge_lps_vector_time] = Gillespie_4_state_5_rate_memory(new_LPS_community, prot_rate_coeff, LPS, startTime, stopTime, mu_val);
            [prot_low_re_lps_timeCommunity, rechallenge_low_lps_vector_time] = Gillespie_4_state_5_rate_memory(new_low_LPS_community, prot_rate_coeff, LPS, startTime, stopTime, mu_val);

            rechallenged_species=prot_re_lps_timeCommunity(1,:);
            low_rechallenged_species=prot_low_re_lps_timeCommunity(1,:);
            %Take 24 items - discard first item beause that is the memory
            rechallenged_species_proportion=rechallenged_species(2:stopTime+1)./sum(cells);
            low_rechallenged_species_proportion=low_rechallenged_species(2:stopTime+1)./sum(cells);
            %
            plotAllEvolution = [challenged_species_proportion rechallenged_species_proportion];
            plotLowAllEvolution = [low_challenged_species_proportion low_rechallenged_species_proportion];


            %calculating difference between experimental and simulated value     
            mean_squared_sum = (challenged_species_proportion(:,tpts(1)) - e2chall4)^2;
            mean_squared_sum = mean_squared_sum + (challenged_species_proportion(:,tpts(2)) - e2chall8)^2;
            mean_squared_sum = mean_squared_sum + (challenged_species_proportion(:,tpts(3)) - e2chall12)^2;

            mean_squared_sum = mean_squared_sum + (low_rechallenged_species_proportion(:,tpts(1)) - e2loTol4)^2;
            mean_squared_sum = mean_squared_sum + (low_rechallenged_species_proportion(:,tpts(2)) - e2loTol8)^2;
            mean_squared_sum = mean_squared_sum + (low_rechallenged_species_proportion(:,tpts(3)) - e2loTol12)^2;

            mean_squared_sum = mean_squared_sum + (rechallenged_species_proportion(:,tpts(1)) - e2Tol4)^2;
            mean_squared_sum = mean_squared_sum + (rechallenged_species_proportion(:,tpts(2)) - e2Tol8)^2;
            mean_squared_sum = mean_squared_sum + (rechallenged_species_proportion(:,tpts(3)) - e2Tol12)^2;

            root_mean_sq_deviation = (mean_squared_sum/9)^0.5;
            
            %AIC = N*log(RSS/N) + 2k 
            if (states>3)
                k=6;
            else
                k=5;
            end    
            aic = (9*log((mean_squared_sum/9))) + (2*k);
            %selecting solutions that have less than 0.5 rmsd
            if (root_mean_sq_deviation < 0.3) 
                %plot(total_time_in_hours, plotAllEvolution, colorMe)
                plotAllEvolutionMat=cat(1,plotAllEvolutionMat,plotAllEvolution);
                plotLowAllEvolutionMat=cat(1,plotLowAllEvolutionMat,plotLowAllEvolution);
                %add lps_mu, deviation score to final params 
                prot_rate_coeff_score = [prot_rate_coeff, mu_lps(mu), root_mean_sq_deviation, aic];
                final_params=cat(1,final_params,prot_rate_coeff_score);
                %I'm making pies
                primary_16_community=cat(1,primary_16_community,primary_community_at_16_hours);
                primary_community=cat(1,primary_community,new_LPS_community);
                %change time points based on cell type kinetics
                secondary_community=cat(1,secondary_community,transpose(prot_re_lps_timeCommunity(:,sample_pri_sec_challenge_at)));
                %For descriptive stats
                bigChallengeMat=cat(1,bigChallengeMat,challenged_species_proportion);
                bigReChallengeMat=cat(1,bigReChallengeMat, rechallenged_species_proportion);

                bigLowChallengeMat=cat(1,bigLowChallengeMat, low_challenged_species_proportion);
                bigReLowChallengeMat=cat(1,bigReLowChallengeMat, low_rechallenged_species_proportion);
         
            end
            counter=counter+1;
            if counter/rows > 0.49 && counter/rows < 0.50
                fprintf('50pc complete - %s, states - %d\n',cytokine_formatted, states)
            end 
        end
    end    
    if (isempty(final_params)==false)
        meta_data_for_file_name = strcat('exp_',cytokine,'_',type,'_',num2str(states),'_states');
        filterValue=minimum_deviation;
        filtered_params = fun_plot_fit_data(final_params, plotAllEvolutionMat, bigChallengeMat, bigReChallengeMat,...
                                    plotLowAllEvolutionMat, bigLowChallengeMat, bigReLowChallengeMat,...
                                    primary_16_community, secondary_community,...
                                    cytokine_formatted, total_time_in_hours, sendToPlottingFunction,...
                                                                    meta_data_for_file_name, search, filterValue, tpts);

    end
end

