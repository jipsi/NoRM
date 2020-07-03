%%Model basis/assumptions:Cells make inflammatory proteins upon LPS stimulus. Single cells respond
%%switching on protein production (p_plus) and switching it off or by stopping the
%%response by not making anymore. This can be either by a stochastic
%%switching back to a negative state(p_minus) or a switch to a
%%non-responsive state (p_nr)
%%which can then be the basis for hypo-responses seen in double LPS
%%challenge etc. These non-responsive cells can then further become
%%non-responsive forever (p_nr_minus) or non-responsive state (p_nr) 
function [vector_community_time_evolution, vector_time] = Gillespie_4_state_5_rate_memory(cell_in_soup, rate_coeff, LPS, startTime, stopTime)
    %%
    %Override combinations here
    time_to_finish = stopTime; %Gillespie reactions run for x hours 
    %initialise individual cell matrix
    %macrophage_matrix=mac_mat;
    
    %total cells
    total_cells = sum(cell_in_soup);
    %Set history
    p_plus=cell_in_soup(1);    
    p_minus=cell_in_soup(2);
    p_nr=cell_in_soup(3);
    p_nr_minus=cell_in_soup(4);   

    time=startTime;
    vector_p_minus = zeros(1,time_to_finish+1);
    vector_p_plus = zeros(1,time_to_finish+1);
    vector_p_nr = zeros(1,time_to_finish+1);
    vector_p_nr_minus = zeros(1,time_to_finish+1);
    vector_time = zeros(1,time_to_finish+1);
    %initiate time index
    i=0;
    %initialise vectors
    vector_p_minus(1,i+1) = p_minus;
    vector_p_plus(1,i+1) = p_plus;
    vector_p_nr(1,i+1) = p_nr;
    vector_p_nr_minus(1,i+1) = p_nr_minus;
    %[alpha. p_minus_TO_p_plus beta. p_plus_TO_p_minus gamma. p_plus_TO_p_nr delta. p_plus_TO_p_nr_minus beta2. p_nr_TO_p_minus]
    %co-efficient affecting rate of appearance of TNF+ cells
    alpha=rate_coeff(1);
    %rate of disappearance of p+ cells
    beta=rate_coeff(2);
    %rate of switching to p_nr
    gamma=rate_coeff(3);
    %rate of switching to p_nr_minus
    gamma2=rate_coeff(4);
    %rate of switching to p_minus from p_nr
    beta2=rate_coeff(5);
    %rate at which LPS degrades
    lps_decay_rate=0.5;
    %LPS dose at start
    LPS_deltaT = LPS;
    %effect of lps
    lps_mu=100;
    %compute reaction rates (propensity*number of species)
    while time <= time_to_finish+5

      rate_p_plus = alpha*(1+(lps_mu*LPS_deltaT));
       rate_p_minus = beta;
       rate_p_nr = gamma;
       rate_p_nr_minus = gamma2;
       rate_p_minus_from_nr = beta2;
       
       %propensity = [1. p_minus_TO_p_plus 2. p_plus_TO_p_minus 3. p_plus_TO_p_nr 4. p_plus_TO_p_nr_minus 5. p_nr_TO_p_minus] 
       propensity = [rate_p_plus*p_minus rate_p_minus*p_plus rate_p_nr*p_plus rate_p_nr_minus*p_nr rate_p_minus_from_nr*p_nr];
       %Combined reaction Ao or propensity
       propAo = sum(propensity);
       %%Time_jump
       random1 = rand;
       while random1==0;
           random1 = rand;
       end
       %Calculate time_to_reaction 
       time_to_reaction = (1/propAo)*(log(1/random1));
       
       %time after reaction
       lastTime = time;
       time = time + time_to_reaction;
       deltaT = time - lastTime;

       %%Determine next reaction
       counter = 1; Mu = 0; Amu = 0; random2 = rand;
       while Amu < random2*propAo,
           Mu = Mu + 1;
           Amu = Amu + propensity(counter);
           counter = counter + 1;   
       end
       %events

        % draw a random row from replacement
%         mac_index = randsample(1:length(macrophage_matrix), 1);
%         individual_mac = macrophage_matrix(mac_index,:);
        %Reaction Key: %x_t;t_x;t_t6;t6_t;t_tn;tn_t;
        if Mu == 1 %p_minus->p_plus
            p_plus = p_plus + 1;
            p_minus = p_minus - 1;
        elseif Mu == 2 %p_plus->p_minus
            p_plus = p_plus - 1;
            p_minus = p_minus + 1;
        elseif Mu == 3 %p_plus->p_nr    
            p_plus = p_plus - 1;
            p_nr = p_nr + 1;
        elseif Mu == 4 %p_nr->p_nr_minus  
            p_nr = p_nr - 1;
            p_nr_minus = p_nr_minus + 1;
        elseif Mu == 5 %p_nr->p_minus  
            p_nr = p_nr - 1; 
            p_minus = p_minus + 1;
        end
        %store/output time and species
        if time >= i
            i=i+1;
            vector_p_minus(1,i) = p_minus;
            vector_p_plus(1,i) = p_plus;
            vector_p_nr(1,i) = p_nr;
            vector_p_nr_minus(1,i) = p_nr_minus;
            vector_time(1,i) = i;
        end
        LPS_deltaT = LPS_deltaT - lps_decay_rate*time_to_reaction*LPS_deltaT;
        if LPS_deltaT < 0
            LPS_deltaT = 0;
        end
    end
    %community arranged with TNF+ve cells first and in that order.
    %vector_community = [vector_p_plus(end) vector_p_minus(end) vector_p_nr(end) vector_p_nr_minus(end)];

    vector_community_time_evolution = [vector_p_plus; vector_p_minus; vector_p_nr; vector_p_nr_minus];
%     simCommunity = cat(1,simCommunity,vector_community_time_evolution);
%     matCommunity = cat(1,matCommunity,vector_community);

end
