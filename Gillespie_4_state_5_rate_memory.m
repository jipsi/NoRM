%%Model basis/assumptions:Cells make inflammatory proteins upon LPS stimulus. Single cells respond
%%switching on protein production (p_plus) and switching it off or by stopping the
%%response by not making anymore. This can be either by a stochastic
%%switching back to a negative state(p_minus) or a switch to a
%%non-responsive state (p_nr)
%%which forms the basis for hypo-responses seen in double LPS
%%challenge etc. These non-responsive cells can then further become
%%non-responsive permanently (p_nr_minus) or be non-responsive state (p_nr)
%%temporarily.
function [vector_community_time_evolution, vector_time] = Gillespie_4_state_5_rate_memory(cell_in_soup, rate_coeff, LPS, startTime, stopTime)
    %%
    %Override combinations here
    time_to_finish = stopTime; %Gillespie reactions run for x hours 
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
    %rate of appearance of p+ cells
    alpha=rate_coeff(1);
    %rate of disappearance of p+ cells
    beta=rate_coeff(2);
    %rate of switching to p_nr
    gamma=rate_coeff(3);
    %rate of switching to p_nr_minus
    delta=rate_coeff(4);
    %rate of switching to p_minus from p_nr
    beta2=rate_coeff(5);
    %rate at which LPS degrades
    lps_decay_rate=0.5;
    %LPS concentratoin at start
    LPS_deltaT = LPS;
    %LPS effect multiplier
    lambda = 100;
    %compute reaction rates (propensity*number of species)
    while time <= time_to_finish+5
       %Adding 1 then makes a small effect at 0 ng/ml
       rate_p_plus = alpha*(1+(lambda*LPS_deltaT));
       rate_p_minus = beta;
       rate_p_nr = gamma;
       rate_p_nr_minus = delta;
       rate_p_minus_from_nr = beta2;
       
       %propensity = [1. p_minus_TO_p_plus 2. p_plus_TO_p_minus 3. p_plus_TO_p_nr 4. p_plus_TO_p_nr_minus 5. p_nr_TO_p_minus] 
       propensity = [rate_p_plus*p_minus rate_p_minus*p_plus rate_p_nr*p_plus rate_p_nr_minus*p_nr rate_p_minus_from_nr*p_nr];
       %Combined reaction Ao or propensity
       propAo = sum(propensity);
       %%Time_jump
       random1 = rand;
       while random1==0
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
       while Amu < random2*propAo
           Mu = Mu + 1;
           Amu = Amu + propensity(counter);
           counter = counter + 1;   
       end
       %events
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
 
    vector_community_time_evolution = [vector_p_plus; vector_p_minus; vector_p_nr; vector_p_nr_minus];

end
