rng('default')

% theta0=[0.05305;	0.15;	0.25;	0.01];
theta0=[0.05305;	0.15;	0.25;	0.01];

tspan = [0 24];

delta = 0.5;

c0=[0 500 0 0];
%first stim
c0(1)=1000;

mu=100; 

[t,y] = ode45(@(t,y) lps_dynamics_3state2(t,y,delta,theta0, mu), tspan, c0);

%plot(t,y(:,3),'-o')
% plot(t,y(:,2),'-o',t,y(:,3),'-o',t,y(:,4),'-o')

c0 = y(end,:);
%second stim
c0(1) = 1000;
tspan = [0 24];
[t1,y1] = ode45(@(t,y) lps_dynamics_3state2(t,y,delta,theta0, mu), tspan, c0);



figure(1)

plot([t; t1+24],[y(:,3); y1(:,3)],'-')
hold on
%plot(t1+24,y1(:,3),'-')
%hold off
% 
%hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cells=[0 500 0 0];
prot_rate_coeff=[0.05305;	0.15;	0.25;	0; 0.01];
low_LPS=1000;
LPS=1000;
LPS_second_dose=1000;
startTime=0;
stopTime=24;


%%%%%%%%%%%%%%%BLOCK1 CHALLENGE 1%%%%%%%%%%%%%%%%%%%
%%%%%Gillespie 1: For low dose and high dose
[prot_lps_timeCommunity, challenge_lps_vector_time] = Gillespie_4_state_5_rate_memory(cells, prot_rate_coeff, LPS, startTime, stopTime, mu);
%[prot_low_lps_timeCommunity, challenge_low_lps_vector_time] = Gillespie_4_state_5_rate_memory_test(cells, prot_rate_coeff, low_LPS, startTime, stopTime);

challenged_species=prot_lps_timeCommunity(1,:);
%low_challenged_species=prot_low_lps_timeCommunity(1,:);
%%%%%%%%%%%%%%%%BLOCK2 CHALLENGE 2%%%%%%%%%%%%%%%%%%%
%take the snapshot of the challenge community at 24hours%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LPS=LPS_second_dose;
%pick current state of community after first dose
%primary_community_at_16_hours=transpose(prot_lps_timeCommunity(:,16));
new_LPS_community=transpose(prot_lps_timeCommunity(:,stopTime));
%new_low_LPS_community=transpose(prot_low_lps_timeCommunity(:,stopTime));
%Gillespie 2: High dose only
[prot_re_lps_timeCommunity, rechallenge_lps_vector_time] = Gillespie_4_state_5_rate_memory(new_LPS_community, prot_rate_coeff, LPS, startTime, stopTime, mu);
%[prot_low_re_lps_timeCommunity, rechallenge_low_lps_vector_time] = Gillespie_4_state_5_rate_memory_test(new_low_LPS_community, prot_rate_coeff, LPS, startTime, stopTime);

rechallenged_species=prot_re_lps_timeCommunity(1,:);
%low_rechallenged_species=prot_low_re_lps_timeCommunity(1,:);
%figure(2)

% plot(0:24,challenged_species(1:25),'-.')
% hold on
% plot(24:47,rechallenged_species(2:25),'-.')
% hold on


%%%run it 5 times
for x=1:3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%BLOCK1 CHALLENGE 1%%%%%%%%%%%%%%%%%%%
    %Gillespie 1: For low dose and high dose
    [prot_lps_timeCommunity, challenge_lps_vector_time] = Gillespie_4_state_5_rate_memory_all_time(cells, prot_rate_coeff, LPS, startTime, stopTime, mu);
    %[prot_low_lps_timeCommunity, challenge_low_lps_vector_time] = Gillespie_4_state_5_rate_memory_test(cells, prot_rate_coeff, low_LPS, startTime, stopTime);

    challenged_species=prot_lps_timeCommunity(1,:);
    %low_challenged_species=prot_low_lps_timeCommunity(1,:);
    %%%%%%%%%%%%%%%%BLOCK2 CHALLENGE 2%%%%%%%%%%%%%%%%%%%
    %take the snapshot of the challenge community at 24hours%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LPS=LPS_second_dose;
    %pick current state of community after first dose
    new_LPS_community=transpose(prot_lps_timeCommunity(:,end));
    %Gillespie 2: High dose only
    [prot_re_lps_timeCommunity, rechallenge_lps_vector_time] = Gillespie_4_state_5_rate_memory_all_time(new_LPS_community, prot_rate_coeff, LPS, startTime, stopTime, mu);
    %[prot_low_re_lps_timeCommunity, rechallenge_low_lps_vector_time] = Gillespie_4_state_5_rate_memory_test(new_low_LPS_community, prot_rate_coeff, LPS, startTime, stopTime);

    rechallenged_species=prot_re_lps_timeCommunity(1,:);
    %low_rechallenged_species=prot_low_re_lps_timeCommunity(1,:);

    %figure(2)


    plot([challenge_lps_vector_time challenge_lps_vector_time(:,end)+rechallenge_lps_vector_time],[challenged_species rechallenged_species],'-.')
    hold on
%     plot(challenge_lps_vector_time(:,end)+rechallenge_lps_vector_time,rechallenged_species,'-.')
%     hold on
end
labels={'ODE 1st+2nd Stimulus','SSA_1 1st+2nd Stimulus','SSA_2 1st+2nd Stimulus','SSA_3 1st+2nd Stimulus'};
xlabel('time in hours')
ylabel('P+ cells (count)')
legend(labels)
hold off