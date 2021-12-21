function dC=lps_dynamics_3state2(t,c,delta,theta, mu)
    dcdt=zeros(4,1);
    %theta definitions
    %alpha/theta(1): +P-N; beta(theta(2): -P+N
    %gamma/theta(3): +NR; beta2/theta(4): +P-NR
    %rate at which LPS degrades
    %rate of change of LPS in time
    dcdt(1)=-delta*c(1);
    %theta(1) value based on LPS
    alpha = theta(1)*(1+mu*(c(1)));
    %rate of change of species 'N' or dcdt(2)
    dcdt(2)= theta(2)*c(3) - alpha*c(2) + theta(4)*c(4);
    %rate of change of species 'P'
    %theta(1) is affected by LPS concentration
    dcdt(3)=alpha*c(2) - theta(2)*c(3)-theta(3)*c(3);
    %rate of change of species 'NR'
    dcdt(4)=theta(3)*c(3) -theta(4)*c(4);

    dC=dcdt;
end


