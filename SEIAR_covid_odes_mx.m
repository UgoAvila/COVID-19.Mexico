function dydt = SEIAR_covid_odes_mx(t,y,beta0I,beta1I,tau_beta,...
    delta0,delta1,tau_delta,gamma0,gamma1,tau_gamma,beta0A,beta1A)
%Equations of the model


betaI=beta0I*exp(-t/tau_beta)+beta1I;
betaA=beta0A*exp(-t/tau_beta)+beta1A;
delta=delta0*exp(-t/tau_delta)+delta1;
gamma=gamma1./(1+exp(-t+tau_gamma))+gamma0;

w=1/4.;
p=0.12;

S=y(1);
I=y(2);
RI=y(3);
RA=y(4);
E=y(5);
A=y(6);

dS_dt = -S.*(betaI*I + betaA*A)./(S+E+I+A+RI+RA);          %Susceptible
dE_dt = S.*(betaI*I + betaA*A)./(S+E+I+A+RI+RA) - w*E;     %Exposed
dI_dt = p*w*E - gamma*I-delta*I;                    %Infected (symptomatic)
dA_dt = (1-p)*w*E - gamma*A;                %Asymptomatic
dRI_dt = gamma*I;               %Recovered from symptomatic infection
dRA_dt = gamma*A;               %Recovered from asymptomatic infection
dydt = [dS_dt; dI_dt; dRI_dt; dRA_dt; dE_dt; dA_dt];
%Deaths are not included since they can be computed from the other variables.

end