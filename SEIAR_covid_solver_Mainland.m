function [t,I,RI,D,RA,E,A]=SEIAR_covid_solver_Mainland(x,t,S0,I0,RI0,RA0,E0,A0,N)

beta0=x(1);
beta1=x(2);
tau_beta=x(3);
delta0=x(4);
delta1=x(5);
tau_delta=x(6);
w=x(7);
p=x(8);

%% Solving the model
options=odeset('NonNegative',(1:6));
[t,y] = ode45(@(t,y)SEIAR_covid_odes_Mainland(t,y,beta0,beta1,tau_beta,delta0,delta1,tau_delta,N,w,p),t,[S0; I0; RI0; RA0; E0; A0],options);
death_mod=N-y(:,1)-y(:,2)-y(:,3)-y(:,4)-y(:,5)-y(:,6);
I=y(:,2);
RI=y(:,3);
D=death_mod;
RA=y(:,4);
E=y(:,5);
A=y(:,6);
end