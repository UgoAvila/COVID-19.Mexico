function [t,I,RI,D,RA,E,A]=SEIAR_covid_solver_mx(x,t,S0,I0,RI0,RA0,E0,A0,N)
% This function computes the solutions of the model with parameter vector
% x, initial values (S0,I0,RI0,RA0,E0,A0), and total population size N,
% over the time range t
beta0I=x(1);
beta1I=x(2);
tau_beta=x(3);
delta0=x(4);
delta1=x(5);
tau_delta=x(6);
gamma0=x(7);
gamma1=x(8);
tau_gamma=x(9);
beta0A=x(10);
beta1A=x(11);

%% Solving the model
options=odeset('NonNegative',(1:6));
[t,y] = ode45(@(t,y)SEIAR_covid_odes_mx(t,y,beta0I,beta1I,tau_beta,...
    delta0,delta1,tau_delta,gamma0,gamma1,tau_gamma,beta0A,beta1A),t,[S0; I0; RI0; RA0; E0; A0],options);
death_mod=N-y(:,1)-y(:,2)-y(:,3)-y(:,4)-y(:,5)-y(:,6);
I=y(:,2);
RI=y(:,3);
D=death_mod;
RA=y(:,4);
E=y(:,5);
A=y(:,6);
end