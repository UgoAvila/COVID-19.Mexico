function dydt = SEIAR_covid_odes_Mainland(t,y,beta0,beta1,tau_beta,delta0,delta1,tau_delta,N,w,p)

ld = 0;         %We assume that beta and delta are decreasing since day 0
 if t<=ld
     beta=beta0;
     delta=delta0;
 else
    beta=beta0*exp((-t+ld)/(tau_beta))+beta1;
    delta=delta0*exp((-t+ld)/tau_delta)+delta1;

 end

gamma=0.0828./(1+exp(-t+11.2885))+ 0;

S=y(1);
I=y(2);
RI=y(3);
RA=y(4);
E=y(5);
A=y(6);

dS_dt = -beta*S*( (I+A)/(S+E+I+A+RI+RA) );
dE_dt = beta*S*( (I+A)/(S+E+I+A+RI+RA) ) - w*E;
dI_dt = p*w*E - gamma*I-delta*I;
dA_dt = (1-p)*w*E - gamma*A;
dRI_dt = gamma*I;               %Recovered from symptomatic infectious
dRA_dt = gamma*A;               %Recovered from asymptomatic
dydt = [dS_dt; dI_dt; dRI_dt; dRA_dt; dE_dt; dA_dt];

end