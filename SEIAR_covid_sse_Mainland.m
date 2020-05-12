function sset=SEIAR_covid_sse_Mainland(x,t,S0,I0,RI0,RA0,E0,A0,C_exp,D_exp,Acumul_exp,N)

[t,I,RI,D,RA,E,A]=SEIAR_covid_solver_Mainland(x,t,S0,I0,RI0,RA0,E0,A0,N);

C_mod = I+RI+D;      %Cumulative infected, calculated by the model
Acumul_mod = A+RA;   %Cumulative asymptomatic, calculated by the model
sse_c=0;
sse_d=0;
sse_ac=0;
for i=1:length(C_exp)
    sse_c=sse_c+(C_mod(i)-C_exp(i))^2;
    sse_d=sse_d+(D(i)-D_exp(i))^2;
    sse_ac=sse_ac+(Acumul_mod(i)-Acumul_exp(i))^2;
end

sset = 10*sse_c + 10*sse_d + sse_ac; %The coefficients are used to compensate the order of magnitude of difference between the SSEs. Better normalization techniques could be used 


end