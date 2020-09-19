function sset=SEIAR_covid_sse_mx(x,t,S0,I0,RI0,RA0,E0,A0,C_exp,R_exp,D_exp,N)
% This function calculates the Sum of Squared Errors with the solutions of
% the SEIARD model with parameter vector x, initial conditions
% (S0,I0,RI0,RA0,E0,A0), experimental data (C_exp,R_exp,D_exp) and
% total population size N

[t,I,RI,D,RA,E,A]=SEIAR_covid_solver_mx(x,t,S0,I0,RI0,RA0,E0,A0,N);

Cmod = I+RI+D;  %Cumulative number of infectives (calculated with the model)
sse_c=0;
sse_r=0;
sse_d=0;
for i=1:length(C_exp)
    sse_c=sse_c+(Cmod(i)-C_exp(i))^2;
    sse_d=sse_d+(D(i)-D_exp(i))^2;
    if i>1
        if R_exp(i)>R_exp(i-1)
            sse_r=sse_r+(RI(i)-R_exp(i))^2;
        end
    end
end

sset = 20*sse_c + sse_r + 10*sse_d;
%The coefficients are used to compensate the order of magnitude of difference between the SSEs.


end