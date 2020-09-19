%Step 2
% This code computes the solutions of the model for each of the 2000 sets
% of parameters (after the burn-in period) obtained in Step 1

load('Yu.mat','Yu')
p=0.12;
% t=linspace(0,360,361);
dataestp_C=zeros(length(t),NN-mm);
dataestp_I=dataestp_C; dataestp_RI=dataestp_C;
dataestp_D=dataestp_C; dataestp_A=dataestp_C;
dataestp_R0=dataestp_C;
dataestp_betaI=dataestp_C; dataestp_betaA=dataestp_C; dataestp_delta=dataestp_C; dataestp_gamma=dataestp_C;
for pp=1:NN-mm
    PARAS=Yu(pp,:);
    [t,I,RI,D,RA,E,A]=SEIAR_covid_solver_mx(PARAS,t,S0,I0,RI0,RA0,E0,A0,N);
    C = I + RI + D;
    betaI_t=zeros(1,length(t));
    delta_t=zeros(1,length(t));
    gamma_t=zeros(1,length(t));
    betaA_t=zeros(1,length(t));
    for i=1:length(t)
        betaI_t(i)=PARAS(1)*exp( -t(i)/x(3) ) + PARAS(2);
        delta_t(i)=PARAS(4)*exp( -t(i)/PARAS(6) ) + PARAS(5);
        gamma_t(i) = ( PARAS(8)./(1+exp(-t(i)+PARAS(9))) + PARAS(7) );
        betaA_t(i)=PARAS(10)*exp( -t(i)/x(3) ) + PARAS(11);
    end
    Rdt = betaI_t.*p./(delta_t + gamma_t) + betaA_t.*(1-p)./gamma_t;
    dataestp_C(:,pp)=C; dataestp_I(:,pp)=I; dataestp_RI(:,pp)=RI;
    dataestp_D(:,pp)=D; dataestp_A(:,pp)=A;
    dataestp_R0(:,pp)=Rdt;
    dataestp_betaI(:,pp)=betaI_t; dataestp_betaA(:,pp)=betaA_t;
    dataestp_delta(:,pp)=delta_t; dataestp_gamma(:,pp)=gamma_t;
end
estp_C=C;estp_I=I; estp_RI=RI; estp_D=D; estp_A=A; estp_R0=Rdt;
estp_betaI=betaI_t; estp_betaA=betaA_t;
estp_delta=delta_t; estp_gamma=gamma_t;
save('estp_C.mat','C')
save('estp_I.mat','I')
save('estp_RI.mat','RI')
save('estp_D.mat','D')
save('estp_A.mat','A')
save('estp_R0.mat','Rdt')
save('estp_betaI.mat','betaI_t')
save('estp_betaA.mat','betaA_t')
save('estp_delta.mat','delta_t')
save('estp_gamma.mat','gamma_t')

save('dataestp_C.mat','dataestp_C')
save('dataestp_I.mat','dataestp_I')
save('dataestp_RI.mat','dataestp_RI')
save('dataestp_D.mat','dataestp_D')
save('dataestp_A.mat','dataestp_A')
save('dataestp_R0.mat','dataestp_R0')
save('dataestp_betaI.mat','dataestp_betaI')
save('dataestp_betaA.mat','dataestp_betaA')
save('dataestp_delta.mat','dataestp_delta')
save('dataestp_gamma.mat','dataestp_gamma')