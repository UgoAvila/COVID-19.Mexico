%Step 3
% Plots the simulation curves of the model and the confidence intervals
data_long=dateshift(date(1),'start','day',0:length(t)-1);
load('par.mat','M')
load('qujian.mat','FangchaM')
load('dataestp_C.mat','dataestp_C')
load('dataestp_I.mat','dataestp_I')
load('dataestp_RI.mat','dataestp_RI')
load('dataestp_D.mat','dataestp_D')
load('dataestp_A.mat','dataestp_A')
load('dataestp_R0.mat','dataestp_R0')
load('dataestp_R0.mat','dataestp_R0')
load('dataestp_betaI.mat','dataestp_betaI')
load('dataestp_betaA.mat','dataestp_betaA')
load('dataestp_delta.mat','dataestp_delta')
load('dataestp_gamma.mat','dataestp_gamma')
C1=1.96*std(dataestp_C,0,2);
I1=1.96*std(dataestp_I,0,2);
RI1=1.96*std(dataestp_RI,0,2);
D1=1.96*std(dataestp_D,0,2);
A1=1.96*std(dataestp_A,0,2);
R01=1.96*std(dataestp_R0,0,2);

BetaI1=1.96*std(dataestp_betaI,0,2);
BetaA1=1.96*std(dataestp_betaA,0,2);
Delta1=1.96*std(dataestp_delta,0,2);
Gamma1=1.96*std(dataestp_gamma,0,2);

% figure(1)
% plot(data_long,estp_R0','-k',data_long,estp_R0'+R01,':r',data_long,estp_R0'-R01,'-.g','LineWidth',2)
% title('Daily reproduction number')
% xlabel('date')
% ylabel('R_d(t)')
% xlim([data_long(1) data_long(end)])
% ylim([0 8])
% legend('R0','95% CI upper estimate','95% CI lower estimate');
% set(gca,'FontSize',13)
% hold on
% plot(data_long,1+t-t,'--','LineWidth',2)
% 
% figure(2)
% plot(data_long,estp_C,'-k',data_long,estp_C+C1,':r',data_long,estp_C-C1,'-.g','LineWidth',3)
% hold on
% plot(date,C_exp,'*')
% title('Cumulative number of symptomatic cases')
% legend({'Best fit solutions','95% CI upper estimate','95% CI lower estimate','Reported data'},'Location','northwest')
% xlabel('Date')
% ylabel('People')
% set(gca,'FontSize',13)
% grid on
% 
% figure(3)
% plot(data_long,estp_D,'-k',data_long,estp_D+D1,':r',data_long,estp_D-D1,'-.g','LineWidth',3)
% hold on
% plot(date,D_exp,'*')
% ylim([0 3.5e5])
% title('Death toll (D(t))')
% legend({'Best fit solutions','95% CI upper estimate','95% CI lower estimate','Reported data'},'Location','northwest')
% xlabel('Date')
% ylabel('People')
% set(gca,'FontSize',13)
% grid on
% 
% figure(4)
% plot(data_long,estp_I,'-k',data_long,estp_I+I1,':r',data_long,estp_I-I1,'-.g','LineWidth',3)
% hold on
% plot(date,I_exp,'*')
% title('Infected (I(t))')
% legend({'Best fit solutions','95% CI upper estimate','95% CI lower estimate','Reported data'},'Location','northwest')
% xlabel('Date')
% ylabel('People')
% set(gca,'FontSize',13)
% grid on
% 
% figure(5)
% plot(data_long,estp_RI,'-k',data_long,estp_RI+RI1,':r',data_long,estp_RI-RI1,'-.g','LineWidth',3)
% hold on
% plot(date,R_exp,'*')
% title('Recovered (R_I(t))')
% legend({'Best fit solutions','95% CI upper estimate','95% CI lower estimate','Reported data'},'Location','northwest')
% xlabel('Date')
% ylabel('People')
% set(gca,'FontSize',13)
% grid on
% 
% figure(6)
% plot(data_long,estp_A,'-k',data_long,estp_A+A1,':r',data_long,estp_A-A1,'-.g','LineWidth',3)
% title('Asymptomatic (A(t))')
% legend({'Best fit solutions','95% CI upper estimate','95% CI lower estimate'},'Location','northwest')
% xlabel('Date')
% ylabel('People')
% set(gca,'FontSize',13)
% grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(7)
plot(t,estp_betaI','-k',t,estp_betaI'+BetaI1,':r',t,estp_betaI'-BetaI1,'-.g','LineWidth',3)
ylim([0 0.9])
title('\beta_I(t)')
legend({'Best fit rate','95% CI upper estimate','95% CI lower estimate'},'Location','northwest')
xlabel('Days')
ylabel('Infection rate (1/day)')
set(gca,'FontSize',13)
grid on

figure(8)
plot(t,estp_betaA','-k',t,estp_betaA'+BetaA1,':r',t,estp_betaA'-BetaA1,'-.g','LineWidth',3)
ylim([0 0.9])
title('\beta_A(t)')
legend({'Best fit rate','95% CI upper estimate','95% CI lower estimate'},'Location','northwest')
xlabel('Days')
ylabel('Infection rate (1/day)')
set(gca,'FontSize',13)
grid on

figure(9)
plot(t,estp_delta','-k',t,estp_delta'+Delta1,':r',t,estp_delta'-Delta1,'-.g','LineWidth',3)
title('\delta(t)')
legend({'Best fit rate','95% CI upper estimate','95% CI lower estimate'},'Location','northwest')
xlabel('Days')
ylabel('Death rate (1/day)')
set(gca,'FontSize',13)
grid on

figure(10)
plot(t,estp_gamma','-k',t,estp_gamma'+Gamma1,':r',t,estp_gamma'-Gamma1,'-.g','LineWidth',3)
title('\gamma(t)')
legend({'Best fit rate','95% CI upper estimate','95% CI lower estimate'},'Location','northwest')
xlabel('Days')
ylabel('Recovery rate (1/day)')
set(gca,'FontSize',13)
grid on