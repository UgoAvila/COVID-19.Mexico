% Data From Johns Hopkins CSSE
url_Cglobe=('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv');
urlwrite(url_Cglobe,'Cglobe.csv');
url_Dglobe=('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv');
urlwrite(url_Dglobe,'Dglobe.csv');
url_Rglobe=('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv');
urlwrite(url_Rglobe,'Rglobe.csv');

C_globe=readtable('Cglobe.csv','ReadRowNames',true);
date=datetime(C_globe{1,4:end},'Format','MM/dd/yy');

D_globe=readtable('Dglobe.csv','ReadRowNames',true);

R_globe=readtable('Rglobe.csv','ReadRowNames',true);

C_exp_mx=str2double(C_globe{{'Row179'},4:end});    %Row 179 corresponds to Mexico
D_exp_mx=str2double(D_globe{{'Row179'},4:end});
R_exp_mx=str2double(R_globe{{'Row166'},4:end});    %Row 166 corresponds to Mexico

initdate=38;                           %Date since which we want to fit the data
date=date(initdate:end);

C_exp=C_exp_mx(initdate:end);
D_exp=D_exp_mx(initdate:end);
R_exp=R_exp_mx(initdate:end);
I_exp=C_exp-R_exp-D_exp;
A_exp=(1-0.12)/0.12*I_exp;    %Considering I/(I+A) = 0.12
     
%%
t=linspace(0,length(C_exp)-1,length(C_exp));
%   beta0I  beta1I  tau_beta  delta0  delta1  tau_delta  gamma0  gamma1  tau_gamma  beta0A  beta1A
x0=[0.4646	0.2514  29.0893   0.0097  0.0355  6.4256     0.0895  0.1539  7.7042     0.4646	0.2514];
lb=[0.      0.      0.        0.      0.      0.         0.      0.      0.         0.      0.];
ub=[1.      0.2    200.      0.0541  0.0541  200.       1.      1.      200.       1.       0.2];

N=127090000; % The population of Mexico
%We set the initial conditions:
I0=I_exp(1);
RI0=R_exp(1);
RA0=RI0;
E0=I0;
A0=A_exp(1);

S0 = N-I0-RI0-RA0-E0-A0-D_exp(1);

%% Optimization section
% 1st Optimization with gradient based algorithm
options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'PlotFcn', { @optimplotfval });
[x,fval,exitflag,output,lambda,grad,hessian] =fmincon(@(x)SEIAR_covid_sse_mx(x,t,S0,I0,RI0,RA0,E0,A0,C_exp,R_exp,D_exp,N),x0,[],[],[],[],lb,ub,[],options);

% Optimization with patternsearch
x0=x;
options = psoptimset;
options = psoptimset(options,'Display', 'off');
options = psoptimset(options,'PlotFcns', { @psplotbestf });
%options = psoptimset(options, 'MaxIter', 500);
[x,fval,exitflag,output] = patternsearch(@(x)SEIAR_covid_sse_mx(x,t,S0,I0,RI0,RA0,E0,A0,C_exp,R_exp,D_exp,N),x0,[],[],[],[],lb,ub,[],options);

% 2nd Optimization with gradient based algorithm
x0=x;
options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'PlotFcn', { @optimplotfval });
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x)SEIAR_covid_sse_mx(x,t,S0,I0,RI0,RA0,E0,A0,C_exp,R_exp,D_exp,N),x0,[],[],[],[],lb,ub,[],options);

%% Solving the model with the best fit parameters

finaldate = 290;            %Date until which we want to make the predictions
t=linspace(0,finaldate,finaldate+1);

[t,I,RI,D,RA,E,A]=SEIAR_covid_solver_mx(x,t,S0,I0,RI0,RA0,E0,A0,N);

%% Graphs
clf

%Plot the kinetics of the outbreak
for i=1:length(t)
        delta_t(i)=x(4)*exp( -t(i)/x(6) ) + x(5);
        betaI_t(i)=x(1)*exp( -t(i)/x(3) ) + x(2);
        betaA_t(i)=x(10)*exp( -t(i)/x(3) ) + x(11);
end

gamma = ( x(8)./(1+exp(-t+x(9))) + x(7) )';
p=0.12;

data_long=dateshift(date(1),'start','day',0:t(end));

figure(12)
plot(data_long,I,':r','LineWidth',2)
hold on
plot(date,I_exp,'r*')
plot(data_long,D,'-k','LineWidth',2)
plot(date,D_exp,'ko')
ylabel('Infected and Dead')
hold on
yyaxis right
plot(data_long,RI,'-.g','LineWidth',2)
plot(date,R_exp,'gd')
plot(data_long,A,'--b','LineWidth',2)
% ylim([0 4e6])
ylabel('Recovered and Asymptomatic')
hold off
legend({'Infected','Infected data','Dead','Dead data','Recovered','Recovered data','Asymptomatic'},'Location','northwest')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

figure(13)
plot(data_long,I,':r','LineWidth',3)
hold on
plot(date,I_exp,'r*')
plot(data_long,D,'-k','LineWidth',3)
plot(date,D_exp,'ko')
xlabel('Date')
ylabel('People')
grid on
plot(data_long,RI,'-.g','LineWidth',3)
plot(date,R_exp,'gd')
plot(data_long,A,'--b','LineWidth',3)
ylim([0 1.2e6])
legend({'Infected','Infected data','Dead','Dead data','Recovered','Recovered data','Asymptomatic'},'Location','northwest')

figure(14)
plot(t,betaI_t','--b','LineWidth',4)
hold on
plot(t,betaA_t',':k','LineWidth',4)
xlabel('Days')
ylabel('Infection rate (1/day)')
legend('Symptomatic infection','Asymptomatic infection')
grid on

figure(3)
plot(t,delta_t',':r','LineWidth',4)
hold on
plot(t,gamma','-.g','LineWidth',4)
xlabel('Days')
ylabel('Death and Recovery rate (1/day)')
legend('Death rate', 'Recovery rate')
grid on

figure(4)
Rdt = betaI_t.*p./(delta_t + gamma) + betaA_t.*(1-p)./gamma;
plot(data_long,1+t-t,'--','LineWidth',2)
ylim([0 16])
hold on
plot(data_long,Rdt,'-','LineWidth',2)
title('Daily reproduction number')
xlabel('date')
ylabel('R_d(t)')

C = I + RI + D;
figure(5)
plot(data_long,C,'-','LineWidth',2)
hold on
plot(date,C_exp,'*')
title('Cumulative number of symptomatic cases')
legend({'Model solutions','Reported data'},'Location','northwest')
xlabel('days')
ylabel('People')
grid on

figure(6)
plot(data_long,D,'-','LineWidth',2)
hold on
plot(date,D_exp,'*')
title('Death toll')
legend({'Model solutions','Reported data'},'Location','northwest')
xlabel('days')
ylabel('People')
grid on

figure(7)
plot(data_long,I,'-k','LineWidth',2)
hold on
plot(date,I_exp,'*')
title('Infected')
legend({'Model solutions','Reported data'},'Location','northwest')
xlabel('days')
ylabel('People')
grid on

figure(8)
plot(data_long,I,'-','LineWidth',2)
hold on
plot(data_long,A,'--','LineWidth',2)
title('Infected')
legend({'Symptomatic (I(t))','Asymptomatic (A(t))'},'Location','northwest')
xlabel('days')
ylabel('People')
grid on

figure(9)
plot(data_long,RI,'-k','LineWidth',2)
hold on
plot(date,R_exp,'*')
title('Recovered')
legend({'Model solutions','Reported data'},'Location','northwest')
xlabel('days')
ylabel('People')
grid on

figure(10)
plot(date,A([1:length(C_exp)]),'-')
title('Asymptomatic cases')
xlabel('days')
ylabel('People')
grid on