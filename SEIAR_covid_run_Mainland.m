%% Data From The Ministry of Health of Mexico with the number of infections and death by state (We retrieve Mexican data)
url_Iglobe=('https://raw.githubusercontent.com/UgoAvila/COVID-19.Mexico/master/Confirmed_Infected_Mexico_And_By_State.csv');
urlwrite(url_Iglobe,'Iglobe.csv');
I_globe=readtable('Iglobe.csv','ReadRowNames',true);
     date=datetime(I_globe{1,1:end},'Format','dd/MM/yy');
 
     url_Dglobe=('https://raw.githubusercontent.com/UgoAvila/COVID-19.Mexico/master/Confirmed_Death_Mexico_And_By_State.csv');
     urlwrite(url_Dglobe,'Dglobe.csv');
     D_globe=readtable('Dglobe.csv','ReadRowNames',true);
 
    
 %% We evalute first the outbreak in the largest city of Mexico, it's capital    
     C_exp_mx=str2double(I_globe{{'Nacional'},1:end});       %Cumulative number of infected cases
     D_exp_mx=str2double(D_globe{{'Nacional'},1:end});  %Number of deaths
     date=date(1:end);               %Day 1 corresponds to 12/03/20
     C_exp=C_exp_mx(1:end);
     D_exp=D_exp_mx(1:end);
     Acumul_exp=C_exp*8.6;                   %Approx. CUMULATIVE number of asymptomatic cases
     
%%
t=linspace(0,140,141);
%x0 contains the best fit parameters from the Mexican outbreak from (nuestro artículo). We also assume that 0<w<0.5, 0.01<p<0.20
x0=[0.4668 0.0100 27.8602 0 0.0829 11.2885 0.5 1./9.];
lb=[0.04,0,0,0,0,0,0,0.01];
ub=[1,0.01,30, 0.0541, 0.0541,100,0.5,0.20]; %the upper and the lower limits can affect the optimization results: physical based constraints should be used

N=127792286; % The population of the state of Mexico City, roughly
I0=C_exp(1);
RI0=0;
RA0=RI0;
E0=I0;
A0=Acumul_exp(1);

S0 = N-I0-RI0-RA0-E0-A0-D_exp(1);

%% Optimization section
% 1Â° Optimization with gradient based algorithm
options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'PlotFcn', { @optimplotfval });
[x,fval,exitflag,output,lambda,grad,hessian] =fmincon(@(x)SEIAR_covid_sse_Mainland(x,t,S0,I0,RI0,RA0,E0,A0,C_exp,D_exp,Acumul_exp,N),x0,[],[],[],[],lb,ub,[],options);

% Optimization with patternsearch
x0=x;
options = psoptimset;
options = psoptimset(options,'Display', 'off');
options = psoptimset(options,'PlotFcns', { @psplotbestf });
%options = psoptimset(options, 'MaxIter', 500);
[x,fval,exitflag,output] = patternsearch(@(x)SEIAR_covid_sse_Mainland(x,t,S0,I0,RI0,RA0,E0,A0,C_exp,D_exp,Acumul_exp,N),x0,[],[],[],[],lb,ub,[],options);

% 2nd Optimization with gradient based algorithm
x0=x;
options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'PlotFcn', { @optimplotfval });
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x)SEIAR_covid_sse_Mainland(x,t,S0,I0,RI0,RA0,E0,A0,C_exp,D_exp,Acumul_exp,N),x0,[],[],[],[],lb,ub,[],options);

%% Solving the model with the best fit parameters

[t,I,RI,D,RA,E,A]=SEIAR_covid_solver_Mainland(x,t,S0,I0,RI0,RA0,E0,A0,N);

%% Graphs
clf

%Plot the kinetics of the outbreak
for i=1:length(t)
        delta_t(i)=x(4)*exp( -t(i)/x(6) )+x(5);
        beta_t(i)=x(1)*exp( -t(i)/x(3) );
end

gamma=0.0828./(1+exp(-t+11.2885))+ 0;

data_long=dateshift(date(1),'start','day',0:t(end));

figure(2)
plot(data_long,I,'-r',data_long,RI,'-g')
hold on
ylabel('Infected and Recovered people')
hold on
yyaxis right
plot(data_long,D,'-k',date,D_exp,'ko',data_long,A,'-b')
 ylim([0 150000])
ylabel('Dead and Asymptomatic')
hold off
legend({'Infected','Recovered','Dead','Dead data','Asymptomatic'},'Location','northwest')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

figure(3)
plot(t,beta_t','LineWidth',2)
ylabel('Infection rate [1/day]')
xlabel('Time [day]')
hold on
yyaxis right
plot(t,delta_t','LineWidth',2)
hold on
plot(t,gamma','LineWidth',2)
ylabel('Death and Recovery rate [1/day]')
hold off
legend('Infection rate','Death rate', 'Recovery rate')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
grid on

figure(4)
p=x(8);
g=gamma';
Rdt = beta_t*p./(delta_t + g) + beta_t*(1-p)./g;
plot(data_long,1+t-t,'--','LineWidth',2)
ylim([0 16])
hold on
plot(data_long,Rdt,'-','LineWidth',2)
title('Daily reproduction number')
xlabel('days')
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
title('Infected')
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
 plot(date,RI([1:length(C_exp)]),'-k')
 title('Recovered')
 xlabel('days')
 ylabel('People')

 figure(10)
 plot(date,A([1:length(C_exp)]),'-')
 title('Asymptomatic cases')
 xlabel('days')
 ylabel('People')

 

AC = A + RA;
 figure(11)
 plot(data_long,AC,'-',date,Acumul_exp,'*')
 title('Cumulative number of asymptomatic cases')
 xlabel('days')
 ylabel('People')
 figure (12)
 data_long=dateshift(date(1),'start','day',0:t(end));
plot(data_long,I,'-r',data_long,RI,'-g')
hold on
ylabel('Infected and Recovered people')
hold on
yyaxis right
plot(data_long,D,'-k',date,D_exp,'ko',data_long,A,'-b')
ylim([0 10000])
ylabel('Dead and Asymptomatic')
hold off
legend({'Infected','Recovered','Dead','Dead data','Asymptomatic'}, 'Location','northwest')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
 