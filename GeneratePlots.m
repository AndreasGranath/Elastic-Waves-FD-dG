close all

figure
% Plot trigonometric errors
Hq=[1/20,1/30,1/40,1/50,1/60,1/70];
load("Trigonometric_EOC_4thOrder.mat")

loglog(Hq,Error(1:6),'-o','LineWidth',1.5,'Color','red')
hold on
loglog(Hq,1305.*Hq.^4,'-.','LineWidth',1.5,'Color','red')

load("Trigonometric_EOC_6thOrder.mat")
loglog(Hq,Error,'-x','LineWidth',1.5,'Color','blue')
hold on
loglog(Hq,10000.*Hq.^5,'--','LineWidth',1.5,'Color','blue')
legend("4th order","$h^4$","6th order","$h^5$","Interpreter","latex","Location","southeast")
xlabel("h","Interpreter","latex")
ylabel("Discrete error")
fontsize(gcf,20,"points")

%%
figure
% Plot Stoneley errors
Hq=[1/20,1/30,1/40,1/50,1/60,1/70];
load("Stoneley_EOC_4thOrder.mat")

loglog(Hq,Error,'-o','LineWidth',1.5,'Color','red')
hold on
loglog(Hq,4e4.*Hq.^4,'-.','LineWidth',1.5,'Color','red')

load("Stoneley_EOC_6thOrder.mat")
loglog(Hq,Error,'-x','LineWidth',1.5,'Color','blue')
hold on
loglog(Hq,4e5.*Hq.^5,'--','LineWidth',1.5,'Color','blue')
legend("4th order","$h^4$","6th order","$h^6$","Interpreter","latex","Location","southeast")
xlabel("h","Interpreter","latex")
ylabel("Discrete error")
fontsize(gcf,20,"points")
%%
figure
load("ShiftedEigenvalues.mat")
delta=[1e-6 5e-6 1e-5 5e-5 5e-5 1e-4 5e-4 1e-3 5e-3 1e-2 5e-2];
semilogx(delta,EIGS(1,1:end-1),'-x','LineWidth',1.5)
hold on
semilogx(delta,EIGS(2,1:end-1),'-o','LineWidth',1.5)
semilogx(delta,EIGS(3,1:end-1),'-^','LineWidth',1.5)
semilogx(delta,EIGS(4,1:end-1),'-v','LineWidth',1.5)
xlabel("Shift parameter $\delta$","Interpreter","latex")
ylabel("$\max_i|e_i|$","Interpreter","latex")
legend("$h=1/10$","$h=1/15$","$h=1/20$","$h=1/30$","Interpreter","latex","location","best")
fontsize(gcf,15,"points")