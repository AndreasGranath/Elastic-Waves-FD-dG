%% Add "Data" folder to path before running this code

close all
addpath('NewData\')
%% Plot discrete errors for the trigonometric experiment in a continuous material
figure

%Hq=[1/20,1/30,1/40,1/50,1/60,1/70];
load("ContinuousExperiment4thOrder.mat")

loglog(Hq(1:6),Error(1:6),'-o','LineWidth',1.5,'Color','red')
hold on
loglog(Hq(1:6),1305.*Hq(1:6).^4,'-.','LineWidth',1.5,'Color','red')

load("ContinuousExperiment6thOrder.mat")
loglog(Hq(1:6),Error(1:6),'-x','LineWidth',1.5,'Color','blue')
hold on
loglog(Hq(1:6),10000.*Hq(1:6).^5,'--','LineWidth',1.5,'Color','blue')
legend("4th order","$h^4$","6th order","$h^5$","Interpreter","latex","Location","southeast")
xlabel("$h$","Interpreter","latex")
ylabel("Discrete error")
fontsize(gcf,15,"points")
%% Plot discrete errors for the trigonometric experiment with curved boundary

figure

%Hq=[1/20,1/30,1/40,1/50,1/60,1/70];
load("CurvedBdryErrors4thOrder.mat")

loglog(Hq(1:6),Error(1:6),'-o','LineWidth',1.5,'Color','red')
hold on
loglog(Hq(1:6),90.*Hq(1:6).^4,'-.','LineWidth',1.5,'Color','red')

load("CurvedBdryErrors6thOrder.mat")
loglog(Hq(1:6),Error(1:6),'-x','LineWidth',1.5,'Color','blue')
hold on
loglog(Hq(1:6),200.*Hq(1:6).^5,'--','LineWidth',1.5,'Color','blue')
legend("4th order","$h^4$","6th order","$h^5$","Interpreter","latex","Location","southeast")
xlabel("$h$","Interpreter","latex")
ylabel("Discrete error")
fontsize(gcf,15,"points")
%% Plot discrete errors for the fine FD mesh example
figure

Hq=1./(10*[2:1:13]);
load("FineFDExperiment4thOrder.mat","Error")

loglog(Hq(7:12),Error(7:12),'-o','LineWidth',1.5,'Color','red')
hold on
loglog(Hq(7:12),15.*Hq(7:12).^4,'-.','LineWidth',1.5,'Color','red')

load("FineFDExperiment6thOrder.mat","Error")
loglog(Hq(7:12),Error(7:12),'-x','LineWidth',1.5,'Color','blue')
hold on
loglog(Hq(7:12),10.*Hq(7:12).^5,'--','LineWidth',1.5,'Color','blue')
legend("4th order","$h^4$","6th order","$h^5$","Interpreter","latex","Location","southeast")
xlabel("$h$","Interpreter","latex")
ylabel("Discrete error")
fontsize(gcf,15,"points")
%% Plot the number of DOFs in the FD and dG domains when FD and dG nodes align
load("ContinuousExperiment4thOrder.mat","SystemSizedG","SystemSizeFD")
Hq=1./(10*[2:1:13]);

loglog(Hq(7:12),SystemSizedG(7:12),'-o','LineWidth',1.5)
hold on
loglog(Hq(7:12),SystemSizeFD(7:12),'-x','LineWidth',1.5)
load("ContinuousExperiment6thOrder.mat","SystemSizedG","SystemSizeFD")
loglog(Hq(7:12),SystemSizedG(7:12),'--','LineWidth',1.5)
hold on
loglog(Hq(7:12),SystemSizeFD(7:12),'-.','LineWidth',1.5)

%% Plot discrete errors from the Stoneley experiment
figure
% Plot Stoneley errors

load("Stoneley_EOC_4thOrder.mat")

loglog(Hq(7:12),Error(7:12),'-o','LineWidth',1.5,'Color','red')
hold on
loglog(Hq(7:12),3e4.*Hq(7:12).^4,'-.','LineWidth',1.5,'Color','red')

load("Stoneley_EOC_6thOrder.mat")
loglog(Hq(7:12),Error(7:12),'-x','LineWidth',1.5,'Color','blue')
hold on
loglog(Hq(7:12),4.5e5.*Hq(7:12).^5,'--','LineWidth',1.5,'Color','blue')
legend("4th order","$h^4$","6th order","$h^5$","Interpreter","latex","Location","southeast")
xlabel("$h$","Interpreter","latex")
ylabel("Discrete error")
fontsize(gcf,15,"points")
%% Plot the eigenvalues of the semidiscretizations with shifted node

figure
load("ShiftedSpectrumTest.mat","shift","Spectrum")
%delta=[1e-6 5e-6 1e-5 5e-5 5e-5 1e-4 5e-4 1e-3 5e-3 1e-2 5e-2];
semilogx(shift,Spectrum(1,:),'-x','LineWidth',1.5)
hold on
semilogx(shift,Spectrum(2,:),'-o','LineWidth',1.5)
semilogx(shift,Spectrum(3,:),'-^','LineWidth',1.5)
semilogx(shift,Spectrum(4,:),'-v','LineWidth',1.5)
xlabel("Shift parameter $\delta$","Interpreter","latex")
ylabel("$\max_i|e_i|$","Interpreter","latex")
legend("$h=1/10$","$h=1/15$","$h=1/20$","$h=1/30$","Interpreter","latex","location","best")
axis([0 0.1 1.6e4 2.2e4])
fontsize(gcf,15,"points")

%% Plot Gaussian pulse saved at six different time steps in U, U(1,:),U(2,:),U(5,:) and U(6,:) used in manuscript

load("GaussianTestData.mat")
GammaCoord=1/5;

u3=U(6,:);
surf(X,Y,reshape(sqrt(u3(1:nx*ny).^2+u3(nx*ny+1:2*nx*ny).^2),[ny,nx]))
hold on
yline(GammaCoord,'LineWidth',2)
    
trisurf(T_linear',full(P_linear(1,:)),full(P_linear(2,:)),full(sqrt(u3(2*nx*ny+p1_indices).^2+u3(2*nx*ny+NDoFs+p1_indices).^2)))                 
set(gca, 'color', 'none');
set(gca,'LooseInset',get(gca,'TightInset'));
           
axis off   
axis equal
           
shading interp
colormap jet
caxis([-5 155])
grid off
hold off
view(2)

%colorbar

%%
close all

load("UpperSurvey.mat")
NS=[length(nonzeros(US(2,:)));length(nonzeros(US(4,:)));length(nonzeros(US(6,:)));length(nonzeros(US(8,:)))]
figure
plot(US(2,6:6:NS(1)),US(1,6:6:NS(1)),'LineWidth',1.5)
hold on
plot(US(4,8:8:NS(2)),US(3,8:8:NS(2)),'LineWidth',1.5)
plot(US(6,12:12:NS(3)),US(5,12:12:NS(3)),'LineWidth',1.5)
plot(US(8,16:16:NS(4)),US(7,16:16:NS(4)),'LineWidth',1.5)
legend("h=1/60","h=1/80","h=1/120","h=1/160","location","northwest")
xlabel("time $t$","interpreter","latex")
ylabel("$|u(0.1,1.8,t)|$","Interpreter","latex")
ax=gca;
ax.FontSize=16;

load("LowerSurvey.mat")
NS=[length(nonzeros(US(2,:)));length(nonzeros(US(4,:)));length(nonzeros(US(6,:)));length(nonzeros(US(8,:)))]
figure
plot(US(2,6:6:NS(1)),US(1,6:6:NS(1)),'LineWidth',1.5)
hold on
plot(US(4,8:8:NS(2)),US(3,8:8:NS(2)),'LineWidth',1.5)
plot(US(6,12:12:NS(3)),US(5,12:12:NS(3)),'LineWidth',1.5)
plot(US(8,(NS(4)-1)/2:16:NS(4)),US(7,(NS(4)-1)/2:16:NS(4)),'LineWidth',1.5)

legend("h=1/60","h=1/80","h=1/120","h=1/160","location","northwest")
xlabel("time $t$","interpreter","latex")
ylabel("$|u(0.1,0.9,t)|$","Interpreter","latex")
ax=gca;
ax.FontSize=16;