%% Add "Data" folder to path before running this code

close all
%% Plot discrete errors for the trigonometric experiment in a continuous material
figure

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

%% Plot discrete errors from the Stoneley experiment
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
%% Plot the eigenvalues of the semidiscretizations with shifted node

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
load("SeismicSurvey2.mat")
NS=[length(nonzeros(US(2,:)));length(nonzeros(US(4,:)));length(nonzeros(US(6,:)));length(nonzeros(US(8,:)))]
figure
%plot(US(2,6:6:NS(1)),US(1,6:6:NS(1)),'LineWidth',1.5)
hold on
plot(US(4,8:8:NS(2)),US(3,8:8:NS(2)),'LineWidth',1.5)
plot(US(6,12:12:NS(3)),US(5,12:12:NS(3)),'LineWidth',1.5)
plot(US(8,16:16:NS(4)),US(7,16:16:NS(4)),'LineWidth',1.5)

US=[];
load("n200GaussianSurvey.mat")
NS=length(nonzeros(US(2,:)));
plot(US(2,20:20:NS(1)),US(1,20:20:NS(1)),'LineWidth',1.5)
legend("h=1/80","h=1/120","h=1/160","h=1/200","exact")
xlabel("time $t$","interpreter","latex")
ylabel("$|u(0.9,0.4,t)|$","Interpreter","latex")
ax=gca;
ax.FontSize=16;

%%
figure
plot(US(2,(NS(1)-1)/2:6:NS(1)),US(1,(NS(1)-1)/2:6:NS(1)),'LineWidth',1.5)
hold on
plot(US(4,(NS(2)-1)/2:8:NS(2)),US(3,(NS(2)-1)/2:8:NS(2)),'LineWidth',1.5)
plot(US(6,(NS(3)-1)/2:12:NS(3)),US(5,(NS(3)-1)/2:12:NS(3)),'LineWidth',1.5)
plot(US(8,(NS(4)-1)/2:16:NS(4)),US(7,(NS(4)-1)/2:16:NS(4)),'LineWidth',1.5)

legend("h=1/60","h=1/80","h=1/120","h=1/160","exact")
xlabel("time $t$","interpreter","latex")
ylabel("$|u(0.1,1.8,t)|$","Interpreter","latex")
ax=gca;
ax.FontSize=16;