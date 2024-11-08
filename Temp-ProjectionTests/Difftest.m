clear

ndG=20; 

xdG=AddLagrangeNodes([0:1/(ndG-1):1],3);
xFD=AddLagrangeNodes([0:1/(3*(ndG-1)):1],3);

[Pf2g_old,Pg2f_old]=glue_projections(xFD,xdG,3);

[Pf2g_new,Pg2f_new]=NewGlueProjections(ndG,3,1);

%%
x_c=[0,0.6];
x_f=[0,1/10,4/10,5/10,0.6];
X_f=AddLagrangeNodes(x_f,order);
X_c=AddLagrangeNodes(x_c,order);
[Pc2f,Pf2c]=SingleBlockCoarseToFine(x_c,x_f,order)
%%
Pf2g_diff=Pf2g_old-Pf2g_new;
Pg2f_diff=Pg2f_old-Pg2f_new;

disp(['Maximal difference f2g is ' num2str(max(max(Pf2g_diff)))])
disp(['Maximal difference g2f is ' num2str(max(max(Pg2f_diff)))])
%%
close all

% subplot(1,2,1)
% plot(Pf2g_old*xFD'-xdG')
% title('f2g old')
% subplot(1,2,2)
% plot(Pg2f_old*xdG'-xFD')
% title('g2f old')




figure
subplot(1,2,1)
plot(Pf2g_new*xFD'-xdG')
title('f2g new')
subplot(1,2,2)
plot(Pg2f_new*xdG'-xFD')
title('g2f new')