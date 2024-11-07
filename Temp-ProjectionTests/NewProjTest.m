
addpath('..\FD-dG Projection operators\')
addpath('..\FD-dG utility codes\')

% Currently set to different number of nodes to show the error which
% arises, if they are set to be equal it behaves as expectedly

for q=1:1
na=314*q+1; % FD nodes

nb=172*q+1; % DG nodes, adds higher order nodes between these later

[H,~,~,~]=SBP4(na,1/(na-1)); %SBP quadrature

xa=[0:1/(na-1):1]; % FD DoFS
xb=[0:1/(nb-1):1]; % dG element boundaries

order=3;

% Construct operators from FD to piecewise polynomial, good and bad
% operators
[Pf2a_b,Pa2f_g,Pf2a_g,Pa2f_b]= make_projection_FD(na,order); 


% Set up glue grid

TOL=10*eps;
xg=sort(uniquetol([xa xb],TOL));


% Number of DoFS per element depending on order
N=order+1;


Na=length(xa)-1; Nb=length(xb)-1;
Ng=length(xg)-1;

% Initialize mass matrices
Mg=zeros(N*Ng,N*Ng);
Ma=zeros(N*Na,N*Na);
Mb=zeros(N*Nb,N*Nb);


% Construct projection operators from FD polynomial grid to glue grid and
% from dG grid with internal DoFS to glue grid

[Pa2g,Pg2a]=FineToCoarseProj(xa,xg,order);
[Pb2g,Pg2b]=FineToCoarseProj(xb,xg,order);


% Construct mass matrices
for j=1:Ng
    h_g=xg(j+1)-xg(j);
    Mg(1+(j-1)*N:N*j,1+(j-1)*N:N*j)=h_g/2*NewoneDimMassMatrix(@LagrangeRbf,order+1);
end

for j=1:Na
    h_a=xa(j+1)-xa(j);
    Ma(1+(j-1)*N:N*j,1+(j-1)*N:N*j)=h_a/2*NewoneDimMassMatrix(@LagrangeRbf,order+1);
end

for j=1:Nb
    h_b=xb(j+1)-xb(j);
    Mb(1+(j-1)*N:N*j,1+(j-1)*N:N*j)=h_b/2*NewoneDimMassMatrix(@LagrangeRbf,order+1);
end


% Use norm compatibility to construct operators from glue 
 %Pg2a=Ma\((Mg*Pa2g)');
 %Pg2b=Mb\((Mg*Pb2g)');

 % assemble the operators without using norm compatibility
%Pg2a=Pa2g\eye(N*Ng);
%Pg2b=Pb2g\eye(N*Ng);


% Construct full operators from FD DoFS to dG nodes

Pa2b=Pg2b*Pa2g*Pf2a_g;
Pb2a=Pa2f_b*Pg2a*Pb2g;

% Construct the piecewise polynomial grids
XB=AddLagrangeNodes(xb,order);
XA=AddLagrangeNodes(xa,order);
XG=AddLagrangeNodes(xg,order);
%% Figure
close all
plot(XB,Pa2b*xa.^0'-XB.^0','LineWidth',1.5)
hold on
plot(XB,Pa2b*xa'-XB','LineWidth',1.5)
plot(XB,Pa2b*xa.^2'-XB.^2','LineWidth',1.5)
plot(XB,Pa2b*xa.^3'-XB.^3','LineWidth',1.5)
legend('Constant','Linear','Quadratic','Cubic')
title('a to b')

figure
plot(xa,Pb2a*XB.^0'-xa.^0','LineWidth',1.5)
hold on
plot(xa,Pb2a*XB'-xa','LineWidth',1.5)

plot(xa,Pb2a*XB.^2'-xa.^2','LineWidth',1.5)
plot(xa,Pb2a*XB.^3'-xa.^3','LineWidth',1.5)
legend('Constant','Linear','Quadratic','Cubic')
title('b to a')

%% Display errors
disp(['Max linear error from a to glue ' num2str([max(Pa2g*XA'-XG')])])
disp(['Max linear error from b to glue ' num2str([max(Pb2g*XB'-XG')])])
disp(['Max constant error from a to b ' num2str([max(Pa2b*xa.^0'-XB.^0')])])
disp(['Max linear error from a to b ' num2str([max(Pa2b*xa'-XB')])])
disp(['Max quadratic error from a to b ' num2str([max(Pa2b*xa.^2'-XB.^2')])])
disp(['Max cubic error from a to b' num2str([max(Pa2b*xa.^3'-XB.^3')])])
disp(['Max constant error from b to a ' num2str([max(Pb2a*XB.^0'-xa.^0')])])
disp(['Max linear error from b to a ' num2str([max(Pb2a*XB'-xa')])])
disp(['Max quadratic error from b to a ' num2str([max(Pb2a*XB.^2'-xa.^2')])])
disp(['Max cubic error from b to a' num2str([max(Pb2a*XB.^3'-xa.^3')])])



       
max_const_err_a2b(q)=[max(Pa2b*xa.^0'-XB.^0')];
max_const_err_bra(q)=[max(Pb2a*XB.^0'-xa.^0')];
end

%% ONE CELL TEST


