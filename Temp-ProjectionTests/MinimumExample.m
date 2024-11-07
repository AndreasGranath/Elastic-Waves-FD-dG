clear
order=3; 
Na=2; Nb=4; N=order+1;

TOL=100*eps;

xa=AddLagrangeNodes([0:1/(Na-1):1],order)';
xb=AddLagrangeNodes([0:1/(Nb-1):1],order)';

% Construct vandermonde matrix
Va=Vandermonde(xa,order);
Vb=Vandermonde(xb,order);

Ma=0.5/(Na-1)*kron(eye(Na-1),NewoneDimMassMatrix(@LagrangeRbf,order+1));
Mb=0.5/(Nb-1)*kron(eye(Nb-1),NewoneDimMassMatrix(@LagrangeRbf,order+1));

% Here norm compatibility works
Pa2b=Vb/Va; Pb2a=Ma\(Pa2b'*Mb);




%% Close all
close all

figure
plot(Pa2b*ones((Na-1)*N,1)-ones((Nb-1)*N,1),'--')
hold on
plot(Pa2b*xa-xb)
plot(Pa2b*xa.^2-xb.^2)
plot(Pa2b*xa.^3-xb.^3)
figure
plot(Pb2a*ones((Nb-1)*N,1)-ones((Na-1)*N,1),'-.')
hold on
plot(Pb2a*xb-xa)
plot(Pb2a*xb.^2-xa.^2)


