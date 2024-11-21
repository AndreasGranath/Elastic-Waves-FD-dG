clear


order=3;
x=[0:0.05:1];
y=[0:1/81:1];  
g=uniquetol([x,y],10*eps);

[Px2g,Pg2x]=FineToCoarseProj(x,g,order);
[Py2g,Pg2y]=FineToCoarseProj(y,g,order);

X=AddLagrangeNodes(x,order);
Y=AddLagrangeNodes(y,order);

Px2y=Pg2y*Px2g;
Py2x=Pg2x*Py2g;

max(Px2y*X'-Y')
%%
clear
TOL=10*eps;
order=3;
N=300;
x=[0:1/N:1];
y=[0:1/N:1]; y(5)=y(5)-0.1/N;

g=uniquetol([x,y],10*eps);

[Px2g,Pg2x]=FineToCoarseProj(x,g,order);
[Py2g,Pg2y]=FineToCoarseProj(y,g,order);

X=AddLagrangeNodes(x,order);
Y=AddLagrangeNodes(y,order);

Px2y=Pg2y*Px2g;
Py2x=Pg2x*Py2g;

max(Px2y*X'-Y')

