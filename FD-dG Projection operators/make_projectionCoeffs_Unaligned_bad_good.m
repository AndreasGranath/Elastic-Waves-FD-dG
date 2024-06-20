% Limit for working is (Nb1,Nb2,Nh)=(4,24,20)

na=17; nb=na; order=3;

% na points on FD side
% nb nodes, so order*nb points
xa=[-1:2/(na-1):1];
xg=[-1:2/(nb-1):1];
xg=AddLagrangeNodes(xg,3);
% Set up parameters for projection operator
Nb1=4; Nb2=24; Nh=28;
No=0.5*(8*Nb1-Nh+4*(nb-na));

H=sym(SBP4(na,2/(na-1)));
M=1/((nb-1))*sym(kron(eye(nb-1),oneDimMassMatrix(@LagrangeRbf,order+2)));

% Initialize variables
X=sym('X',[Nb1*Nb2+Nh 1]);

% Calculate the number of columns
Nr=na-2*Nb1;

% Initialize projection operator
Pa2b=sym(zeros(4*(nb-1),na));
Pa2b(1:Nb2,1:Nb1)=rot90(reshape(X(1:Nb2*Nb1),[Nb2 Nb1]),2);
Pa2b(end-(Nb2-1):end,end-(Nb1-1):end)=reshape(X(1:Nb1*Nb2),[Nb2 Nb1]);

for i=1:Nr+1
    Pa2b(No+4*(i-1)+1:No+Nh+4*(i-1),Nb1+i)=reshape(X(Nb1*Nb2+1:end), [Nh,1]);
end

% Set up equations
x1=sym(ones(1,na))'; x2=sym(xa)'; x3=sym(xa).^2'; x4=sym(xa).^3';
xg1=sym(ones(1,4*(nb-1)))'; xg2=sym(xg)'; xg3=sym(xg).^2'; xg4=sym(xg).^3';

eqn1=Pa2b*x1-xg1; eqn2=Pa2b*x2-xg2; eqn3=Pa2b*x3-xg3; eqn4=Pa2b*x4-xg4;

Pb2a=H\((M*Pa2b)');

eqn5=Pb2a*xg1-x1; eqn6=Pb2a*xg2-x2; eqn7=Pb2a*xg3-x3; eqn8=Pb2a*xg4-x4;

EQN1=[eqn1;eqn2];
%EQN2=[eqn5;eqn6];

tic 
S=solve([EQN1;EQN2],X,'ReturnConditions',true);
toc 

vars=S.parameters;
S=struct2cell(S);



% Assemble into projection operator
%%
PA2B=sym(zeros(4*(na-1),na));
PA2B(1:Nb2,1:Nb1)=rot90(reshape(S(1:Nb1*Nb2),[Nb2 Nb1]),2);
PA2B(end-(Nb2-1):end,end-(Nb1-1):end)=reshape(S(1:Nb1*Nb2),[Nb2 Nb1]);

for i=1:Nr
    PA2B(No+4*(i-1)+1:No+Nh+4*(i-1),Nb1+i)=reshape(S(Nb1*Nb2+1:end-2), [Nh,1]);
end

PB2A=H\((M*PA2B)');




%PA2B=subs(Pa2b,vars,Y); Pb2a=subs(Pb2a,vars,Y);

Pf1=matlabFunction(PA2B);
Pf2=matlabFunction(PB2A);

InitGuess=zeros(length(vars),1);

f=fminsearch(@(x)l2error_bad_good(Pf1,Pf2,x,xa',xg'),InitGuess);

% Insert optimzed coefficients
PA2B=double(subs(PA2B,vars,f'));
PB2A=double(subs(PB2A,vars,f'));

bg_Block1=PA2B(1:Nb2,1:Nb1);
bg_Block2=PA2B(end-Nb2+1:end,end-Nb1+1:end);
Interior=nonzeros(PA2B(No+1:No+Nh,Nb1+1));
save("bad_good_Coeffs_21.mat","bg_Block1","bg_Block2","Interior","Nb1","Nb2","Nh","No","Nr");

% Plot pointwise errors
figure

subplot(2,2,1)
plot(xa,PB2A*xg.^0'-xa.^0')
title('Constant')
subplot(2,2,2)
plot(xa,PB2A*xg.^1'-xa.^1')
title('x^1')
subplot(2,2,3)
plot(xa,PB2A*xg.^2'-xa.^2')
title('x^2')
subplot(2,2,4)
plot(xa,PB2A*xg.^3'-xa.^3')
title('x^3')

figure
title('Point-wise errors for good operator PFD2dG')
subplot(2,2,1)
plot(xg,PA2B*xa.^0'-xg.^0','-r')
title('Constant')
subplot(2,2,2)
plot(xg,PA2B*xa.^1'-xg.^1','-r')
title('x^1')
subplot(2,2,3)
plot(xg,PA2B*xa.^2'-xg.^2','-r')
title('x^2')
subplot(2,2,4)
plot(xg,PA2B*xa.^3'-xg.^3','-r')
title('x^3')
