% Limit for working is (Nb1,Nb2,Nh)=(4,24,20)

function [PA2B,PB2A,error_acc,error_spec,eigen_quant_acc,eigen_quant_spec]=make_projectionCoeffs_good_bad_Higher_Order(na,order)


xa=[-1:2/(na-1):1];
xg=AddLagrangeNodes(xa,order);
% Set up parameters for projection operator
Nb1=6; Nb2=48; Nh=36;
% Nb1 tested: 6,8
%Nh tested:
% 38,46 with Nb2=48
% -- Nb2 tested --
% Nb2=53; with Nb1=8

No=0.5*(2*(order+1)*Nb1-Nh);
[D,HI]=diagonal_sbp(order+2,na-1);
H=2/(na-1)*sym(inv(HI));

M=1/((na-1))*sym(kron(eye(na-1),NewoneDimMassMatrix(@LagrangeRbf,order+1)));

% Initialize variables
X=sym('X',[Nb1*Nb2+Nh 1]);

% Calculate the number of columns
Nr=na-2*Nb1;

% Initialize projection operator
Pa2b=sym(zeros((order+1)*(na-1),na));
Pa2b(1:Nb2,1:Nb1)=rot90(reshape(X(1:Nb1*Nb2),[Nb2 Nb1]),2);
Pa2b(end-(Nb2-1):end,end-(Nb1-1):end)=reshape(X(1:Nb1*Nb2),[Nb2 Nb1]);

for i=1:Nr
    Pa2b(No+(order+1)*(i-1)+1:No+Nh+(order+1)*(i-1),Nb1+i)=reshape(X(Nb1*Nb2+1:end), [Nh,1]);
end

% Set up equations
x1=sym(ones(1,na))'; x2=sym(xa)'; x3=sym(xa).^2'; x4=sym(xa).^3';
x5=sym(xa).^4'; 

xg1=sym(ones(1,(order+1)*(na-1)))'; xg2=sym(xg)'; xg3=sym(xg).^2'; xg4=sym(xg).^3';
xg5=sym(xg).^4'; 

eqn1=Pa2b*x1-xg1; eqn2=Pa2b*x2-xg2; eqn3=Pa2b*x3-xg3; eqn4=Pa2b*x4-xg4;
eqn9=Pa2b*x5-xg5; 

Pb2a=H\((M*Pa2b)');

eqn5=Pb2a*xg1-x1; eqn6=Pb2a*xg2-x2; eqn7=Pb2a*xg3-x3; eqn8=Pb2a*xg4-x4;
eqn11=Pb2a*xg5-x5;

EQN1=[eqn1;eqn2;eqn3;eqn4;eqn9(Nb2+1:end-Nb2)];
EQN2=[eqn5;eqn6;eqn7;eqn8(Nb1+1:end-Nb1);eqn11(Nb1+1:end-Nb1)];

tic
S=solve([EQN1;EQN2],X,'ReturnConditions',true);
vars=S.parameters;
S=struct2cell(S);
toc



% Assemble into projection operator
%%
PA2Borg=sym(zeros((order+1)*(na-1),na));
PA2Borg(1:Nb2,1:Nb1)=rot90(reshape(S(1:Nb1*Nb2),[Nb2 Nb1]),2);
PA2Borg(end-(Nb2-1):end,end-(Nb1-1):end)=reshape(S(1:Nb1*Nb2),[Nb2 Nb1]);

for i=1:Nr
    PA2Borg(No+(order+1)*(i-1)+1:No+Nh+(order+1)*(i-1),Nb1+i)=reshape(S(Nb1*Nb2+1:end-2), [Nh,1]);
end

% Note PA2B=Pw2uB, PB2A=Pu2wG; In paper, want Pu2wG=H\((M*Pw2uB))

PB2Aorg=H\((M*PA2Borg)');



%PA2B=subs(Pa2b,vars,Y); Pb2a=subs(Pb2a,vars,Y);

Pf1=matlabFunction(PA2Borg);
Pf2=matlabFunction(PB2Aorg);

ll=length(vars);
InitGuess=zeros(ll,1);

options=optimset('MaxFunEvals',1e6,'MaxIter',1e5,'TolX',1e-7)
f=fminsearch(@(x)l2error_good_bad(Pf1,Pf2,x,xa',xg',order),InitGuess,options);

PA2B=double(subs(PA2Borg,vars,f'));
PB2A=double(subs(PB2Aorg,vars,f'));

Block1=PA2B(1:Nb2,1:Nb1);
Block2=PA2B(end-Nb2+1:end,end-Nb1+1:end);
Interior=PA2B(No+1:No+Nh,Nb1+1);

%eigen_quant_acc=norm(PA2B*PB2A,"fro")+norm(PB2A*PA2B,"fro");


g=fminsearch(@(x)Eigen_Value_opt(Pf1,Pf2,x'),f,options);

% Insert optimzed coefficients
PA2B=double(subs(PA2Borg,vars,g'));
PB2A=double(subs(PB2Aorg,vars,g'));
error_spec=norm(PA2B*xa.^4'-xg.^4')+norm(PB2A*xg'.^3-xa'.^3);

Block1=PA2B(1:Nb2,1:Nb1);
Block2=PA2B(end-Nb2+1:end,end-Nb1+1:end);
Interior=PA2B(No+1:No+Nh,Nb1+1);

save("new_good_bad_Coeffs_4th_order_AccuracyEigfminsearch.mat","Block1","Block2","Interior","Nb1","Nb2","Nh","No")

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
plot(xa,PB2A*xg.^4'-xa.^4')
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
plot(xg,PA2B*xa.^4'-xg.^4','-r')
title('x^3')
end
