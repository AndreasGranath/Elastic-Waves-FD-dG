%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION: Constructs coefficients to use in the FD-FD glue projection
%              operator pair (P_{u2v}^g,P_{v2u}^b) given an order of
%              accuracy that these should follow. Note that the
%              optimization method used and the order in which the accuracy
%              and spectrum is optimized must be changed manually!
%              
%              Saves down the coefficients in separate files.
%
% INPUT : - N: number of FD nodes uniformily distributed on [-1,1]
%
%         - order: order of accuracy for the operators
%
% OUTPUT: - Pu2v: good projection operator from FD grid to piecewise
%                 discontinuous grid x_v having x_u as element boundaries

%         - Pv2u: "bad" projection operator from piecewise discontinuous
%                  grid x_v having x_u as element boundaries to x_u
%
% AUTHOR:  Andreas Granath (andreas.granath@umu.se)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [PA2B,PB2A]=new_make_projectionCoeffs_bad_good(na,order)


% Assembl FD grid on [-1,1]and the associated piecewise discontinuous grid
% xg
xa=[-1:2/(na-1):1];

xg=AddLagrangeNodes(xa,order);

% Set up parameters for projection operator (NOTE THAT THIS IS THE MINIMAL
% WORKIN CASE, DO NOT ALTER THESE TOO MUCH!)

if order==3
    Nb1=4; Nb2=24; Nh=20;
elseif order==4
    Nb1=6; Nb2=48; Nh=36;
end


No=0.5*(2*(order+1)*Nb1-Nh);
% Set up SBP-FD quadrature for norm compatibility
if order==3
    H=sym(SBP4(na,2/(na-1)));
elseif order==4
    H=sym(SBP4(na,2/(na-1)));
end

% Assemble dG mass matrix
M=1/((na-1))*sym(kron(eye(na-1),oneDimMassMatrix(@LagrangeRbf,order+2)));

% Initialize symbolic  variables
X=sym('X',[Nb1*Nb2+Nh 1]);

% Calculate the number of columns in the projection operator
Nr=na-2*Nb1;

% Initialize projection operator

Pa2b=sym(zeros((order+1)*(na-1),na));

% Load closure stencils
Pa2b(1:Nb2,1:Nb1)=rot90(reshape(X(1:Nb1*Nb2),[Nb2 Nb1]),2);
Pa2b(end-(Nb2-1):end,end-(Nb1-1):end)=reshape(X(1:Nb1*Nb2),[Nb2 Nb1]);

% Load repeated stencil into the interior of the projection operator
for i=1:Nr
    Pa2b(No+(order+1)*(i-1)+1:No+Nh+(order+1)*(i-1),Nb1+i)=reshape(X(Nb1*Nb2+1:end), [Nh,1]);
end

% Set up equations from accuracy conditions



if order==3

    x1=sym(ones(1,na))'; x2=sym(xa)'; x3=sym(xa).^2'; x4=sym(xa).^3';
    xg1=sym(ones(1,4*(na-1)))'; xg2=sym(xg)'; xg3=sym(xg).^2'; xg4=sym(xg).^3';
    
    eqn1=Pa2b*x1-xg1; eqn2=Pa2b*x2-xg2; eqn3=Pa2b*x3-xg3; eqn4=Pa2b*x4-xg4;
    
    % Use norm compatibility to set up the equations for Pb2a
    
    Pb2a=H\((M*Pa2b)');
    
    eqn5=Pb2a*xg1-x1; eqn6=Pb2a*xg2-x2; eqn7=Pb2a*xg3-x3; eqn8=Pb2a*xg4-x4;

    EQN1=[eqn1;eqn2;eqn3(Nb2+1:end-Nb2);eqn4(Nb2+1:end-Nb2)];
    EQN2=[eqn5;eqn6;eqn7;eqn8(Nb1+1:end-Nb1)];

else % 4th order

    x1=sym(ones(1,na))'; x2=sym(xa)'; x3=sym(xa).^2'; x4=sym(xa).^3'; x5=sym(xa).^4';
    xg1=sym(ones(1,(order+1)*(na-1)))'; xg2=sym(xg)'; xg3=sym(xg).^2'; xg4=sym(xg).^3'; xg5=sym(xg).^4';

    eqn1=Pa2b*x1-xg1; eqn2=Pa2b*x2-xg2; eqn3=Pa2b*x3-xg3; eqn4=Pa2b*x4-xg4;
    eqn5=Pa2b*x5-xg5;
    
    % Use norm compatibility to set up the equations for Pb2a
    
    Pb2a=H\((M*Pa2b)');
    
    eqn6=Pb2a*xg1-x1; eqn7=Pb2a*xg2-x2; eqn8=Pb2a*xg3-x3; eqn9=Pb2a*xg4-x4;
    eqn10=Pb2a*xg5-x5;

    EQN1=[eqn1;eqn2;eqn3;eqn4(Nb2+1:end-Nb2);eqn5(Nb2+1:end-Nb2)];
    EQN2=[eqn6;eqn7;eqn8;eqn9;eqn10(Nb1+1:end-Nb1)];

end


S=solve([EQN1;EQN2],X,'ReturnConditions',true);
vars=S.parameters;
S=struct2cell(S);



% Assemble into projection operator PA2B
%%
PA2B=sym(zeros(4*(na-1),na));
PA2B(1:Nb2,1:Nb1)=rot90(reshape(S(1:Nb1*Nb2),[Nb2 Nb1]),2);
PA2B(end-(Nb2-1):end,end-(Nb1-1):end)=reshape(S(1:Nb1*Nb2),[Nb2 Nb1]);

for i=1:Nr
    PA2B(No+4*(i-1)+1:No+Nh+4*(i-1),Nb1+i)=reshape(S(Nb1*Nb2+1:end-2), [Nh,1]);
end

% Obtain PV2U by norm compatibility

PB2A=H\((M*PA2B)');


Pf1=matlabFunction(PA2B);
Pf2=matlabFunction(PB2A);

% Start optimization over the free parameters, initial guess zeros
InitGuess=zeros(length(vars),1);

% Currently set to use "fminunc" for both minimizations, change these
% manually
%options=optimoptions('fminsearch','MaxFunctionEvaluations',1e6,'MaxIterations',1e5,'OptimalityTolerance',1e-7);
    options=optimset('MaxFunEvals',1e6,'MaxIter',1e5,'TolX',1e-7)
f=fminsearch(@(x)Eigen_Value_opt(Pf1,Pf2,x'),InitGuess,options);

% Immediately optimize over l^2-error

g=fminsearch(@(x)l2error_bad_good(Pf1,Pf2,x,xa',xg',order),f,options);

% Insert optimzed coefficients
PA2B=double(subs(PA2B,vars,g'));
PB2A=double(subs(PB2A,vars,g'));

% Insert optimzed coefficients
PA2B=double(subs(PA2B,vars,f'));
PB2A=double(subs(PB2A,vars,f'));

Block1=PA2B(1:Nb2,1:Nb1);
Block2=PA2B(end-Nb2+1:end,end-Nb1+1:end);
Interior=nonzeros(PA2B(No+1:No+Nh,Nb1+1));
save("bad_good_Coeffs_4th_order_AccuracyEig_Flattening_Fminsearch.mat","Block1","Block2","Interior","Nb1","Nb2","Nh","No")

% Plot pointwise errors


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
