%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION:  Generates projection operators from  the linear glue grids xa and xb via
%               an intermediary glue grid x_g which contains both of the given ones. Note
%               that these projection operators are created to be exact up to the order
%               given as well as being norm compatible. E.g. if M_a and M_b corresponds
%               to the quadrature matrices on the separate glue grids it holds that
%               M_a*Pb2a=(M_b*Pa2b)'.
%
%                       xa:      x-----x-----x-----x-----x-----x-----x-----x-----x
%
%                       xg:      x---x-x-x---x---x-x-x---x---x-x-x---x---x-x-x---x
%
%                       xb:      x---x---x---x---x---x---x---x---x---x---x---x---x
%
% INPUT : - xa,xb: vectors containing the two grids that we want to project
%                  between
%
%        - order: the order we want our projection operators to be exact up
%                 to
%
%
% OUTPUT: - PA2B, PB2A: projection operators between the two grids xa and
%                 xb which project through a common glue mesh
%
% AUTHOR:  Andreas Granath (andreas.granath@umu.se)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [PA2B,PB2A]=glue_projections(xa,xb,order)







 N=order+1;

xa=[xa(1:N:end),xa(end)]; xb=[xb(1:N:end),xb(end)];

Na=length(xa); Nb=length(xb);
ha=1/(Na-1); hb=1/(Nb-1);
TOL=100*eps;
x_g=sort(uniquetol([xa,xb],TOL));
%x_g(1)=[]; x_g(end)=[];


    n_g=length(x_g)-1;
   
 % Assemble "global" mass matrix for the glue x_g by looping over local
 % contributions

    for i=1:n_g
        x_i_g1=x_g(i); x_i_g2=x_g(i+1);
        h_gi=x_i_g2-x_i_g1;
        inds=(order+1)*(i-1)+1:(order+1)*i;
        M_g(inds,inds)=h_gi/2.*NewoneDimMassMatrix(@LagrangeRbf,order+1);
    end

% Initialize projection operators to the glue

Pa2g=zeros(N*n_g,N*(Na-1));
Pb2g=zeros(N*n_g,N*(Nb-1));

% Construct mass matrices for x_a and x_b needed to assemble the projection
% operators in the reverse direction using norm compatibility

M_a=ha/2*kron(eye(Na-1),NewoneDimMassMatrix(@LagrangeRbf,order+1));
M_b=hb/2*kron(eye(Nb-1),NewoneDimMassMatrix(@LagrangeRbf,order+1));

for k_l=1:Na-1

    % Determine nodes bounding the current cell

    xi1_l=xa(k_l); xi2_l=xa(k_l+1);

    % Determine corresponding internal nestled points in the glue

    k_g=find(all([x_g>=xi1_l-TOL;x_g< xi2_l+TOL]));

    x_g_lag=sym(AddLagrangeNodes(x_g(k_g),order));
    x_a_loc=sym(AddLagrangeNodes([xi1_l,xi2_l],order));

    % Assemble Vandermonde matrices on x_a and x_g

    for j=1:order
        V_a_kl=Vandermonde(x_a_loc,order);
        V_g_kl=Vandermonde(x_g_lag,order);
    end

    if k_l==1 || k_l==Na-1
        %V_g_kl
        %cond(V_g_kl)

    end

    % Construct local projection operator using pseudoinverse of the
    % Vandermonde matrix V_a_kl

    Pa2g_loc=V_g_kl/V_a_kl;

    % Add local contribution to global projection operator
    Pa2g(1+(k_g(1)-1)*N:N*(k_g(end)-1),1+(k_l-1)*N:N*k_l)=Pa2g_loc;
          
end

Pg2a=sym(M_a)\(sym(M_g)*Pa2g)';

Pa2g=double(Pa2g); Pg2a=double(Pg2a);
for k_l=1:Nb-1

    % Determine nodes bounding the current cell

    xi1_l=xb(k_l); xi2_l=xb(k_l+1);

    % Determine corresponding internal nestled points in the glue

    k_g=find(all([x_g< xi2_l+TOL ; x_g>=xi1_l-TOL]));

    x_g_lag=sym(AddLagrangeNodes(x_g(k_g),order));
    x_b_loc=sym(AddLagrangeNodes([xi1_l,xi2_l],order));

    % Assemble Vandermonde matrices on x_b and x_g

    V_b_kl=sym(Vandermonde(x_b_loc,order));
    V_g_kl=sym(Vandermonde(x_g_lag,order));

    % Construct local projection operator using pseudoinverse of the

    % Vandermonde matrix V_b_kl

    Pb2g_loc=V_g_kl/V_b_kl;


   % Add local contribution to global projection operator
   
    Pb2g(1+(k_g(1)-1)*N:N*(k_g(end)-1),1+(k_l-1)*N:N*k_l)=Pb2g_loc;          
end

Pg2b=sym(M_b)\(sym(M_g)*Pb2g)';

Pb2g=double(Pb2g); Pg2b=double(Pg2b);
% Construct operators via the glue

PA2B=Pg2b*Pa2g;
PB2A=Pg2a*Pb2g;



end


