

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION: assembles the pair (P_{u2v}^b,P_{v2u}^g) of norm compatible
%              projection operators between an FD grid on N nodes and the   
%              resulting piecewise discontinuous grid having the N FD nodes
%              as element boundaries.The operators are exact up to an order 
%              as specified by 'order'.  
%              
%            x1 - - - x2 - - - x3 - - - x4 . . . xN-1 - - - xN      (x_u)
%                       |                      ^       
%                       | P_{u2v}              |  P_{v2u}       
%                       v                      |
%      x1 - - - x2 x2 - - - x3 x3 - - - x4 . . . xN-1 XN-1 - - - xN (x_v)
%              
%
% INPUT : -N : a scalar denoting how many equidistribued notes we have in
%              the FD grid
%
%        - order: the desired order of accuracy for the projection
%                 operators
%
% OUTPUT: - Pu2v_g/Pu2v_b: Matrices being "good" and "bad" projection operator between the FD
%                          grid x_u and the piecewise discontinuous grid x_v
%
%         - Pv2u_g/Pv2u_b: Matrices being "good" and "bad" projection operators
%                          between the piecewise discontinuous grid x_v and FD grid x_u
%
% AUTHOR:  Andreas Granath (andreas.granath@umu.se)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Pu2v_b,Pv2u_g,Pu2v_g,Pv2u_b]= make_projection_FD(N,order)

    % Assemble 1D mass matrix along the grid x_v

    M=0.5/((N-1))*kron(eye(N-1),NewoneDimMassMatrix(@LagrangeRbf,order+1));

    % Load pre-assembled coefficients for the projection operators
    % depending on the order

    if order==3
      load('bad_good_Coeffs_3rd_order_AccuracyEigfminsearch.mat');

      % Assemble inverse SBP norm

       [~,HI]= diagonal_sbp(4,N-1);

    elseif order==4
      load('good_bad_Coeffs_4th_order_EigAccuracy.mat');

      % Assemble inverse SBP norm
       [~,HI]= diagonal_sbp(6,N-1);
    end

    % Obtain SBP norm 
    H=1/(N-1)*(HI\eye(size(HI)));

    % Initialize projection between x_u and x_v
    Pu2v_b=zeros((order+1)*(N-1),N);
 
    % Allocate boundary stencils
    Pu2v_b(1:Nb2,1:Nb1)=Block1;
    Pu2v_b(end-Nb2+1:end,end-Nb1+1:end)=Block2;

    % Allocate repeated interior stencil

     Nr=N-2*Nb1;
    for i=1:Nr
        Pu2v_b(No+(order+1)*(i-1)+1:No+Nh+(order+1)*(i-1),Nb1+i)=Interior;
    end

    % Construct projection between x_v and x_u using norm-compatibility
    Pv2u_g=H\((M*Pu2v_b)');
    if order==3
      load('good_bad_Coeffs_3rd_order_AccuracyEigfminsearch.mat');
 
    elseif order==4
      load('bad_good_Coeffs_4th_order_EigAccuracy.mat');
  
    end

    Pu2v_g=zeros((order+1)*(N-1),N);
    Pu2v_g(1:Nb2,1:Nb1)=Block1;
    Pu2v_g(end-Nb2+1:end,end-Nb1+1:end)=Block2;
    Nr=N-2*Nb1;
    for i=1:Nr
        Pu2v_g(No+(order+1)*(i-1)+1:No+Nh+(order+1)*(i-1),Nb1+i)=Interior;
    end

    Pv2u_b=H\((M*Pu2v_g)');

   

end