%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION:  calculate quantity used to minimize spectrum of the FD-dG
%               projection operators Pf1,Pf2 for a given set of free
%               parameters x. This quantity is either a 'flattening norm'
%               or the Fröbenius norm.

% INPUT :   - Pf1,Pf2: FD-dG and dG-FD projection operators in matrix form
%                      but as functions of a set of free parameters
%
%           - x: a vector containing a choice of free parameters for the
%                projection oeprators
%
%
% OUTPUT: - err: a scalar quantity relating to the spectrum of the
%                operators
%
% AUTHOR:  Andreas Granath (andreas.granath@umu.se)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function err=Eigen_Value_opt(Pf1,Pf2,x)

    % Convert numerical choice of free parameters
    xcell=num2cell(x);

    % Calculate product of operators, this is what we want to approximate
    % an identity matrix

    P21=full(Pf2(xcell{:})*Pf1(xcell{:}));
    P12=full(Pf1(xcell{:})*Pf2(xcell{:}));

    % For flattening norm, calculate and sort eigenvalues of the products
    e1=sort(real(eig(P12))); e2=sort(real(eig(P21)));
   
    err=sqrt(norm(diff(e1))^2+norm(diff(e2))^2);

    % For Fröbenius norm
  %  err=sqrt(norm(P12,"fro")^2+norm(P21,"fro")^2);
    
end