function V=Vandermonde(v,order)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION: This code generates a Vandermonde matrix of a specified  
%               order for the function, used to construct projection
%               operators
%
% INPUT :- v: vector containing a grid evaluation of some function
%
%        - order: scalar defining the order of the Vandermonde matrix
%
%
% OUTPUT: - V: a Vandermonde matrix corresponding to the grid evaluation v
%
% AUTHOR: Andreas Granath (andreas.granath@umu.se)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    V=zeros(length(v),order+1);
    for i=1:order+1
        V(:,i)=v.^(i-1);
    end
end
