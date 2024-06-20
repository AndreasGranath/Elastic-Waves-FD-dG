function [S,dSdr,dSds]=P1Shapes(r,s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION: P^1 dG shape functions on the isogeometric reference element
%             with vertices at (0,0), (1,0) and (0,1)
%
% INPUT : - r,s: coordinates in the reference element where we want to 
%                evaluate the shape functions
%
% OUTPUT: - S: A vector containing the values of the shape functions at the
%              point (r,s).
%
%        -dSds,dSdr: vectors containing the derivatives of the shape      
%                    functions w.r.t r and s at the point (r,s).
%
% AUTHOR:  Andreas Granath (andreas.granath@umu.se)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S=[1-r-s; r; s];
    dSdr=[-1;1;0];
    dSds=[-1;0;1];
end