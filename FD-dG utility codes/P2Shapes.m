function [S,dSdr,dSds]=P2Shapes(r,s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION: P^2 dG shape functions on the isogeometric reference element
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

    S=[1-3*r-3*s+2*r^2+2*s^2+4*r*s;
        2*r^2-r;
        2*s^2-s;
        4*r*s;
        4*s-4*r*s-4*s^2;
        4*r-4*r^2-4*r*s];
    dSdr=[-3+4*r+4*s;4*r-1;0;4*s;-4*s;4-8*r-4*s];
    dSds=[-3+4*r+4*s;0;4*s-1;4*r;4-4*r-8*s;-4*r];
end