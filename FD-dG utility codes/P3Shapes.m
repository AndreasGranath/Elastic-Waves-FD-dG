function [S,dSdr,dSds]=P3Shapes(r,s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION: P^3 dG shape functions on the isogeometric reference element
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

    S=[1-5.5*r-5.5*s+9*r^2+18*r*s+9*s^2-4.5*r^3-13.5*r^2*s-13.5*r*s^2-4.5*s^3;
       r-4.5*r^2+4.5*r^3;
       s-4.5*s^2+4.5*s^3;
       -4.5*r*s+13.5*r^2*s;
       -4.5*r*s+13.5*r*s^2;
       -4.5*s+4.5*r*s+18*s^2-13.5*s^2*r-13.5*s^3;
       9*s-22.5*r*s-22.5*s^2+13.5*r^2*s+13.5*s^3+27*s^2*r;
       9*r-22.5*r^2-22.5*r*s+13.5*r^3+13.5*s^2*r+27*r^2*s;
       -4.5*r+18*r^2+4.5*r*s-13.5*r^3-13.5*r^2*s;
       27*r*s-27*r^2*s-27*s^2*r];

    dSdr=[-5.5+18*r+18*s-13.5*r^2-27*r*s-13.5*s^2;
          1-9*r+13.5*r^2;
          0;
          -4.5*s+27*r*s;
          -4.5*s+13.5*s^2;
          4.5*s-13.5*s^2;  
          -22.5*s+27*r*s+27*s^2;
          9-45*r-22.5*s+40.5*r^2+13.5*s^2+54*r*s;
          -4.5+36*r+4.5*s-40.5*r^2-27*r*s;
          27*s-54*r*s-27*s^2];

    dSds=[-5.5+18*r+18*s-13.5*r^2-27*r*s-13.5*s^2;
          0;
          1-9*s+13.5*s^2;
          -4.5*r+13.5*r^2;
          -4.5*r+27*r*s;
          -4.5+4.5*r+36*s-27*r*s-40.5*s^2;
          9-22.5*r-45*s+13.5*r^2+40.5*s^2+54*r*s;
          -22.5*r+27*s*r+27*r^2;
          4.5*r-13.5*r^2;
          27*r-27*r^2-54*r*s];

    
end