function [S,dSdx,dSdy,detJ]=Isopmap(xV,yV,r,s,shapefcn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION: evaluates the set of dG basis functions and collects them
%              in a vector S. Also calculates the x- and y- derivatives and
%              the determinant of the Jacobian corresponding to the
%              reference mapping of the current element
%
% INPUT : - xV,yV: vectors with x and y coordinates in physical space for
%                  the current element
%         - r,s:   points in reference element
%         - shapefcn : anonymous function corresponding to basis functions
%
%
% OUTPUT: -
%
% AUTHOR:  Andreas Granath (andreas.granath@umu.se)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Evaluate reference basis functions at the point (r,s) 
    [S,dSdr,dSds]=shapefcn(r,s);

    % Determine coefficients in map from reference to physical element
    rc=[xV',yV',ones(3,1)]\[0;1;0];
    sc=[xV',yV',ones(3,1)]\[0;0;1];

    % Calculate the derivative of basis functions in physical space

    dSdx=dSdr*rc(1)+dSds*sc(1);
    dSdy=dSdr*rc(2)+dSds*sc(2);

    % Calculate the Jacobian
    detJ=1/(rc(1)*sc(2)-rc(2)*sc(1));
    
end