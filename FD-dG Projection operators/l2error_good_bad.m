
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION: calculates the discrete l^2 projection error to the leading
%              order as dictated by the input 'order' when projecting          
%              between xa and xb with a good operator Pf1 and between xb and
%              xa using a bad operator Pf2.
%              
%
% INPUT : - xa,xb: vectors containing the two grids that we want to project
%                  between
%
%        - order: the polynomial order of the dG method used
%
%        - Pf2,Pf2: projection operators between xa and xb as functions of  
%                   a set of free parameters
%
%         - x: a vector containing a choice of free parameters to use in
%             the projection operators
%
% OUTPUT: - err: a scalar, the l^2 error
%
% AUTHOR:  Andreas Granath (andreas.granath@umu.se)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err=l2error_good_bad(Pf1,Pf2,x,xa,xb,order)
    xcell=num2cell(x);

    switch order
        case 3
             err=sqrt(sum((Pf1(xcell{:})*xa.^3-xb.^3).^2)+sum((Pf2(xcell{:})*xb.^2-xa.^2).^2));
       case 4
             err=sqrt(sum((Pf1(xcell{:})*xa.^4-xb.^4).^2)+sum((Pf2(xcell{:})*xb.^3-xa.^3).^2));      
    end
end