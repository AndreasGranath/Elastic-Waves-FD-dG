%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION:  defines a map from [-1,1] to a given edge in the
%               isoparametric element with vertices at (0,0), (1,0) and     
%               (0,1), used to map Gaussian quadratures
%
% INPUT : - edgeNumber: edge in the isoparametric element
%
%
% OUTPUT: - r,s: functions mapping to the x and y dimensions respectively 
%                in the reference element 
%
% AUTHOR:  Andreas Granath (andreas.granath@umu.se)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r,s]=GaussianMapping(edgeNumber)

    switch edgeNumber
    
        case 1
            r=@(t) 0.5.*(-t+1); s=@(t) 0.5.*(1+t);
                        
        case 2
            r=@(t) 0.*t; s=@(t) 0.5.*(1-t);        
           
        case 3
            r=@(t) 0.5.*(t+1); s=@(t) 0.*t;
          
    end

   

end