function [x,y]=bdry_fnc(bs,s)
 GammaCoord=1/5;
 switch nargin
     case 0
         x=4;
     case 1
         A=[0 1 2 3
            1 2 3 4
            1 1 1 1
            0 0 0 0];
         x=A(:,bs);
     case 2
         x=(s>=0 & s<1)+(2-s).*(s>=1 & s<2)-(3-s).*(s>=3 & s<=4);
         y=GammaCoord.*s.*(s>=0 & s<1)+GammaCoord.*(s>=1 & s<2)+...
          GammaCoord.*(3-s).*(s>=2 & s<3)+0.4.*(s-3).*(s-4).*exp(-10.*(s-3.5).^2).*(s>=3 & s<=4);
 end
end