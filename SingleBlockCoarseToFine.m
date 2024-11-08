
    
function [Pc2f,Pf2c]=SingleBlockCoarseToFine(xc,xf,order)
    N=order+1;



    nf=length(xf)-1;
    
    XC=sym(AddLagrangeNodes(xc,order));
    XF=sym(AddLagrangeNodes(xf,order));
    
    MC=sym(1/2*NewoneDimMassMatrix(@LagrangeRbf,order+1));
    MF=sym(zeros(N*nf,N*nf));
    for j=1:nf
        MF(1+N*(j-1):j*N,1+N*(j-1):j*N)=sym(1/2*NewoneDimMassMatrix(@LagrangeRbf,order+1));
    end
    
    
    pc2f=sym(Vandermonde(XF,order))/sym(Vandermonde(XC,order))
    
    pf2c=MC\((pc2f'*MF));

    Pc2f=double(pc2f); Pf2c=double(pf2c);
end
