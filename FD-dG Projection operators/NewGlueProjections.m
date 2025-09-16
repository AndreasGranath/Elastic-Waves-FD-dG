function [PFD2dG,PdG2FD]=NewGlueProjections(ndG,order,flag)
    
      if flag==0
         nFD=ndG;
         scaling=1;
      elseif flag==1
         nFD=order*(ndG-1)+1;
         scaling=order;
      end

    x1=sym([0,1]);
    x2=sym([0:1/scaling:1]);
    
    X1=sym(AddLagrangeNodes(x1,order));
    X2=sym(AddLagrangeNodes(x2,order));
    
    M1=sym(1/2*NewoneDimMassMatrix(@LagrangeRbf,order+1));
    M2=sym(1/(2*scaling)*kron(eye(scaling),NewoneDimMassMatrix(@LagrangeRbf,order+1)));
    
    p12loc=sym(Vandermonde(X2,order))/sym(Vandermonde(X1,order));
    
    p21loc=M1\((p12loc'*M2));

    P12loc=double(p12loc); P21loc=double(p21loc);
    hdG=1/(ndG-1); hFD=1/(nFD-1);
    XdG=AddLagrangeNodes([0:1/(ndG-1):1],order);
    XFD=AddLagrangeNodes([0:1/(nFD-1):1],order);
    
    PdG2FD=kron(eye(ndG-1),P12loc);
    PFD2dG=kron(eye(ndG-1),P21loc);

end