

function [Pc2f,Pf2c]=FineToCoarseProj(xa,xg,order)


TOL=10*eps;

N=order+1;
Na=length(xa)-1; 
Ng=length(xg)-1;


Pc2f=zeros(N*Ng,N*Na);
Pf2c=zeros(N*Na,N*Ng);



for i=1:Na
    x_l=xa(i); x_r=xa(i+1);
    
    k_g=find([xg>=x_l-TOL & xg<=x_r+TOL]);
    
    % Construct operator from coarse element to fine element

    x_c=AddLagrangeNodes([x_l,x_r],order);
    x_f=AddLagrangeNodes(xg(k_g),order);
    Ng=length(k_g);

    V_c=Vandermonde(x_c,order);
    V_f=Vandermonde(x_f,order);
    Pc2f_loc=V_f/V_c;

    M_c_loc=(x_r-x_l)/2*NewoneDimMassMatrix(@LagrangeRbf,order+1);
    M_f_loc=zeros((Ng-1)*N,(Ng-1)*N);
    for j=1:Ng-1
        h_f=x_f(N*j)-x_f(1+N*(j-1));
        M_f_loc(1+N*(j-1):N*j,1+N*(j-1):N*j)=h_f/2*NewoneDimMassMatrix(@LagrangeRbf,order+1);  
    end

    Pf2c_loc=M_c_loc\((M_f_loc*Pc2f_loc)');

    diffMax=max(Pf2c_loc*x_f'-x_c');
    if diffMax>1e-10
        1
    end

    ind_c=1+(i-1)*N;
    ind_f=(k_g(1)-1)*N+1;

    Pf2c(ind_c:ind_c+N-1,ind_f:ind_f+N*(Ng-1)-1)=Pf2c_loc;
    Pc2f(ind_f:ind_f+N*(Ng-1)-1,ind_c:ind_c+N-1)=Pc2f_loc;

           
end
2

end