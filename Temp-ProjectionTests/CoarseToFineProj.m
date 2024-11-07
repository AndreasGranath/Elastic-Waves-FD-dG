

function [Pa2g,Pg2a]=CoarseToFineProj(xa,xg,order)


TOL=10*eps;

N=order+1;
Na=length(xa)-1; 
Ng=length(xg)-1;


Pa2g=zeros(N*Ng,N*Na);
Pg2a=zeros(N*Na,N*Ng);



for i=1:Na
    x_l=xa(i); x_r=xa(i+1);
    
    k_g=find([xg>=x_l-TOL & xg<=x_r+TOL]);
    
    % Construct operator from coarse element to fine element

    x_c=AddLagrangeNodes([x_l,x_r],order);
    x_f=AddLagrangeNodes(xg(k_g),order);
    Ng=length(k_g);

    Pa2g_loc=zeros((Ng-1)*N,N);
    Mg_loc=zeros(N*(Ng-1),N*(Ng-1));
    Ma=(x_l-x_r)/2*NewoneDimMassMatrix(@LagrangeRbf,order+1);

    for j=1:Ng-1
        x_f_loc=x_f(1+N*(j-1):N*j);       
        Pa2g_loc(1+N*(j-1):N*j,1:N)=Vandermonde(x_f_loc,order)/Vandermonde(x_c,order);   
        Mg_loc(1+N*(j-1):N*j,1+N*(j-1):N*j)=(xg(k_g(j))-xg(k_g(j+1)))/2*NewoneDimMassMatrix(@LagrangeRbf,order+1);
    end

    ind_N=(i-1)*N+1;
    ind_g=(k_g(1)-1)*N+1;
    Pa2g(ind_g:ind_g+N*(Ng-1)-1,ind_N:ind_N+N-1)=Pa2g_loc;
    Pg2a(ind_N:ind_N+N-1,ind_g:ind_g+N*(Ng-1)-1)=Ma\((Mg_loc*Pa2g_loc)');


            
end

end