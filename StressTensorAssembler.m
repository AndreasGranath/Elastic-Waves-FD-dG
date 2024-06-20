function [s11,s12,s22,s11_tt,s12_tt,s22_tt,f,f_tt]=StressTensorAssembler(u,v,mu,lambda,rho)
    syms x y t
    % Assemble normal

    u=sym(u);
    v=sym(v);
 
    u_x=diff(u,x); u_y=diff(u,y); v_x=diff(v,x); v_y=diff(v,y);

    % Use parametrized normal N(t)=[-h'(t),1]/sqrt(1+h'(t)^2) for top surface
    
    
    s11=(2*mu+lambda).*u_x+lambda.*v_y;
    s22=(2*mu+lambda).*v_y+lambda.*u_x;
    s12=mu.*(u_y+v_x);
    s11_tt=diff(s11,t,2); s12_tt=diff(s12,t,2); s22_tt=diff(s22,t,2);
    
   % S=[s11,s12;s12,s22];
   % S=matlabFunction(S);
   
   s11=matlabFunction(s11); s22=matlabFunction(s22); s12=matlabFunction(s12);
   s11_tt=matlabFunction(s11_tt); s22_tt=matlabFunction(s22_tt); s12_tt=matlabFunction(s12_tt);

    % Assemble load
    f1=rho*diff(u,t,2)-mu*(diff(u,x,2)+diff(u,y,2))-(mu+lambda)*(diff(u,x,2)+diff(v_x,y,1));
    f2=rho*diff(v,t,2)-mu*(diff(v,x,2)+diff(v,y,2))-(mu+lambda)*(diff(u_x,y,1)+diff(v,y,2));
    f=matlabFunction([f1;f2]);
    f_tt=matlabFunction(diff([f1;f2],t,2));
 


   



end