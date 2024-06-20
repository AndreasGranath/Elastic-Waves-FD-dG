% Assembles the Stoneley interface solution for the situation with two
% blocks of width x_l meeting along the x-axis. We follow the notation from
% the paper:A HIGH ORDER FINITE DIFFERENCE METHOD FOR THE
%           ELASTIC WAVE EQUATION IN BOUNDED DOMAINS WITH
%            NONCONFORMING INTERFACES

function [u,v,u1,u2,v1,v2,c,v_tt]=TwoBlockStoneley(c_0,mu1,mu2,lambda1,lambda2,rho1,rho2,y_l,k)
   
   
   % Solve for the interface wave speed c_s given the material parameters
    cs1=sqrt(mu1/rho1); cs2=sqrt(mu2/rho2);

    cp1=sqrt((lambda1+2*mu1)/rho1); cp2=sqrt((lambda2+2*mu2)/rho2);

    s2=@(c)k*sqrt(1-(c/cs2)^2); s1=@(c) k*sqrt(1-(c/cs1)^2);
   r2=@(c) k*sqrt(1-(c/cp2)^2); r1=@(c) k*sqrt(1-(c/cp1)^2);


%     D=@(c) real([2*(1-c^2/cp1^2)^0.5 2*mu2/mu1*(1-c^2/cp2^2)^0.5 (2-c^2/cs1^2) mu2/mu1*(2-c^2/cs2^2);
%             (c^2/cs1^2-2) mu2/mu1*(2-c^2/cs2^2) -2*(1-c^2/cs1^2)^0.5 2*mu2/mu1*(1-c^2/cs2^2)^0.5;
%             1 -1 (1-c^2/cs1^2)^0.5 -(1-c^2/cs2^2)^0.5;
%             (1-c^2/cp1^2)^0.5 (1-c^2/cp2^2)^0.5 1 1]);
           
    D=@(c) real([k,-s1(c),-k,-s2(c);
        -r1(c),k,-r2(c),-k;
        (r1(c).^2*(2*mu1+lambda1)-k^2.*lambda1),-s1(c).*2*mu1*k,-((2*mu2+lambda2)*r2(c)^2-k^2.*lambda2),-s2(c)*2*mu2*k;
     -2*mu1*r1(c)*k, mu1*(s1(c)^2+k^2), -2*mu2*r2(c)*k, -mu2*(s2(c)^2+k^2)]);
     dmat_f=@(c) det(D(c));

    c=fzero(dmat_f,c_0);


    DMAT=D(c);

    CMAT=DMAT(2:end,2:end);
    Cload=DMAT(2:end,1);
    
    tempvec=-CMAT\Cload;
    %tempvec=null(DMAT);
    A=1; B=tempvec(1); C=tempvec(2); D=tempvec(3);

    % Substitute in the obtained wavespeed

    r1=r1(c); r2=r2(c); s1=s1(c); s2=s2(c);

    % Assemble the fields and associated tractions

    syms x y t 

    u1=(A.*k.*exp(-r1.*(y-y_l))-s1.*B.*exp(-s1.*(y-y_l))).*cos(k.*(x-c.*t));
    u2=(-r1.*A.*exp(-r1.*(y-y_l))+B.*k.*exp(-s1.*(y-y_l))).*sin(k.*(x-c.*t));
 
    u=matlabFunction([u1;u2]);
   
    v1=(C.*k.*exp(r2.*(y-y_l))+s2.*D.*exp(s2.*(y-y_l))).*cos(k.*(x-c.*t));
    v2=(r2.*C.*exp(r2.*(y-y_l))+D.*k.*exp(s2.*(y-y_l))).*sin(k.*(x-c.*t));
 
    v_tt=diff([v1;v2],t,2);
    v=matlabFunction([v1;v2]);
    v_tt=matlabFunction(v_tt);
    
    
    u1=matlabFunction(u1);
    u2=matlabFunction(u2);

    v1=matlabFunction(v1);
    v2=matlabFunction(v2);



    
end