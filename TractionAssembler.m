

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION: Assembles symbolic traction and load on a square domain given
%              analytic fields u,v and specified constant material
%              parameters
%              
% INPUT : - u,v: functions @(t,x,y) displacement fields in x and y directions
%
%         - mu,lambda,rho: scalars descirbing the Lamé parameters and
%                          density
%
% OUTPUT: - TN,TS,TE,TW: functions @(t,x,y) describing the traction along
%                        the north, south, east and west boundaries
%
%         - TN_tt,TS_tt,TE_tt,TW_tt: second time derivatives of the
%                                    traction along the north, south, east and west boundaries
%
%         - f,f_tt: the load resulting from the chosen fields u,v and
%                   material parameters
%
% AUTHOR:  Andreas Granath (andreas.granath@umu.se)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TN,TN_tt,TS,TS_tt,TE,TE_tt,TW,TW_tt,f,f_tt]=TractionAssembler(u,v,mu,lambda,rho)
    syms x y t
    % Assemble normal

 
    u=sym(u);
    v=sym(v);

 

    e1_tt=matlabFunction(diff(u,t,2)); e2_tt=matlabFunction(diff(v,t,2));


    u_x=diff(u,x); u_y=diff(u,y); v_x=diff(v,x); v_y=diff(v,y);

    % Use parametrized normal N(t)=[-h'(t),1]/sqrt(1+h'(t)^2) for top surface
    
    
    s11=(2*mu+lambda)*u_x+lambda*v_y;
    s22=(2*mu+lambda)*v_y+lambda*u_x;
    s12=mu*(u_y+v_x);

    %s11_S=(2*mu+lambda)*us_x+lambda*vs_y;
    %s22_S=(2*mu+lambda)*vs_y+lambda*us_x;
    %s12_S=mu*(us_y+vs_x);

    % Traction on northern interface
    N1=0;
    N2=1;
    t1=s11*N1+s12*N2;
    t2=s12*N1+s22*N2;
    TN=[t1;t2];
    TN_tt=diff(TN,t,2);

   % Nms=[N1;N2];

    % Traction on southern interface


    N1=0;
    N2=-1;


    t1=s11*N1+s12*N2;
    t2=s12*N1+s22*N2;
    TS=[t1;t2];
    TS_tt=diff(TS,t,2);

        Nms(1:2,2)=[N1;N2];

    % Traction on western interface
    N1=-1;
    N2=0;
    t1=s11*N1+s12*N2;
    t2=s12*N1+s22*N2;
    TW=[t1;t2];
    TW_tt=diff(TW,t,2);

            %Nms(1:2,3)=[N1;N2];

    % Traction on eastern interface
    N1=1;
    N2=0;
    t1=s11*N1+s12*N2;
    t2=s12*N1+s22*N2;
    TE=[t1;t2];
    TE_tt=diff(TE,t,2);

            %Nms(1:2,4)=[N1;N2];



    TN=matlabFunction(TN); TN_tt=matlabFunction(TN_tt);
    TS=matlabFunction(TS); TS_tt=matlabFunction(TS_tt);
    TW=matlabFunction(TW); TW_tt=matlabFunction(TW_tt);
    TE=matlabFunction(TE); TE_tt=matlabFunction(TE_tt);
   

    % Assemble load
    f1=rho*diff(u,t,2)-mu*(diff(u,x,2)+diff(u,y,2))-(mu+lambda)*(diff(u,x,2)+diff(v_x,y,1));
    f2=rho*diff(v,t,2)-mu*(diff(v,x,2)+diff(v,y,2))-(mu+lambda)*(diff(u_x,y,1)+diff(v,y,2));
    f=matlabFunction([f1;f2]);
    f_tt=matlabFunction(diff([f1;f2],t,2));
    f_tttt=matlabFunction(diff([f1;f2],t,4));



   



end