addpath('FD-dG Projection operators\')
addpath('FD-dG utility codes\')

%% Construct FD geometry
clear

%% Construct dG geometry
GammaCoord=1/5;
order=3;
US=zeros(9,30000);
US2=zeros(9,30000);

% Stoneley experiments performed
% tau=35*sqrt(2*mu+lambda), dt=0.011*h

 x_l=1; y_l=2;

% Points to monitor, n1=50, n2= 325 (n2 is where pulse is initiated)
rect=[3 4 0 x_l x_l 0 GammaCoord GammaCoord 0 0 zeros(1,16)]';

 rect1=[2 4 0.25 0.35 0.35 0.25 0.163 0.163 0.0793 0.0793 zeros(1,16)]';
 rectMid=[2 4 0.275 0.325 0.325 0.275 0.163 0.163 0.095 0.095 zeros(1,16)]';
 MidCirc=[1 0.3 0.1043 0.05/2 0 0 0 0 0 0 zeros(1,16)]';

   Bcirc1=[4 0.3 0.0793 0.025 0.05 pi/2 zeros(1,20)]';

  rect2=[2 4 0.65 0.75 0.75 0.65 0.163 0.163 0.0793 0.0793 zeros(1,16)]';
 rectMid2=[2 4 1-0.275 1-0.325 1-0.325 1-0.275 0.163 0.163 0.095 0.095 zeros(1,16)]';
 MidCirc2=[1 0.7 0.1043 0.05/2 0 0 0 0 0 0 zeros(1,16)]';

Bcirc2=[4 0.7 0.0793 0.025 0.05 pi/2 zeros(1,20)]';

 Mpoly=[2 12 0.4435 0.4623 0.5 0.537 0.5565 0.5565 0.5339 0.5339 0.5 0.4661 0.4661 0.4435 ...
     0.163 0.163 0.125 0.163 0.163 0.0498 0.0498 0.1252 0.0875 0.1252 0.0498 0.0498]';



 gd=[rect,rect1,Bcirc1,rect2,Bcirc2,Mpoly,rectMid,rectMid2];
 ns=char('rect','rect1','rect2','Bcirc1','Bcirc2','Mpoly','rectMid','rectMid2');

 ns = ns';

 %sf = '((rect-Bcirc1-Bcirc2-Mpoly-ellipsetest)-rect1-rect2+rectMid+rectMid2-MidCirc2+MidCirc2+MidCirc)';
 %ssf='rectMid+rectMid2';
 sf = 'rect-(rect1+Bcirc1)-(rect2+Bcirc2)-Mpoly+(rectMid+rectMid2)';
 g = decsg(gd,sf,ns);

 model = createpde(1);
 gg=geometryFromEdges(model,g);


Qs=[20];
 Error=zeros(1,10);

 % Define problem type
 ProblemType = 'Gaussian';

 for q=1:1

      
 % Generate P^1 mesh 
 %nx=21;
 nx=10*Qs(q)+1;
 hFD=x_l/(nx-1); hdG=hFD;
 hy=hFD;
 ny=length([GammaCoord:hy:y_l]);
 mesh_FEM=generateMesh(model,'Hmax',hdG,'GeometricOrder','linear');

 % Set penalty parameter


% Initialize physical parameters

mu1=0.1; lambda1=1; rho1=1;
mu2=0.01; lambda2=0.1; rho2=1;


mat_param=[mu1,lambda1,rho1,mu2,lambda2,rho2];

 tau=60;



 

 tend=3; dt=0.02*hdG; NumTimeSteps=round(tend/dt);



e1=@(t,x,y) 550.*exp(-1000.*((x-0.4).^2+(y-0.135).^2)+0.01.*t)+550.*exp(-1000.*((x-0.6).^2+(y-0.135).^2)+0.01.*t); 
e2=@(t,x,y) 550.*exp(-1000.*((x-0.4).^2+(y-0.135).^2)+0.01.*t)+550.*exp(-1000.*((x-0.6).^2+(y-0.135).^2)+0.01.*t);

uex=@(t,x,y) [e1(t,x,y);e2(t,x,y)];


disp("Determine analytic traction")

[TN_FD,TN_FD_tt,TS_FD,~,TE_FD,TE_FD_tt,~,~,f_FD,f_FD_tt]=TractionAssembler(e1,e2,mu1,lambda1,rho1);
[TN_dG,~,TS_dG,TS_dG_tt,TE_dG,TE_dG_tt,~,~,f_dG,f_dG_tt]=TractionAssembler(e1,e2,mu2,lambda2,rho2);




[TdG,ABulk,xx,X,yy,Y,P,Hx,Hy,HH,Mgamma,M_EW,MI,T,Ex,EN,ES,EdG_E,EdG_W,EdG_S,EdG,NDoFs,tdG,pdG,Mbdry] = Coupled_Elastic_Solver(hdG,hFD,hy,mat_param,tau,GammaCoord,mesh_FEM,order,x_l,y_l,0,ProblemType);
Survey_index=intersect(find(yy==0.4),find(xx==0.9));





Mbdry=kron(eye(2),Mbdry);

%u2=real([S_FD(dt,xx,yy);S_dG(dt,P(1,:)',P(2,:)')]); 
%u1=real([S_FD(0,xx,yy);S_dG(0,P(1,:)',P(2,:)')]);

Tr_FD= @(t) real(Hx\(Ex*TE_FD(t,xx,yy))+Hy\((EN)*TN_FD(t,xx,yy)));
 Tr_dG= @(t) real(MI*((EdG_E-EdG_W)*Mbdry*TE_dG(t,P(1,:)',P(2,:)')+(EdG_S)*Mbdry*TS_dG(t,P(1,:)',P(2,:)')));

 Tr_FD_tt= @(t) real(Hx\(Ex*TE_FD_tt(t,xx,yy))+Hy\((EN)*TN_FD_tt(t,xx,yy)));
 Tr_dG_tt= @(t) real(MI*((EdG_E-EdG_W)*Mbdry*TE_dG_tt(t,P(1,:)',P(2,:)')+(EdG_S)*Mbdry*TS_dG_tt(t,P(1,:)',P(2,:)')));



fload=@(t) real(0*[1/rho1*Tr_FD(t);1/rho2*Tr_dG(t)]+0.*[1/rho1*f_FD(t,xx,yy);1/rho2*f_dG(t,P(1,:)',P(2,:)')]);
 fd_load=@(t) real(0*[1/rho1*Tr_FD_tt(t);1/rho2*Tr_dG_tt(t)]+0.*[1/rho1*f_FD_tt(t,xx,yy);1/rho2*f_dG_tt(t,P(1,:)',P(2,:)')]);



u2=real([uex(dt,xx,yy);uex(dt,P(1,:)',P(2,:)')]); 
u1=real([uex(0,xx,yy);uex(0,P(1,:)',P(2,:)')]);

% for i=1:length(tdG(1,:))
%     p1_indices(1+3*(i-1):3*i)=[1:3]+10*(i-1);
% end
% 
% for i=1:length(T(1,:))
%     P_linear(:,1+3*(i-1):3*i)=full(P(:,1+10*(i-1):3+10*(i-1)));
%     T_linear(1:3,i)=[1:3]+3*(i-1);
% end

U_survey=zeros(1,NumTimeSteps-1);
TT=zeros(1,NumTimeSteps-1);





disp("Starting time stepping")
for i=1:NumTimeSteps-1
     tt=i*dt;
 
        k1=speye(size(ABulk))*u2+dt^2/12*ABulk*u2;
        k2=speye(size(ABulk))*fload(tt)+dt^2/12*ABulk*fload(tt);

      u3=2*u2-u1+dt.^2.*ABulk*k1+dt^2.*k2+dt^2*(dt^2/12*speye(size(ABulk)))*fd_load(tt);
      u1=u2; u2=u3;
      U_survey(i)=sqrt(u3(Survey_index)^2+u3(nx*ny+Survey_index)^2);
     
        TT(i)=tt+dt;

        if q==4           
            US(9,i)=sqrt(e1(tt+dt,xx(Survey_index),yy(Survey_index))^2+e2(tt+dt,xx(Survey_index),yy(Survey_index))^2);
            %US2(9,i)=sqrt(e1(tt+dt,xx(Survey_index2),yy(Survey_index2))^2+e2(tt+dt,xx(Survey_index2),yy(Survey_index2))^2);
        end
    

         

end
    
      US(2*q-1,1:NumTimeSteps-1)=U_survey;
      US(2*q,1:NumTimeSteps-1)=TT;       
      NTS(q)=NumTimeSteps-1;
     %q
   
end



%%
 %save("surfacesurvey1305.mat","mat_param","NTS","tau","US")
  % save("EOC__4th_order_Stoneley_0422.mat","EOC","Error","mat_param","tau","dt")h
