% Set up a rectangular domain [0,x_l] x [0, y_l] partitioned into the two
% subdomains [0,x_l] x [0,GammaCoord] and [0,x_l] x [GammaCoord,1] 

GammaCoord=1/5;
 x_l=1; y_l=1;


% Initialize dG geometry

model = createpde;
gg=geometryFromEdges(model,@bdry_fnc);
mesh_FEM=generateMesh(model);
% Set material parameters and densities
mu1=1; lambda1=1; rho1=1;
mu2=1; lambda2=1; rho2=1;

mat_param=[mu1,lambda1,rho1,mu2,lambda2,rho2];

% Determine wave speeds
v_s=max([sqrt((2*mu1+lambda1)/rho2),sqrt((2*mu2+lambda2)/rho2)]);

% Set penalty term used both in the interior of the dG scheme and in the
% semi-discretizations at the interface, note that it is scaled by the
% maximal wavespeed

 tau=150/v_s;

 % Initialize vector for errors
Num_Refinements=6; order=4;
Error=zeros(1,Num_Refinements);

Hq_dG=zeros(1,Num_Refinements); Hq=zeros(1,Num_Refinements);

dGL2Error=zeros(1,Num_Refinements); 

FDError=zeros(1,Num_Refinements);

EOC=zeros(1,Num_Refinements-1);

 % Loop over mesh refinements
 for q=1:Num_Refinements

 % Set mesh sizes, discetize FD side with nx nodes, currently matching 
 % element sizes

 nx=10*(q+1)+1;

 hFD=x_l/(nx-1); hdG=hFD;

 hy=(y_l-GammaCoord)/(nx-1);

 % Determine number of nodes along y-direction on FD side
 ny=length(GammaCoord:hy:y_l);

 % Assemble dG mesh using the mesh size hdG

 mesh_FEM=generateMesh(model,'Hmax',hdG,'GeometricOrder','linear');

% Call on coupled elastic solver 

[TdG,ABulk,xx,X,yy,Y,P,Hx,Hy,HH,Mgamma,M_EW,MI,T,Ex,EN,ES,EdG_E,EdG_W,EdG_S,EdG,NDoFs,tdG,pdG,Mbdry,Nglob] = Coupled_Elastic_Solver(hdG,hFD,hy,mat_param,tau,GammaCoord,mesh_FEM,order,x_l,y_l,0);
Mbdry=kron(eye(2),Mbdry);
% Assemble analytic Stoneley solution and obtain the wave speed c

w1=@(t,x,y) sin(3.*x+2.*y-t); w2=@(t,x,y) sin(2.*x+3.*y-2.*t);

uex=@(t,x,y) [w1(t,x,y);w2(t,x,y)];

% Set time stepping properties
 tend=1; dt=0.015*hdG; NumTimeSteps=round(tend/dt);



% Determine analytic traction and forcing

disp("Determine analytic traction")

 
 
%%---------------------- Assemble stress tensor ------------------------%

[su_11,su_12,su_22,su_11_tt,su_12_tt,su_22_tt,~,~]=StressTensorAssembler(w1,w2,mu1,lambda1,rho1);
[sv_11,sv_12,sv_22,sv_11_tt,sv_12_tt,sv_22_tt,~,~]=StressTensorAssembler(w1,w2,mu2,lambda2,rho2);

N1=Nglob(1:length(Nglob)/2); N2=Nglob(length(Nglob)/2+1:end);


Sv=@(t) Ebdry_v*([diag(sv_11(t,P(1,:),P(2,:)))*N1+diag(sv_12(t,P(1,:),P(2,:)))*N2;
        diag(sv_12(t,P(1,:),P(2,:)))*N1+diag(sv_22(t,P(1,:),P(2,:)))*N2]);

Sv_tt=@(t) Ebdry_v*([diag(sv_11_tt(t,P(1,:),P(2,:)))*N1+diag(sv_12_tt(t,P(1,:),P(2,:)))*N2;
        diag(sv_12_tt(t,P(1,:),P(2,:)))*N1+diag(sv_22_tt(t,P(1,:),P(2,:)))*N2]);

% Set up numerical traction and load for time stepping
%u2=real([S_FD(dt,xx,yy);S_dG(dt,P(1,:)',P(2,:)')]); 
%u1=real([S_FD(0,xx,yy);S_dG(0,P(1,:)',P(2,:)')]);

[TN_FD,TN_FD_tt,TS_FD,~,TE_FD,TE_FD_tt,~,~,f_FD,f_FD_tt]=TractionAssembler(w1,w2,mu1,lambda1,rho1);

Tr_FD= @(t) real(Hx\(Ex*TE_FD(t,xx,yy))+Hy\((EN)*TN_FD(t,xx,yy)));
%Tr_dG=  @(t) real(MI*((EdG_E-EdG_W)*Mbdry*TE_dG(t,P(1,:)',P(2,:)')+(EdG_S)*Mbdry*TS_dG(t,P(1,:)',P(2,:)')));
    Tr_dG=@(t) MI*(Mbdry*Sv(t));

 Tr_FD_tt= @(t) real(Hx\(Ex*TE_FD_tt(t,xx,yy))+Hy\((EN)*TN_FD_tt(t,xx,yy)));
 %Tr_dG_tt= @(t) real(MI*((EdG_E-EdG_W)*Mbdry*TE_dG_tt(t,P(1,:)',P(2,:)')+(EdG_S)*Mbdry*TS_dG_tt(t,P(1,:)',P(2,:)')));
Tr_dG_tt=@(t) MI*(Mbdry*Sv_tt(t));



fload=@(t) real([1/rho1*Tr_FD(t);1/rho2*Tr_dG(t)]+[1/rho1*f_FD(t,xx,yy);1/rho2*f_dG(t,P(1,:)',P(2,:)')]);
 fd_load=@(t) real([1/rho1*Tr_FD_tt(t);1/rho2*Tr_dG_tt(t)]+[1/rho1*f_FD_tt(t,xx,yy);1/rho2*f_dG_tt(t,P(1,:)',P(2,:)')]);



u2=real([S_FD(dt,xx,yy);S_dG(dt,P(1,:)',P(2,:)')]); 
u1=real([S_FD(0,xx,yy);S_dG(0,P(1,:)',P(2,:)')]);

% Obtain p1 nodes from the dG triangulation for plotting purposes only
% for i=1:length(tdG(1,:))
%     p1_indices(1+3*(i-1):3*i)=[1:3]+10*(i-1);
% end
% 
% for i=1:length(T(1,:))
%     P_linear(:,1+3*(i-1):3*i)=full(P(:,1+10*(i-1):3+10*(i-1)));
% end

% Perform explicit time stepping 

disp("Starting time stepping")
% 
for i=1:NumTimeSteps-1

    % Set current time step
     tt=i*dt;
 

      k1=ABulk*u2+dt^2/12*ABulk*(ABulk*u2);

      k2=fload(tt)+dt^2/12*(ABulk*fload(tt)+fd_load(tt));

      u3=2*u2-u1+dt.^2.*k1+dt^2*k2;
      u1=u2; u2=u3;

      if max(u3)>1.1
        disp("Solution has diverged")
        break
      end

%  Plot the solution att time tt+dt

%  surf(X,Y,reshape(u3(1:nx*ny),[ny,nx]))
%
%  hold on
%
%  yline(GammaCoord,'LineWidth',2)
%     
%  trisurf(T_linear',full(P_linear(1,:)),full(P_linear(2,:)),full(u3(2*nx*ny+p1_indices)))
%
%  xlim([0 1])
%  ylim([0 1])
%  set(gca, 'color', 'none');
%  set(gca,'LooseInset',get(gca,'TightInset'));
%  set(gca,'ytick',[])
%  shading interp
%  colormap jet
%  view(2)
%  drawnow

                 
                 
end

% Save current mesh size for plotting and calculating EOC
Hq(q)=hFD; Hq_dG(q)=hdG;

 

% Calculate discrete FD and dG errors, combine into a total discrete error

FDError(q) = sqrt(real((u3(1:2*nx*ny)-S_FD(tt+dt,xx,yy)))'*HH*(real(u3(1:2*nx*ny)-S_FD(tt+dt,xx,yy))));
dGL2Error(q)=sqrt(NewL2Error(P,T,u3(2*nx*ny+1:2*nx*ny+NDoFs),u3(2*nx*ny+NDoFs+1:end),@(x,y) vdG1(tt+dt,x,y),@(x,y) vdG2(tt+dt,x,y),order));
Error(q)=sqrt(FDError(q)^2+dGL2Error(q)^2);

% From second refinement and onwards, calculate the EOC of the method

   if q>1    
     EOC(q)=log(Error(q-1)/Error(q))/log(Hq_dG(q-1)/Hq_dG(q));
   
    end
 end

%save("CurvedEOC_4thOrder.mat","mat_param","tau","dt","tend","Error","EOC")

