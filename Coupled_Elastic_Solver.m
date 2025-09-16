function  [TractiondG,ABulk,xx,X,yy,Y,P,Hx,Hy,HH,Mgamma,M_EW,MI,T,Ex,EN,ES,EdG_E,EdG_W,EdG_S,EdG,NDoFs,tdG,pdG,Mbdry,Nglob,Pu2wG] = Coupled_Elastic_Solver(hdG,hFD,hy,mat_param,tau,GammaCoord,mesh_FEM,order,x_l,y_l,shift,ProblemType)
    mu=mat_param(1); lambda=mat_param(2);  rho1=mat_param(3);
    mu2=mat_param(4); lambda2=mat_param(5); rho2=mat_param(6);


    h1=hFD; h2=hy;


    xFD=0:h1:x_l; 

    NxFD=length(xFD); 

    yFD=GammaCoord:hy:y_l;
    ny=length(yFD); nx=NxFD;
    NyFD=ny;
    xx=kron(xFD',ones(NyFD,1)); yy=kron(ones(NxFD,1),yFD');
    X=reshape(xx,[NyFD,NxFD]); Y=reshape(yy,[NyFD,NxFD]);
    NFD=length(xx);

    % Initialize identity matrices
    Ix=speye(NxFD); Iy=speye(NyFD);
    Ixy=kron(Ix,Iy);

    % Construct quadrature matrix and differential operators, since we use full
    % compatibility, boundary derivative operator S is not needed!
    if order==3
          [hx,d1x,d2x,sx]=SBP4(NxFD,h1);
          [hy,d1y,d2y,~]=SBP4Full(NyFD,h2);
    else
         [hx,d1x,d2x,sx]=SBP6(NxFD,h1);
         [hy,d1y,d2y,~]=SBP6Full(NyFD,h2);
    end

    Sx=kron(speye(2),kron(sx,Iy)); Sy=kron(speye(2),kron(Ix,d1y));
    D1x=kron(speye(2),kron(d1x,Iy)); D2x=kron(speye(2),kron(d2x,Iy));
    D1y=kron(speye(2),kron(Ix,d1y)); D2y=kron(speye(2),kron(Ix,d2y));
    Dxy=kron(speye(2),kron(d1x,d1y));

    Hx=kron(speye(2),kron(hx,Iy)); Hy=kron(speye(2),kron(Ix,hy));
    HH=Hx*Hy;

    hxi=sparse(diag(diag(1./hx))); hyi=sparse(diag(diag(1./hy)));

    HI=kron(speye(2),kron(hxi,hyi));

    % Initialize boundary-finding matrices on FD size

    e0y=zeros(1,length(yFD)); e0y(1,1)=1; eS=e0y'*e0y;
    eny=zeros(1,length(yFD)); eny(1,end)=1; eN=eny'*eny;
    e0x=zeros(1,length(xFD)); e0x(1,1)=1; eW=e0x'*e0x;
    enx=zeros(1,length(xFD)); enx(1,end)=1; eE=enx'*enx;



    EW=sparse(kron(speye(2),kron(eW,speye(NyFD)))); EE=sparse(kron(speye(2),kron(eE,speye(NyFD))));
    EN=sparse(kron(speye(2),kron(speye(NxFD),eN)));
    ES=sparse(kron(speye(2),kron(speye(NxFD),eS)));

    Ex=EE-EW; Ey=EN;
    %% 
    %%%%%%%%%%%%%%%%%%%%% Initialize material matrices %%%%%%%%%%%%%%%%%%%%

    L1=[(2*mu+lambda)*Ixy, 0*Ixy;0*Ixy, mu*Ixy];
    L3=[mu*Ixy, 0*Ixy;0*Ixy, (2*mu+lambda)*Ixy];
    L2=[0*Ixy, lambda*Ixy;mu*Ixy,0*Ixy];

    %%%%%%%%%%%%%%%%%%%%%% FD discretization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    disp('Asembling FD bulk')

    DD=L1*D2x+L3*D2y+(L2+L2')*Dxy;

    Tx=(L1*Sx+L2*D1y);
    Ty=(L3*D1y+L2'*D1x);


    % Construct SAT terms for the traction along the boundaries of FD domain

    SATb1=-(Hx\(Ex*Tx));
    SATb2=-(Hy\(Ey*Ty));


    AFD=DD+SATb1+SATb2;

    chi1=2*mu2+lambda2; chi2=2*mu+lambda;
    alpha=tau/(chi1); beta=tau/(chi2);


    % Construct analytical geometry objects so that it is possible to e.g.
    % remove sections of the mesh if desired


      [pdG,edG,tdG]=meshToPet(mesh_FEM);
      tdG(4,:)=[];
    
    if shift~=0
        InterfaceNumber=intersect(find(pdG(1,:)==0.8),find(pdG(2,:)==GammaCoord));
        pdG(1,InterfaceNumber)=pdG(1,InterfaceNumber)-shift*hdG;
    end

    %%%%%%%%%%%%%%%%%%%% Assemble dG bulk matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assemble bulk matrix A=-K+F where K is the stiffness matrix and
    % F is the flux matrix. Also assemble 2xNumDoFs matrix P of nodal values
    % and triangulation matrix T modified according to the desired order

    disp('Asembling dG bulk')


    [AdG,MdG,P,T,TractiondG,Mbdry,Nglob]=dG_data(pdG,tdG,mu2,lambda2,hdG,alpha,order);

    NDoFs=length(P(1,:));
    npts=(order+1)*(order+2)/2;

     for i=1:length(MdG(1,:))/npts
         MI(npts*(i-1)+1:npts*i,1+npts*(i-1):npts*i)=inv(MdG(npts*(i-1)+1:npts*i,npts*(i-1)+1:npts*i));
     end

    %MI=inv(MdG);
    AdG=MI*AdG;

    % Modify edge indices depending of what experiment you want to perform

    
        if ProblemType == "Gaussian"
            EdgeIndices=[20,1,2,19]; % Gaussian pulse in UMU geometry
        elseif ProblemType == "Trigonometric" 
            EdgeIndices=[1,2,3,4];   % Trigonometric example
        elseif ProblemType == "Stoneley"
            EdgeIndices=[1,2,3,4];   % Stoneley example
        elseif ProblemType ==  "CurvedBdry"
            EdgeIndices=[2,1,4,3];   % Curved geometry example
        else
            disp('Undefined experiment')
 
        end

    % Construct matrices picking up DOFS along the boundaries
    [EdG,EdG_t,permInds,indicesOnInterface,EdG_E,EdG_S,EdG_W,bdryIndices,n_E,n_S]=AssembleBoundaryFindingMatrix(edG,tdG,P,T,order,EdgeIndices);

    EdG_E=kron(eye(2),EdG_E); EdG_W=kron(eye(2),EdG_W); 
    EdG_S=kron(eye(2),EdG_S);

    EdG=kron(eye(2),EdG);
    EdG_t=kron(eye(2),EdG_t);
    EFD=ES;

 % Extract indices for boundaries
    i_E=nonzeros(bdryIndices(1,:)); i_S=nonzeros(bdryIndices(2,:));
    i_W=nonzeros(bdryIndices(3,:));
   % Assemble quadrature matrices  

   n_interface=length(indicesOnInterface)/(order+1);

    mE=hdG/2*kron(eye(n_E),NewoneDimMassMatrix(@LagrangeRbf,order+1));
    mgamma=hdG/2*kron(eye(n_interface),NewoneDimMassMatrix(@LagrangeRbf,order+1));
    mS=hdG/2*kron(eye(n_S),NewoneDimMassMatrix(@LagrangeRbf,order+1));


     Mgamma=sparse(NDoFs,NDoFs); 
     Mgamma(indicesOnInterface(permInds),indicesOnInterface(permInds))=mgamma;
     Mgamma(i_S,i_S)=mS;
     Mgamma=kron(eye(2),Mgamma);
    
     M_EW=sparse(NDoFs,NDoFs);
     M_EW(i_E,i_E)=mE; M_EW(i_W,i_W)=mE;
     M_EW=kron(eye(2),M_EW);

    % Construct projection operators
   
    [pw2uG,pw2uB,pu2wG,pu2wB]=ConstructFullProjs(0:h1:x_l,AddLagrangeNodes(0:hdG:x_l,order),order);

    
    pw2uG=pw2uG/x_l;pw2uB=pw2uB/x_l;pu2wG=pu2wG/x_l;pu2wB=pu2wB/x_l;

    %Set up full size projection operators
    Pw2uG=sparse(NDoFs,NFD);Pw2uB=sparse(NDoFs,NFD);
    Pu2wG=sparse(NFD,NDoFs); Pu2wB=sparse(NFD,NDoFs);

    Pu2w=sparse(NFD,NDoFs); Pw2u=sparse(NDoFs,NFD);

    Ptemp=sparse(NDoFs,NDoFs);

    Pw2uG(indicesOnInterface(permInds),1:NyFD:(NxFD-1)*NyFD+1)=pw2uG;
    Pw2uB(indicesOnInterface(permInds),1:NyFD:(NxFD-1)*NyFD+1)=pw2uB;
   


    Pu2wG(1:NyFD:(NxFD-1)*ny+1,indicesOnInterface(permInds))=pu2wG;
    Pu2wB(1:NyFD:(NxFD-1)*ny+1,indicesOnInterface(permInds))=pu2wB;

    Pw2uG=kron(eye(2),Pw2uG); Pw2uB=kron(eye(2),Pw2uB);
    Pu2wG=kron(eye(2),Pu2wG); Pu2wB=kron(eye(2),Pu2wB);


    % Assemble traction operators
    TFDy=-ES*Ty; TractiondG=TractiondG*EdG_t;
   
    % Determine min_K h_K
    hdGmin=mesh_FEM.MinElementSize;

    % Assemble SAT and numerical flux terms
    SAT_dG_c=0.5*HI*(-(TFDy)'*Hx*Pu2wG*EdG+2*alpha*chi1/hdGmin*EFD*Hx*Pu2wB*EdG+...
        +2*beta*chi2/h1*EFD*Hx*Pu2wG*EdG);

    SAT_FD_c=0.5*HI*((TFDy)'*Hx*EFD-2*alpha*chi1/hdGmin*EFD*Hx*Pu2wB*Pw2uG*EFD+...
        -2*beta*chi2/h1*EFD*Hx*EFD);

    FLUX_dG_c=0.5*MI*((TractiondG)'*Mgamma*EdG-2*alpha*chi1/hdGmin*EdG*Mgamma*EdG...
        -2*beta*chi2/h1*EdG*Mgamma*Pw2uB*Pu2wG*EdG);

    FLUX_FD_c=0.5*MI*(-(TractiondG)'*Mgamma*Pw2uG*EFD+2*alpha*chi1/hdGmin*EdG*Mgamma*Pw2uG*EFD+2*beta*chi2/h1*EdG*Mgamma*Pw2uB*EFD);

    BULK_c1=[1/rho1*SAT_FD_c,  1/rho1*SAT_dG_c; 1/rho2*FLUX_FD_c, 1/rho2*FLUX_dG_c];

    %-------------------------------------------%

    SAT_dG_t=-0.5*HI*(EFD*Hx*Pu2wB*TractiondG);

    SAT_FD_t=-0.5*HI*(EFD*Hx*TFDy);

    FLUX_dG_t=0.5*MI*(EdG*Mgamma*TractiondG);

    FLUX_FD_t=-0.5*MI*(EdG*Mgamma*Pw2uB*TFDy);


    BULK_t=[1/rho1*SAT_FD_t, 1/rho1*SAT_dG_t; 1/rho2*FLUX_FD_t, 1/rho2*FLUX_dG_t];
   
    % Assemble total bulk matrix for semidiscretization

    ABulk=[1/rho1*AFD, sparse(2*NxFD*NyFD,2*NDoFs); sparse(2*NDoFs,2*NxFD*NyFD), 1/rho2*AdG]+BULK_c1+BULK_t;

end