function [A,M,P,T,Traction,Mbdry,Nglob]=dG_data(p,tdG,mu,lambda,h,tau,order)


% Determine neighbours
Nbrs=Tri2Tri(p,tdG);

% Convert P1 regular FEM grid to dG grid


[P,T]=ChangeOrderFromLinear(order,p,tdG);

% Determine number of points on each element as well as the total number of
% points and elements

nPts=length(T(:,1));
nn=length(P(1,:));
nt=length(T(1,:));


% Initialize 1D and 2D quadratures

%eqpts=[-sqrt(3/5) 0 sqrt(3/5)]; % Gaussian edge quadratures, exact up to order 5
%eqweights=[5/9 8/9 5/9];
%eqpts=[-sqrt(3/7+2/7*sqrt(6/5)),-sqrt(3/7-2/7*sqrt(6/5)),sqrt(3/7-2/7*sqrt(6/5)),sqrt(3/7+2/7*sqrt(6/5))];
%eqweights=[(18-sqrt(30))/36,(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36];

eqpts=[-1/3*sqrt(5+2*sqrt(10/7)) -1/3*sqrt(5-2*sqrt(10/7)) 0 1/3*sqrt(5-2*sqrt(10/7)) 1/3*sqrt(5+2*sqrt(10/7))];
eqweights=[(322-13*sqrt(70))/900 (322+13*sqrt(70))/900 128/225 (322+13*sqrt(70))/900 (322-13*sqrt(70))/900];

neqpts=length(eqpts); % Number of quadrature points along edges

[qwgts,rspts]=quadratureTriSymmetric(8); % 2D quadrature

% Initialize bulk matrices
K=sparse(2*nn,2*nn);
M=sparse(2*nn,2*nn);
F=sparse(2*nn,2*nn);
Nglob=sparse(2*nn,1);
b=sym(zeros(2*nn,1));
Traction=sparse(2*nn,2*nn);

Mbdry=sparse(nn,nn);
%mu=1; lambda=1;
D=[2*mu+lambda lambda 0;lambda 2*mu+lambda 0;0 0 mu];

% Edge to node lookup table (order dependent,currently for P2 elements)
%edge2node=[2 4 3;3 5 1;1 6 2];
%edge2node=[2 3;3 1;1 2];

if order==3
    edge2node=[2,4,5,3;3,6,7,1;1,8,9,2];
    x_edge=[1,2/3,1/3,0;0,0,0,0;0,1/3,2/3,1];
    y_edge=[0 1/3 2/3 1;1 2/3 1/3 0;0 0 0 0];
elseif order==4
    edge2node=[2,4,5,6,3;3,7,8,9,1;1,10,11,12,2];
    x_edge=[1,3/4,2/4,1/4,0;0 0 0 0 0;0,1/4,2/4,3/4,1];
    y_edge=[0,1/4,2/4,3/4,1;1,3/4,2/4,1/4,0;0,0,0,0,0];
end

for i=1:nt % Total element loop

    KK=sparse(2*nPts,2*nPts);
    MK=sparse(nPts,nPts);

    % local-to-global map of the points at the current element
    l2g=T(:,i);

    % Determine nodes bounding the current element
    xv=P(1,1+nPts*(i-1):3+nPts*(i-1)); yv=P(2,1+nPts*(i-1):3+nPts*(i-1));

    % Construct map to reference element
    Xcoeffs=[ones(3,1),xv',yv']\[0;1;0];

    Ycoeffs=[ones(3,1),xv',yv']\[0;0;1];
  


    for q=1:length(qwgts)

        % Determine 2D quadrature points as well as the basis functions and
        % derivatives at those points

        r=rspts(q,1); s=rspts(q,2);  
        if order==3
             [S,dSdx,dSdy,detJ]=Isopmap(xv,yv,r,s,@P3Shapes);
        elseif order==4
              [S,dSdx,dSdy,detJ]=Isopmap(xv,yv,r,s,@P4Shapes);
        end
        wxarea=qwgts(q)*detJ/2;
    
        % Initialize strain matrix to use for stiffness matrix
        B=sparse(3,2*nPts);
        B(1,1:nPts)=dSdx'; 
        B(2,nPts+1:end)=dSdy'; B(3,1:nPts)=dSdy';
        B(3,nPts+1:end)=dSdx';

        % Compute local stiffness and mass matrices using the Voight
        % formulation

        KK=KK+(B'*D)*B.*wxarea;
        MK=MK+(S*S').*wxarea;

       
    end

    % Total degrees of freedom, including the contribution of the negative
    % element

    dofs=[l2g' [l2g']+nn];
    
    % Add local mass contribution to the total mass matrix
    M(l2g,l2g)=M(l2g,l2g)+MK;
    M(l2g+nn,l2g+nn)=M(l2g+nn,l2g+nn)+MK;

    % Add local stiffness contribution to the global stiffness matrix
    K(dofs,dofs)=K(dofs,dofs)+KK;

    for j=1:3 % Loop over edges bounding element K_i

        % Determine neighbouring element K_nbr on edge E_j of K_i
        nbr=Nbrs(i,j);

        % Find the node indices from lookuptable for edge E_j
        e2n=edge2node(j,:);
        xe=x_edge(j,:); ye=y_edge(j,:);

        % Check if E_j is internal edge or domain boundary, in which case
        % nbr<0

        if nbr<0 % Treat the case when E_j is on the domain boundary
           
            nbr=i;

            % Initialize local traction
            TractionLoc=sparse(2*nPts,2*nPts);
    

                

           % Determine normal of E_j and its length
            dr=P(:,e2n(end)+nPts*(i-1))-P(:,e2n(1)+nPts*(i-1));
           edgeLength=sqrt(dr'*dr);
           nloc=[dr(2),-dr(1)]./edgeLength;

          % Initialize normal finding matrix for Voight formulation
           N=[nloc(1) 0 nloc(2);0 nloc(2) nloc(1)];
          

           % Determine traction matrix Tu (note, not boundary value!) on element boundary and add it to global
           % traction
           
           for q=1:length(e2n)
             
               if order==3
                  [~,dSdx,dSdy,~]=Isopmap(xv,yv,xe(q),ye(q),@P3Shapes);  
                  Bloc=[dSdx',sparse(1,nPts);sparse(1,nPts),dSdy';dSdy',dSdx'];
               elseif order==4
                  [~,dSdx,dSdy,~]=Isopmap(xv,yv,xe(q),ye(q),@P4Shapes);  
                  Bloc=[dSdx',sparse(1,nPts);sparse(1,nPts),dSdy';dSdy',dSdx'];
               end

                   Traction_Temp=N*D*Bloc;
                   TractionLoc(e2n(q),1:nPts)=Traction_Temp(1,1:nPts);
                   TractionLoc(e2n(q),nPts+1:2*nPts)=Traction_Temp(1,nPts+1:2*nPts);
                   TractionLoc(nPts+e2n(q),1:nPts)=Traction_Temp(2,1:nPts);
                   TractionLoc(e2n(q)+nPts,nPts+1:2*nPts)=Traction_Temp(2,nPts+1:end);

                   % Bdry flux
                           
                                      
           end



           % Assemble local traction integral for boundary value \int
           % \sigma(u)\cdot n ds (vanishing traction)

           
           Traction(l2g,l2g)=TractionLoc(1:nPts,1:nPts);
           Traction(l2g,nn+l2g)=TractionLoc(1:nPts,nPts+1:end);
           Traction(nn+l2g,l2g)=TractionLoc(nPts+1:end,1:nPts);
           Traction(nn+l2g,nn+l2g)=TractionLoc(nPts+1:end,nPts+1:end);
           
           %   BdryMat(dofs,dofs)=BdryMat(dofs,dofs)+BdryTractionLoc;

          [r,s]=GaussianMapping(j);       
          rr=r(eqpts); ss=s(eqpts);
      
           b1=zeros(nPts,1); b2=zeros(nPts,1);
           Mbdry_loc=zeros(nPts,nPts);

           for q=1:neqpts
         
         % total weight for integration
         weight=eqweights(q)*edgeLength/2;

         % Evaluate basis functions and its derivatives on positive and
         % negative element

         if order==3
            [S,~,~,~]=Isopmap(xv,yv,rr(q),ss(q),@P3Shapes);     
         elseif order==4
            [S,~,~,~]=Isopmap(xv,yv,rr(q),ss(q),@P4Shapes);       
         end

   

        %C=[S',sparse(1,nPts);sparse(1,nPts),S'];
        %b1=b1+S*T1(rr(q),ss(q),nloc(1),nloc(2),t)*weight;
        %b2=b2+S*T2(rr(q),ss(q),nloc(1),nloc(2),t)*weight;

        Mbdry_loc=Mbdry_loc+S*S'*weight;
       
       
     
       

      
           end
           b(l2g,1)=b1; b(l2g+nn,1)=b2;
           Mbdry(l2g,l2g)=Mbdry(l2g,l2g)+Mbdry_loc;

           % Output global Normal vector (Nx2)
           Nglob(l2g)=nloc(1).*ones(length(l2g),1);
           Nglob(l2g+nn)=nloc(2).*ones(length(l2g),1);
        

        end  
        if nbr>i % If an internal edge, compute the local flux contributions

          % Determine degrees of freedom for neighbouring element
          negdofs=T(:,nbr);

           % Determine normal of E_j and its length
           dr=P(:,e2n(end)+nPts*(i-1))-P(:,e2n(1)+nPts*(i-1));
           edgeLength=sqrt(dr'*dr);
           nloc=[dr(2),-dr(1)]./edgeLength;
           
     
           % Determine nodes on the neighbouring element
           xnv=P(1,[1:3]+nPts*(nbr-1));
           ynv=P(2,[1:3]+nPts*(nbr-1));
       
   
          
           % Determine corresponding edge on negative element, start by
           % finding the indices of the common nodes from the negative and
           % positive elements [Maybe a bit slow]   

           
           
            c1=intersect(find(abs(xnv-xv(e2n(1)))<100*eps),find(abs(ynv-yv(e2n(1)))<100*eps));
            c2=intersect(find(abs(xnv-xv(e2n(end)))<100*eps),find(abs(ynv-yv(e2n(end)))<100*eps));
             coords=[c1, c2];

         
           % Use the indices found to determine the edge index jn on the
           % negative element to use in the mapping of weights.

              if sort(coords)==[1,3]
                            jn=2;
             elseif sort(coords)==[1,2]
                            jn=3;
              else
                            jn=1;

              end

     
         % Construct maps from quadratures on [-1,1] to reference elment
         % edge, keeping in mind that there is a difference when mapping
         % from the negative element, resulting in [rn,sn], compared to
         % from the positive elment, resulting in [r,s]. Note that this has
         % to respect the fact that the quadrature points as seen from
         % K_nbr are in the reversed order as when seen from K_j

         [r,s]=GaussianMapping(j);
         [rn,sn]=GaussianMapping(jn);
        
         rn=@(t) rn(-t); sn=@(t) sn(-t);
         rr=r(eqpts); ss=s(eqpts);
         rrn=rn(eqpts);ssn=sn(eqpts);

        % Initialize normal finding matrix for Voight formulation
        N=[nloc(1) 0 nloc(2);0 nloc(2) nloc(1)];

        % Initialize local flux matrices
        FPos1=sparse(2*nPts,2*nPts);  FPos2=sparse(2*nPts,2*nPts);
        FNeg1=sparse(2*nPts,2*nPts);  FNeg2=sparse(2*nPts,2*nPts);

        FPos=zeros(2*nPts,2*nPts); FNeg=zeros(2*nPts,2*nPts);
        for q=1:neqpts
         
         % total weight for integration
         weight=eqweights(q)*edgeLength/2;

         % Evaluate basis functions and its derivatives on positive and
         % negative element

         if order==3
            [SP,dSPdx,dSPdy,~]=Isopmap(xv,yv,rr(q),ss(q),@P3Shapes);
            [SN,dSNdx,dSNdy,~]=Isopmap(xnv,ynv,rrn(q),ssn(q),@P3Shapes);
         elseif order==4
            [SP,dSPdx,dSPdy,~]=Isopmap(xv,yv,rr(q),ss(q),@P4Shapes);
            [SN,dSNdx,dSNdy,~]=Isopmap(xnv,ynv,rrn(q),ssn(q),@P4Shapes);
         end

        % Initialize matrices corresponding to the jump in basis functions,
        % Cpos,Cneg and the local strain matrices B,Bneg.

        Cpos=[SP',sparse(1,nPts);sparse(1,nPts),SP'];
        Cneg=[SN',sparse(1,nPts);sparse(1,nPts),SN'];

        B=[dSPdx',sparse(1,nPts);sparse(1,nPts),dSPdy';dSPdy',dSPdx'];
        Bneg=[dSNdx',sparse(1,nPts);sparse(1,nPts),dSNdy';dSNdy',dSNdx'];


        % Assemble local flux terms using the Voight formulation, where
        % sigma(u)=DBu and sigma(u)*n=NDBu.

        FPos1=FPos1+0.5*(Cpos'*N*D*B+(N*D*B)'*Cpos-2*tau*(2*mu+lambda)/h*Cpos'*Cpos)*weight;
        FNeg1=FNeg1+0.5*(Cpos'*N*D*Bneg-(N*D*B)'*Cneg+2*tau*(2*mu+lambda)/h*Cpos'*Cneg)*weight;
        FPos2=FPos2+0.5*(-Cneg'*N*D*B+(N*D*Bneg)'*Cpos+2*tau*(2*mu+lambda)/h*Cneg'*Cpos)*weight;
        FNeg2=FNeg2-0.5*(Cneg'*N*D*Bneg+(N*D*Bneg)'*Cneg+2*tau*(2*mu+lambda)/h*Cneg'*Cneg)*weight;   
                                                 
        end

        % Determine total degrees of freedom on the current edge E_j
        nbrdofs=[negdofs' negdofs'+nn];

        %Add the local flux contribution to the global flux matrix
        F(dofs,dofs)=F(dofs,dofs)+FPos1;
        F(dofs,nbrdofs)=F(dofs,nbrdofs)+FNeg1;
        F(nbrdofs,dofs)=F(nbrdofs,dofs)+FPos2;
        F(nbrdofs,nbrdofs)=F(nbrdofs,nbrdofs)+FNeg2;

                     
        end
       
    end

end
A=(-K+F);

end