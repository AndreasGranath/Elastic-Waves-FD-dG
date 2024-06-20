
function err=NewL2Error(p,t,fh,fh2,f,f2,order)

    nel=length(t(1,:));               % Number of elements in the mesh
    err2=0 ;                         % Initialize quadratic L2 error
    nPts=length(t(:,1));                    % Number of points per element
    
    
  
    n=8;
    [wvec,cmat]=quadratureTriSymmetric(n);


     for i=1:nel
         l2g=t(:,i);
         x=p(1,l2g); y=p(2,l2g);
         xv=x(1:3);yv=y(1:3);
        
         Mat=[xv(2)-xv(1) xv(3)-xv(1) xv(1);
              yv(2)-yv(1) yv(3)-yv(1) yv(1) ];
        
    
            for j=1:length(wvec)
             r=cmat(j,1); s=cmat(j,2);
             if order==1
                 [S,~,~,detJ]=Isopmap(xv,yv,r,s,@P1Shapes);
             elseif order==3
             [S,~,~,detJ]=Isopmap(xv,yv,r,s,@P3Shapes);
            elseif order==4
                [S,~,~,detJ]=Isopmap(xv,yv,r,s,@P4Shapes);
            end



              globCoords=Mat*[r;s;1];
            Uh1=fh(1+(i-1)*nPts:i*nPts)'*S;
            Uh2=fh2(1+(i-1)*nPts:i*nPts)'*S;
            U1=f(globCoords(1),globCoords(2));
            U2=f2(globCoords(1),globCoords(2));
            S1=(Uh1-U1).^2;
            S2=(Uh2-U2).^2;
      
            wxarea=detJ/2*wvec(j);
                   
             
            err2 = err2 + (S1+S2)*wxarea;
           end
     end
    
   err=err2;
end