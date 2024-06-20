% plots the shape of one or several parametric triangles

% INPUTS
% p : [N x D] matrix of nodal coordinates for the parametrizised
%     surface, the number of nodes N implicitly defines the order,
%     and D is the spatial dimension
% t : [N x M] triangle connectivity matrix where N is the number of
%     triangles and M is the number of nodes per triangle which gives
%     the order of basis functions
% OUTPUTS
% -

function paramTriSolPlot(p,t,tFE,uFE,pResolution)

if (~exist('pResolution', 'var'))
    pResolution = 5;
end

pDim=size(p,2); % spatial dimension
%pOrder=2; % polynomial order of parametrization
%feOrder=2; % FE basis order

if pDim==2
    % creating a refined reference grid for triangle plotting
    [pbary,tbary]=triBaryGrid(pResolution);


hold on;
for n=1:size(t,1)

    % map pbary to Cartesian coordinates
    pElmDoFs=t(n,:);
    xnod=p(pElmDoFs,1);ynod=p(pElmDoFs,2);
    feElmDoFs=tFE(n,:);unod=uFE(feElmDoFs);
    
    x=zeros(size(pbary,1),1); y=x;u=x; 
    %z=x; 
    for i=1:size(pbary,1)
        l1=pbary(i,1); l2=pbary(i,2); 
        [pPHI,~,~]=P3Shapes(l1,l2);
        x(i)=pPHI'*xnod; y(i)=pPHI'*ynod; 
        %z(i)=pPHI*znod;
        fePHI=P3Shapes(l1,l2);
        u(i)=unod*fePHI;
    end

    % plot parametric triangles
    TR=triangulation(tbary,x,y);
    trisurf(tbary,x,y,u,'EdgeColor','none','FaceColor','interp');

    % extract and plot triangle edges
    E=freeBoundary(TR)';
    plot3(x(E),y(E),u(E),'k-','linewidth',1);
    %patchline(x(E),y(E),z(E),'linestyle','-','edgecolor','k','linewidth',1,'edgealpha',0.2);

end
hold off;


axis equal;
%light;

hold off
end