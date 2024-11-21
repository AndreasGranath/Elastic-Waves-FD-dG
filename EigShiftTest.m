% Set up a rectangular domain [0,x_l] x [0, y_l] partitioned into the two
% subdomains [0,x_l] x [0,GammaCoord] and [0,x_l] x [GammaCoord,1] 

GammaCoord=1/5;
 x_l=1; y_l=1;

% Initialize dG geometry

rect=[3 4 0 x_l x_l 0 GammaCoord GammaCoord 0 0]';
gd=rect;
sf = 'rect';
ns = char('rect');
ns = ns';

g = decsg(gd,sf,ns);

model = createpde(1);
gg=geometryFromEdges(model,g);

% Set material parameters and densities
mu1=1; lambda1=1; rho1=1;
mu2=1; lambda2=1; rho2=1;

mat_param=[mu1,lambda1,rho1,mu2,lambda2,rho2];

% Determine wave speeds
v_s=max([sqrt((2*mu1+lambda1)/rho2),sqrt((2*mu2+lambda2)/rho2)]);

% Set penalty term used both in the interior of the dG scheme and in the
% semi-discretizations at the interface, note that it is scaled by the
% maximal wavespeed

 tau=150;
 

 % Initialize vector for errors
Num_Refinements=12;
 SystemSizeFD=zeros(1,Num_Refinements);
 SystemSizedG=zeros(1,Num_Refinements);
Error=zeros(1,Num_Refinements); order=3;

Hq_dG=zeros(1,Num_Refinements); Hq=zeros(1,Num_Refinements);

dGL2Error=zeros(1,Num_Refinements); 

FDError=zeros(1,Num_Refinements);

EOC=zeros(1,Num_Refinements-1);

h=[1/10,1/15,1/20,1/30];
shift=[1e-6,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,1e-1];

Spectrum=zeros(4,11);
 % Loop over mesh refinements

 for q=1:4
  for s=1:11
 % Set mesh sizes, discetize FD side with nx nodes, currently matching 
 % element sizes


 hFD=h(q); hdG=hFD;

 nx=1/h(q)+1;
 hy=(y_l-GammaCoord)/(nx-1);

 % Determine number of nodes along y-direction on FD side


 % Assemble dG mesh using the mesh size hdG

 mesh_FEM=generateMesh(model,'Hmax',hdG,'GeometricOrder','linear');

% Call on coupled elastic solver 

[TdG,ABulk,xx,X,yy,Y,P,Hx,Hy,HH,Mgamma,M_EW,MI,T,Ex,EN,ES,EdG_E,EdG_W,EdG_S,EdG,NDoFs,tdG,pdG,Mbdry,Nglob] = Coupled_Elastic_Solver(hdG,hFD,hy,mat_param,tau,GammaCoord,mesh_FEM,order,x_l,y_l,shift(s));
 
 Spectrum(q,s)=h(q)^2*max(abs(eig(full(ABulk))))
 end
 end