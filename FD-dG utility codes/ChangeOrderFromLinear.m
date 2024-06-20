
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION:  Converts matlab [p,~,t] structure into a dG structure
%               [P,~,T] of polynomial order given by the input 'order'
%              
%
% INPUT : - order: the desired dG polynomial order of the method
%
%         - p,t: linear point data matrix (p) and triangulation matrix (t)      
%                from the matlab [p,e,t] triplet. 
%
% OUTPUT: - P,T: a dG point data matrix (P) and triangulation matrix (T)
%                corresponding to a method of order given by 'order'
%
% AUTHOR:  Andreas Granath (andreas.granath@umu.se)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P,T]=ChangeOrderFromLinear(order,p,t)

  
% Make sure that the 4th row in the linear triangulation matrix is
  % removed

  if length(t(:,1))==4
    t(4,:)=[];
  end

    % Set total number of points per element

    nPoints=(order+1)*(order+2)/2;
 
   % Initialize dG point- and triangulation matrix

    P=zeros(2,length(t(1,:))*nPoints);
    T=zeros(nPoints,length(t(1,:)));


    % Construct reference mapping

    x_ref=[0,1,0]; y_ref=[0,0,1];
    refmat=[ones(3,1),x_ref',y_ref'];

    % Set up reference element depending on order

    if order==1
        X_r=[0,1,0];
        Y_r=[0,0,1];
    elseif order==2
        X_r=[0,1,0,1/2,0,1/2];
        Y_r=[0,0,1,1/2,1/2,0];
    elseif order==3
        X_r=[0,1,0,2/3,1/3,0,0,1/3,2/3,1/3];
        Y_r=[0,0,1,1/3,2/3,2/3,1/3,0,0,1/3];
    elseif order==4

    X_r=[0,1,0,3/4,2/4,1/4,0,0,0,1/4,2/4,3/4,1/4,1/4,2/4];
    Y_r=[0,0,1,1/4,2/4,3/4,3/4,2/4,1/4,0,0,0,2/4,1/4,1/4];
    end

    % Traverse each element

    for i=1:length(t(1,:))

       % Determine vertices on element
        locInds=t(:,i);

        xloc=p(1,locInds); yloc=p(2,locInds);

       % Map to reference element
        c_x=refmat\xloc'; c_y=refmat\yloc';

        P(1,1+nPoints*(i-1):nPoints*i)=([ones(length(X_r),1),X_r',Y_r']*c_x)';
        P(2,1+nPoints*(i-1):nPoints*i)=([ones(length(X_r),1),X_r',Y_r']*c_y)';

         T(1:nPoints,i)=1+(i-1)*nPoints:i*nPoints';        
    end

    T=sparse(T); P=sparse(P);
    

end
