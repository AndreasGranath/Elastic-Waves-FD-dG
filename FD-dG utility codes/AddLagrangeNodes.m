function v=AddLagrangeNodes(x,order)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION:  Pads a vector x with Lagrange nodes depending on the input
%               order
%
% INPUT : - x: Vector x of linear nodes which we want to convert to dG with
%             Lagrange nodes determined by the order
%
%         - order: the order of Lagrange nodes added, greater than 1
%
% OUTPUT: - v: the padded vector
%
% AUTHOR:  Andreas Granath (andreas.granath@umu.se)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    n=length(x);
    weights=1/order:1/order:1-1/order;
    v=zeros(1,(n-1)*(order+1));
    v(1:(order+1):end-order)=x(1:(end-1)); 
    v(order+1:(order+1):end)=x(2:end);
    
    for i=1:n-1
        v(2+(i-1)*(order+1):order+(i-1)*(order+1))=x(i)+weights.*(x(i+1)-x(i));
    end
end