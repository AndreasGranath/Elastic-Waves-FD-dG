
function [pbary,tbary]=triBaryGrid(n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION: creates a triangle grid of a triangle with grid coordinates
%              exoressed in barycentric coordinates (used for plotting)
%
% INPUT :  -n: a scalar defining the number of subdivisions of a triangle
%              side
%
%
% OUTPUT: - pbary,tbary: point matrix (pbary) and triangulation matrix
%           (tbary) of the barycentric triangulation
%
% AUTHOR:  Fredrik Bengzon, Mats G Larson
%          From the book "The Finite Element Method" by Larson and Bengzon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pbary=zeros(0.5*(n+1)*(n+2),3);

% generate barycentric points
c=0;
for i=0:n
    for j=0:n-i
        c = c + 1;
        pbary(c,:) = [(n-i-j)/n,j/n,i/n];
    end
end

% generate triangles
tbary=[];
c=0;
for i=0:n-1
    c=c+1;
    tbary(end+1,:)=[c,c+1,c+1+(n-i)];
    for j=1:n-i-1
        c=c+1;
        tbary(end+1,:)=[c,c+1+(n-i),c+(n-i)];
        tbary(end+1,:)=[c,c+1,c+1+(n-i)];
    end
    c=c+1;
end

end