% GAUSSIAN_QUAD computes the Nth order Gaussian Quadrature Points and Weights
% [x,w] = Gaussian_quad(N)
%
% inputs:
%   N: quadrature order
%
% outputs:
%   x: quadrature points
%   w: quadrature weights
%
% Based on routines given in
%   @BOOK{HesthavenWarburton2008,
%     title = {Nodal Discontinuous {G}alerkin Methods: {A}lgorithms, Analysis, and
%    Applications},
%     publisher = {Springer},
%     year = {2008},
%     author = {Hesthaven, Jan S. and Warburton, Tim},
%     volume = {54},
%     series = {Texts in Applied Mathematics},
%     doi = {10.1007/978-0-387-72067-8}
%   }
function [x,w] = Gaussian_quad(N)
  if(N==0)
    x = 0;
    w = 2;
  else
    h1 = 2*(0:N);
    Je = 2./(h1(1:N) + 2).*(1:N).*(1:N)./sqrt((h1(1:N)+1).*(h1(1:N)+3));
    A = spdiags([Je,0;0,Je]',[-1 1],N+1,N+1);
    [V,x] = eig(full(A));
    x = diag(x);
    w = 2*(V(1,:)').^2;
  end
end
