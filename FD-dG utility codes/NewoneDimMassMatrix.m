% SUMMARY: This code generates the 1D mass matrix along the interval [-1,1]
% using given basis functions as specified by funs, e.g. Lagrange basis
% functions. It uses Gaussian quadrature to evaluate it, so note that the
% variable order has to be adapted so that it is sufficiently high

% INPUT:
% - funs: a column vector function containing the basis functions phi_i so
%         that M_{ij} is obtained by integrating phi_i(x)*phi_j(x). The
%         size of funs will also be used to determine the dimensions of M.
%
% - order: a scalar determining the order of quadrature used. Note that
%          an n-point Gaussian quadrature is exact for polyomials of order 
%          up to 2*n-1. Thus if we have a basis of order b, we integrate a
%          polynomial of order 2*b so we need to choose order=b+1 for the
%          method to be exact. It currently handles order up to 4.


function M=NewoneDimMassMatrix(funs,order)
    v= funs;
    M=zeros(length(v(0,order-1)),length(v(0,order-1)));

    switch order
        case 1
            xg=0; wg=2;
        case 2
            xg=[-1/sqrt(3) 1/sqrt(3)]; wg=[1 1];
        case 3
            xg=[-sqrt(3/5) 0 sqrt(3/5)]; wg=[5/9 8/9 5/9];
        case 4
            xg=[-sqrt(3/7+2/7*sqrt(6/5)) -sqrt(3/7-2/7*sqrt(6/5)) sqrt(3/7-2/7*sqrt(6/5)) sqrt(3/7+2/7*sqrt(6/5))];
            wg=[(18-sqrt(30))/36 (18+sqrt(30))/36 (18+sqrt(30))/36 (18-sqrt(30))/36];
        case 5
            xg=[-1/3*sqrt(5+2*sqrt(10/7)) -1/3*sqrt(5-2*sqrt(10/7)) 0 1/3*sqrt(5-2*sqrt(10/7)) 1/3*sqrt(5+2*sqrt(10/7))];
            wg=[(322-13*sqrt(70))/900 (322+13*sqrt(70))/900 128/225 (322+13*sqrt(70))/900 (322-13*sqrt(70))/900];
        case 6
            xg=[-0.9324695142031521 -0.6612093864662645 -0.2386191860831969 ...
                0.2386191860831969 0.6612093864662645 0.9324695142031521];
            wg=[0.1713244923791704 0.3607615730481386 0.4679139345726910 ...
                 0.4679139345726910 0.3607615730481386 0.1713244923791704];
        case 7
            xg=[-0.9491079123427585 -0.7415311855993945 -0.4058451513773972 0.0000000000000000 ...
                0.4058451513773972 0.7415311855993945 0.9491079123427585];

            wg=[0.1294849661688697 0.2797053914892766 0.3818300505051189 0.4179591836734694 ...
                	0.3818300505051189 0.2797053914892766 0.1294849661688697];

    end
    
    for i=1:order
        M=M+wg(i)*v(xg(i),order-1)*v(xg(i),order-1)';
    end
