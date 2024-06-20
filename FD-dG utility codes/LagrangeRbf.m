

function f=LagrangeRbf(x,order)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION: determine Lagrange basis functions in 1D over the reference
%              domain [-1,1] 
%
% INPUT : - x: a scalar value at which the functions should be evaluated
%
%         - order: the desired order of the Lagrange basis
%
% OUTPUT: - f: a vector of size 'order'+1 with the values of the basis
%              functions at point x
%
% AUTHOR:  Andreas Granath (andreas.granath@umu.se)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if order==1
        % Explicit form of the polynomial
        phi1=@(x) -0.5.*(x-1);
        phi2=@(x) 0.5.*(x+1);

        % Evaluate at x
        f=[phi1(x);phi2(x)];
    elseif order==2
        % Explicit form of the polynomial
        phi1=@(x) 1/2.*x.*(x-1);
        phi2=@(x) -(x.^2-1);
        phi3=@(x) 1/2.*x.*(x+1);
        
        % Evaluate at x
        f=[phi1(x);phi2(x);phi3(x)];

    elseif order==3
          % Explicit form of the polynomial
            phi1=@(x) -9/16.*(x.^2-1/9).*(x-1);
            phi2= @(x) 27/16*(x.^2-1).*(x-1/3); 
            phi3=@(x) -27/16*(x.^2-1).*(x+1/3);
            phi4=@(x) 9/16*(x.^2-1/9).*(x+1);

          % Evaluate at x
          f=[phi1(x);phi2(x);phi3(x);phi4(x)];

    elseif order==4
        % Explicit form of the polynomial
        phi1=@(x) 2/3.*x.*(x-1).*(x.^2-1/4);
        phi2=@(x) -8/3.*x.*(x-1/2).*(x.^2-1);
        phi3=@(x) 4.*(x.^2-1).*(x.^2-1/4);
        phi4=@(x) -8/3.*x.*(x.^2-1).*(x+1/2);
        phi5=@(x) 2/3.*x.*(x+1).*(x.^2-1/4);
   
        % Evaluate at x
        f=[phi1(x);phi2(x);phi3(x);phi4(x);phi5(x)];

    elseif order==5
        % Explicit form of the polynomial
        phi1=@(x) -5^4/(16*6*8).*(x-1).*(x.^2-1/25).*(x.^2-9/25);
        phi2=@(x) 5^5/(16*6*8).*(x.^2-1).*(x.^2-1/25).*(x-3/5);
        phi3=@(x) -5^5/(16*4*6)*(x.^2-1).*(x.^2-9/25).*(x-1/5);
        phi4=@(x) 5^5/(16*4*6)*(x.^2-1)*(x.^2-9/25).*(x+1/5);
        phi5=@(x) -5^5/(16*6*8)*(x.^2-1).*(x.^2-1/25).*(x+3/5);
        phi6=@(x) 5^4/(16*6*8).*(x+1).*(x.^2-1/25).*(x.^2-9/25);
    
        % Evaluate at x
         f=[phi1(x);phi2(x);phi3(x);phi4(x);phi5(x);phi6(x)];

    end

end