function [PA2B,PB2A]= Copy_of_make_projection_good_bad(N,order)
    %x=[-1:2/(N-1):1];
    H=SBP4(N,2/(N-1));
    M=1/((N-1))*kron(eye(N-1),Copy_of_oneDimMassMatrix(@LagrangeRbf,order+2));

    load('good_bad_Coeffs.mat');
    
    Pa2b=zeros((order+1)*(N-1),N);
    Pa2b(1:Nb2,1:Nb1)=gb_Block1;
    Pa2b(end-Nb2+1:end,end-Nb1+1:end)=gb_Block2;
    Nr=N-2*Nb1;
    for i=1:Nr
        Pa2b(No+(order+1)*(i-1)+1:No+Nh+(order+1)*(i-1),Nb1+i)=Interior;
    end
    q=zeros(2*Nb1*Nb2+Nh,1);
    q(1:Nb1*Nb2)=reshape(gb_Block1,[Nb1*Nb2,1]);
    q(Nb1*Nb2+1:Nb1*Nb2+Nh)=Interior;
    q(Nb1*Nb2+Nh+1:end)=reshape(gb_Block2,[Nb1*Nb2,1]);

    Pb2a=H\((M*Pa2b)');

  
    opts = optimset('Algorithm', 'levenberg-marquardt','Display','off');
    x = lsqnonlin(@(x) objective(q,x,Nb1,Nb2,Nh,Nr,M,H,order,N,No) ,q,[],[],opts);

    PB2A=H\(M*PA2B)';
end

function r=objective(q,x,Nb1,Nb2,Nh,Nr,M,H,order,N,No)
   
    Pa2b=zeros((order+1)*(N-1),N);
    Pa2b(1:Nb2,1:Nb1)=reshape(x(1:Nb1*Nb2),[Nb2,Nb1]);
    Pa2b(end-Nb2+1:end,end-Nb1+1:end)=reshape(x(Nb1*Nb2+Nh+1:end),[Nb2,Nb1]);
    
    for i=1:Nr
        Pa2b(No+(order+1)*(i-1)+1:No+Nh+(order+1)*(i-1),Nb1+i)=x(Nb1*Nb2+1:Nb1*Nb2+Nh);
    end

    Pb2a=H\(M*Pa2b)';

    e=sort(real(eig(full(Pa2b*Pb2a))));
    r=diff(e);
     plot(e,'*')    
  drawnow
end