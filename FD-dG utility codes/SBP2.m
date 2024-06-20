%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2:a  ordn. SBP Finita differens         %%%
%%% operatorer framtagna av Ken Mattsson    %%%
%%%                                         %%%
%%% H           (Normen)                    %%%
%%% D1=H^(-1)Q  (approx f?rsta derivatan)   %%%
%%% D2          (approx andra derivatan)    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H,D1,D2,S] = SBP2(N,h)
%N=11;h=1/10;
H = speye(N,N);
H(1,1) = 1/2;
H(end,end) = 1/2;
H = H*h;
D2 = spalloc(N,N,3*N); D1 = spalloc(N,N,3*N);
D2(1,1:3) = [1 -2 1]; D1(1,1:2) = [-1,1];
D2(end,end-2:end) = [1 -2 1]; D1(N,N-1:N) = [-1,1];
for i = 2:N-1
    D2(i,i) = -2;
    D2(i,i-1) = 1;
    D2(i,i+1) = 1;
    D1(i,i-1) = -1/2;
    D1(i,i+1) = 1/2;
end
D2 = D2/h/h; D1 = D1/h;
S = speye(N,N);
S(1,1:3) = [3/2 -2 1/2]*(-1);
S(end,end-2:end) = [1/2 -2 3/2];
S = S/h;
E0 = sparse(N,N); E0(1,1) = 1;
EN = sparse(N,N); EN(N,N) = 1;
B = -E0 + EN;
M = -H*D2+B*S;