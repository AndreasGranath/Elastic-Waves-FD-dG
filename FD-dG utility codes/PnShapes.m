function [S,dSdr,dSds]=PnShapes(rr,ss)
r=[0,1,0,3/4,2/4,1/4,0,0,0,1/4,2/4,3/4,1/4,1/4,2/4];
s=[0,0,1,1/4,2/4,3/4,3/4,2/4,1/4,0,0,0,2/4,1/4,1/4];

f=@(r,s) [1,r,s,r.^2,r.*s,s.^2,r.^3,r.^2.*s,r.*s.^2,s.^3,r.^4,r.^3.*s,r.^2.*s.^2,r.*s.^3,s.^4];

for i=1:length(r)
    V(i,:)=f(r(i),s(i));
end

for i=1:length(r)
    c=zeros(length(r),1);
    c(i)=1;
    coeffs(i,:)=(V\c)';
end

S= coeffs*f(rr,ss)';


dfds= @(r,s) [0,0,1,0,r,2.*s,0,r.^2,2.*r.*s,3.*s.^2,0,r.^3,2.*s.*r.^2,3.*r.*s.^2,4.*s.^3]';
dfdr=@(r,s) [0,1,0,2.*r,s,0,3.*r.^2,2.*r.*s,s.^2,0,4.*r.^3,3.*r.^2.*s,2.*r.*s.^2,s.^3,0]';

dSdr=coeffs*dfdr(rr,ss); dSds=coeffs*dfds(rr,ss);
end