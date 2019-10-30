function y = inv_Qfun(f1,q1)

[~,~,T]=size(q1);

dT=1/(T-1);

v1 = zeros(size(q1));
for j=1:T
        v1(:,:,j)=q1(:,:,j).*norm(q1(:,:,j),'fro')*dT;
end

ExpOfqs = zeros(size(v1));
for i=1:T
        ExpOfqs(:,:,i)=expm(v1(:,:,i));
end

y(:,:,1)=f1;

for i=2:(T+1)
    y(:,:,i)=y(:,:,i-1)*ExpOfqs(:,:,i-1);
end

