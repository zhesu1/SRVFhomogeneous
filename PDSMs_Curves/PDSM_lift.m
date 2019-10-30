function y = PDSM_lift(x)

% lift

[~,~,T]=size(x);

y(:,:,1)=x(:,:,1);


for i=2:T
    y(:,:,i)=y(:,:,i-1)*sqrtm(inv(y(:,:,i-1))*x(:,:,i)...
        *x(:,:,i)'*inv(y(:,:,i-1)'));
end