function q = S2_curveToq(p)

[m,~,T] = size(p);

dT=1/(T-1);

dp(:,:,1)=(p(:,:,2)-p(:,:,1))/dT;
dp(:,:,T)=(p(:,:,T)-p(:,:,T-1))/dT;
for j=2:(T-1)
    dp(:,:,j)=(p(:,:,j+1)-p(:,:,j-1))/(2*dT);
end

v = zeros(size(p));
q = zeros(size(p));
L = zeros(1, T);
for j=1:T
    v(:,:,j)=((p(:,:,j)'*dp(:,:,j))-(p(:,:,j)'*dp(:,:,j))')/2; 
    L(j) = sqrt(norm(v(:,:,j),'fro'));
    if L(j) > 0
        q(:,:,j)=v(:,:,j)/L(j);
    else
        q(:,:,j) = zeros(m, m);
    end
end