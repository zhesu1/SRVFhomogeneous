function q = PDSM_to_q(p)

[m,~,T] = size(p);

dT=1/(T-1);

dp(:,:,1)=(p(:,:,2)-p(:,:,1))/dT;
dp(:,:,T)=(p(:,:,T)-p(:,:,T-1))/dT;
for j=2:(T-1)
    dp(:,:,j)=(p(:,:,j+1)-p(:,:,j-1))/(2*dT);
end

v = zeros(m, m, T);
q = zeros(m, m, T);
L = zeros(1, T);
for i=1:T
    v(:,:,i)=inv(p(:,:,i))*dp(:,:,i); % do left translation 
    v(:,:,i)=v(:,:,i)-trace(v(:,:,i))*eye(m,m)/m; % project onto sl(n,R)
    L(i) = sqrt(norm(v(:,:,i),'fro'));
    if L(i) > 0
        q(:,:,i)=v(:,:,i)/L(i);
    else
        q(:,:,i) = zeros(m, m);
    end
end