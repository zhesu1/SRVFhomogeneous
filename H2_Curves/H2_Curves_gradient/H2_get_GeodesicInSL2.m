function [alphab] = H2_get_GeodesicInSL2(f1,q1,f2,q2,S)

[m ,~ ,T]=size(q1);

q = zeros(m, m , T, S+1);
for j=1:T
    for k=1:S+1
        q(:,:,j,k)=q1(:,:,j)*(1-(k-1)/S)+q2(:,:,j)*(k-1)/S;
    end
end

% shooting vector in G
v=H2_invRieExpOnSLn(f1\f2);

StartingPts = zeros(m, m ,S+1);
for i=1:S+1
    StartingPts(:,:,i)=f1*H2_RieExpOnSLn((i-1)*v/S);
end

% gammab is the geodesic from gamma1 to gamma2
alphab = zeros(m, m , T+1, S+1);
for j=1:(S+1)
        alphab(:,:,:,j) = inv_Qfun(StartingPts(:,:,j),q(:,:,:,j));
end