function plotEllipse(A, C)

N = 50; 

[~,D,V] = svd(A);

a = sqrt(D(1,1));
b = sqrt(D(2,2));
c = sqrt(D(3,3));
[X,Y,Z] = ellipsoid(0,0,0,a,b,c,N);

XX = zeros(N+1,N+1);
YY = zeros(N+1,N+1);
ZZ = zeros(N+1,N+1);
for k = 1:length(X)
    for j = 1:length(X)
        point = [X(k,j) Y(k,j) Z(k,j)]';
        P = V * point;
        XX(k,j) = P(1)+C(1);
        YY(k,j) = P(2)+C(2);
        ZZ(k,j) = P(3)+C(3);
    end
end

mesh(XX,YY,ZZ);
axis equal

