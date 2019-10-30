function y=S2_startPoint(x,n)

[dimM,~]=size(x);

theta=2*pi*rand;
ytheta=[cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];


if eq(x(:,1),-n)
    y(:,:,1)=eye(dimM);
    y(1,1,1)=-1;y(dimM,dimM,1)=-1;
else
    y(:,:,1)=(eye(dimM)-2*(n+x(:,1))*(n+x(:,1))'/...
        (norm(n+x(:,1)))^2)*(eye(dimM)-2*(n*n'))*ytheta;
end