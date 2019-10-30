function znew=H2_vfzTovzInvexp(z,y)

[~,dim]=size(z);
f_z=z*expm(logm(z')-logm(z));

v=logm(inv(f_z)*y);

Delta=0.000001;
Basis=zeros(dim,dim,(dim-1)*dim/2);
vcoefficient = zeros(1, dim*(dim-1)/2);
for i=2:dim
    for j=1:i-1
        Basis(i,j,(i-1)*(i-2)/2+j)=1/sqrt(2);
        Basis(j,i,(i-1)*(i-2)/2+j)=-1/sqrt(2);
        vcoefficient((i-1)*(i-2)/2+j)=v(i,j);
    end
end

z0 = zeros(dim, dim, (dim-1)*dim/2);
Df = zeros(dim, dim, (dim-1)*dim/2);
DF = zeros(dim, dim, (dim-1)*dim/2);
for i=1:(dim-1)*dim/2
    z0(:,:,i)=z*expm(Delta*Basis(:,:,i));
    Df(:,:,i)=inv(f_z)*(z0(:,:,i)*...
        expm(logm(z0(:,:,i)')-logm(z0(:,:,i)))-f_z)/Delta;
    DF(:,:,i)=(Df(:,:,i)-Df(:,:,i)')/2; %project to so(dim,R)
end

Amatrix = zeros((dim-1)*dim/2,(dim-1)*dim/2);
for i=1:(dim-1)*dim/2
    for j=1:(dim-1)*dim/2
        Amatrix(i,j)=trace(DF(:,:,j)*Basis(:,:,i)');
    end
end

prevcoe=inv(Amatrix)*vcoefficient';

prev=zeros(dim,dim);

for i=2:dim
    for j=1:i-1
        prev(i,j)=prevcoe((i-1)*(i-2)/2+j)/sqrt(2);
        prev(j,i)=-prevcoe((i-1)*(i-2)/2+j)/sqrt(2);
    end
end

znew=z*expm(prev);
