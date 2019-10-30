function beta=S2_projection(alpha)

[dimM,~,T]=size(alpha);
n=zeros(1,dimM);n(1,dimM)=1;n=n';

beta = zeros(dimM, T);
for i=1:T
    beta(:,i)=alpha(:,:,i)*n;
end
