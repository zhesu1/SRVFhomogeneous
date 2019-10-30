function [q2_gamma, gamma]=PDSM_get_q_gammaRep(q1,q2)


[dimM,~,T]=size(q1);

tvec=linspace(0,1,T);

% map to R^n*n
qq1 = zeros(dimM^2,T);
qq2n = zeros(dimM^2,T);
for i=1:dimM
    for j=1:dimM
        qq1((i-1)*dimM+j,:)=q1(i,j,:);
        qq2n((i-1)*dimM+j,:)=q2(i,j,:);
    end
end
%
% mex DynamicProgrammingQ.c
[gamma1]=DynamicProgrammingQ(qq1,qq2n,0,0);
gamma=invertGamma(gamma1);
gamma=(gamma-gamma(1))/(gamma(end)-gamma(1));

% figure;
% plot(tvec,gamma);

%% get (q2,gamma)

q2_gamma = zeros(dimM,dimM,T);
qq2ngamma = zeros(dimM^2, T);
for i=1:dimM^2
    qq2ngamma(i,:)=interp1(tvec, qq2n(i,:),gamma,'linear');
end

for i=1:dimM
    for j=1:dimM
    q2_gamma(i,j,:)= qq2ngamma((i-1)*dimM+j,:).*sqrt(gradient(gamma,1/(T-1)));
    end
end