function [q2_gamma,gamma]=S2_get_q_gamma_Rep(q1,q2)


[dimM,~,T]=size(q1);

tvec=linspace(0,1,T);
% reparametrization for P(S^n)
% map to R^n

qq1 = zeros(dimM-1, T);
qq2n = zeros(dimM-1, T);
for i=1:(dimM-1)
    for j=1:T
        qq1(i,j)=q1(i,dimM,j);
        qq2n(i,j)=q2(i,dimM,j);
    end
end
%
%% mex DynamicProgrammingQ.c

[gamma1]=DynamicProgrammingQ(qq1,qq2n,0,0);
gamma=invertGamma(gamma1);
gamma=(gamma-gamma(1))/(gamma(end)-gamma(1));

% figure;
% plot(tvec,gamma);

%% get (q2,gamma)

q2_gamma = zeros(dimM,dimM,T);
qq2ngamma = zeros(dimM-1, T);
for i=1:dimM-1
    qq2ngamma(i,:)=interp1(tvec, qq2n(i,:),gamma,'linear');
end

%%-----------------------------------

for i=1:dimM-1
    q2_gamma(i,dimM,:)= qq2ngamma(i,:).*sqrt(gradient(gamma,1/(T-1)));
    q2_gamma(dimM,i,:)=-qq2ngamma(i,:).*sqrt(gradient(gamma,1/(T-1)));
end