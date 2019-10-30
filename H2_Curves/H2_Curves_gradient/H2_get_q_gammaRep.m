function [q2_gamma,gamma]=H2_get_q_gammaRep(q1,q2)


[dimM,~,T]=size(q1);

tvec=linspace(0,1,T);

% map to R^n
qq1 = zeros(dimM, T);
qq2n = zeros(dimM, T);
for i=1:dimM
    for j=1:T
        qq1(i,j)=q1(i,dimM,j);
        qq2n(i,j)=q2(i,dimM,j);
    end
end
%mex DynamicProgrammingQ.c

[gamma1]=DynamicProgrammingQ(qq1,qq2n,0,0);
gamma=invertGamma(gamma1);
gamma=(gamma-gamma(1))/(gamma(end)-gamma(1));

%%
% figure;
% plot(tvec,gamma);

%% get (q2,gamma)

q2_gamma = zeros(dimM,dimM,T);
qq2ngamma = zeros(dimM, T);
for i=1:dimM
    qq2ngamma(i,:)=interp1(tvec, qq2n(i,:),gamma,'linear');
end

q2_gamma(1,2,:)= qq2ngamma(1,:).*sqrt(gradient(gamma,1/(T-1)));
q2_gamma(2,2,:)= qq2ngamma(2,:).*sqrt(gradient(gamma,1/(T-1)));
q2_gamma(2,1,:)= qq2ngamma(1,:).*sqrt(gradient(gamma,1/(T-1)));
q2_gamma(1,1,:)=-qq2ngamma(2,:).*sqrt(gradient(gamma,1/(T-1)));