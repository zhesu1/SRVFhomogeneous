clear;
T = 100; %spatial discretization of the curve
m = 3;
t = linspace(0,1,T);

A1 = rand(m, m, T);
A2 = rand(m, m, T);

for i = 1: T
    A1(:,:,i) = A1(:,:,i)-diag(diag(A1(:,:,i)));
    A2(:,:,i) = A2(:,:,i)-diag(diag(A2(:,:,i)));
    
    % get alphas in SL(3)
    alpha1(:,:,i) = expm(A1(:,:,i)); 
    alpha2(:,:,i) = expm(A2(:,:,i));
end

%% get curves in PDSMs

for i = 1: T
    beta1(:,:,i) = sqrtm(alpha1(:,:,i)*alpha1(:,:,i)');
    beta2(:,:,i) = sqrtm(alpha2(:,:,i)*alpha2(:,:,i)');
end

%% Calculate Geodesic. 
% -) Iter decides how often dynamic programming and optimization over starting
%    point are iterated
% -) optrigid decides if rigid motions are factored out
% -) doplot decides if the minimizing geodesic is plotted

tic
% [dis, v, qv, gamma_opt] =...
    PDSM_Main_Geodesic(beta1, beta2,'iter',1,'optrigid',0,'doplot',1);
toc
