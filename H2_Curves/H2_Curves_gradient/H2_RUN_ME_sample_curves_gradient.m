clear;
T=100; %spatial discretization of the curve
tvec=linspace(0,1,T);

%% Select example curves
%%%%%%%%%%%%%%%
%Pair 1

% beta1 = [2*tvec+3; 0.5*ones(1,T)];
% beta2 = [2*ones(1,T); 2*tvec+0.5];

%%%%%%%%%%%%%%%%
% Pair 2

% t1 = linspace(0,1/4,251);
% t2 = linspace(1/4,1/2,251);
% t3 = linspace(1/2,1,500);
% 
% A1 = [t1 ; 2*t1];
% A2 = [t2 ; -t2+3/4]; A2(:,[1,26])=[];
% A3 = [t3 ; 2*t3-3/4];
% 
% beta1 = [A1, A2, A3]+5;
% beta2 = [tvec/2+1.5; tvec-1/2]+5;

%%%%%%%%
%Pair 3

beta1 = [2*tvec-1; (2*tvec-1).*(2*tvec).*(2*tvec-2)+3];
beta2 = [2*ones(1,T); 2*tvec+2];


%% Calculate Geodesic. 
% -) Iter decides how often dynamic programming and optimization over starting
%    point are iterated
% -) optrigid decides if rigid motions are factored out
% -) doplot decides if the minimizing geodesic is plotted
tic
% [dis, v, qv, gamma_opt] =...
    H2_Main_Geodesic(beta1,beta2,'iter',1,'optrigid',0,'doplot',1);
toc

