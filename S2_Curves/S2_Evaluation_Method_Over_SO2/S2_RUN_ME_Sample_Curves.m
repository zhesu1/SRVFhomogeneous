clear;
T=100; %spatial discretization of the curve

tvec=linspace(0,1,T);

angle=pi/12;
matrix=[cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1];

beta1=[cos(tvec); zeros(1,T); sin(tvec)];
beta2=matrix*[cos(tvec); sin(tvec); zeros(1,T)];

%% Calculate Geodesic. 
% -) Iter decides how often dynamic programming and optimization over starting
%    point are iterated
% -) optrigid decides if rigid motions are factored out
% -) doplot decides if the minimizing geodesic is plotted

tic
% [dis, v, qv, gamma_opt] =...
    S2_Main_Geodesic_OverSO2(beta1, beta2,'iter',1,'optrigid',0,'doplot',1);
toc


