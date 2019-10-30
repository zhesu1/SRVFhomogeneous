clear;
T=100; %spatial discretization of the curve

tvec=linspace(0,1,T);

beta1=[cos(2*sin(tvec)+cos(tvec)).*cos(sin(tvec)+2*cos(tvec)); cos(2*sin(tvec)+cos(tvec)).*sin(sin(tvec)+2*cos(tvec)); sin(2*sin(tvec)+cos(tvec))];
beta2=[cos(sin(tvec)).*cos(tvec.^2); cos(sin(tvec)).*sin(tvec.^2); sin(sin(tvec))];

%% Calculate Geodesic. 
% -) Iter decides how often dynamic programming and optimization over starting
%    point are iterated
% -) optrigid decides if rigid motions are factored out
% -) doplot decides if the minimizing geodesic is plotted
tic
% [dis, v, qv, gamma_opt] =...
    S2_Main_Geodesic_gradient(beta1, beta2,'iter',10,'optrigid',0,'doplot',1);
toc


