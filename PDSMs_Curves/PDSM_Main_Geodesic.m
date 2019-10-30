%Calculate geodesic between curves beta1 and beta2
% Input:
% beta1,beta2....Input Curves (same number of points is assumed)
% -) Iter decides how often dynamic programming and optimization over starting
%    point are iterated
% -) optrigid decides if rigid motions are factored out
% -) doplot decides if the minimizing geodesic is plotted
% Output: 
%    dist...minimal distance between orbits
%    shooting_vector...intial velocity
%    opt_gamma...optimal reparametrization

function [dis,v,qv,gamma_opt] = PDSM_Main_Geodesic(beta1, beta2,varargin)
max_iter = 1; %number of iteration, if not specified differently 
optRigid = 0; %no optimization over rigid motions, if not specified differently 
doplot = 0; %no plots, if not specified differently 
%% Read initial data
ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii},'char'))
        switch (lower(varargin{ii}))
            case 'optrigid'
                ii = ii + 1;
                optRigid = varargin{ii};  
            case 'iter'
                ii = ii + 1;
                max_iter = varargin{ii};
            case 'doplot'
                ii = ii + 1;
                doplot = varargin{ii};  
            otherwise
                error('Invalid option: ''%s''.',varargin{ii});
        end  
    end
    ii = ii + 1;
end

%% lift

[m, ~, T] = size(beta1);
alpha1=PDSM_lift(beta1);
alpha2=PDSM_lift(beta2);

%% get q 

q1=PDSM_to_q(alpha1);
q2=PDSM_to_q(alpha2);

%% Get Optimal (alpha2(0),q2)

epsilon=0.17;
dis = inf;
Threshold_Value = inf;
ite = 0;
gamma_opt = linspace(0,1,T);

while (Threshold_Value > 0.001 && ite < max_iter)
    
Dis = dis;
if optRigid == 0 %no Rigid Motions
    [alpha2(:,:,1),q2]=...
        PDSM_gradientMethod(alpha1(:,:,1),q1,alpha2(:,:,1),q2,epsilon); % Gradient Method
else %optimizing over rigid motions
    alpha2(:,:,1) = alpha1(:,:,1); %align the starting points
    [q2]=PDSM_gradientMethod_ModG(q1,q2,epsilon); 
end    
       
[q2, gamma]=PDSM_get_q_gammaRep(q1,q2); % reparametrizationH2
dis = PDSM_distance_between_Orbits(alpha1(:,:,1),q1,alpha2(:,:,1),q2);
Threshold_Value = norm(Dis-dis);

gamma_opt = composeGamma(gamma_opt,gamma);
ite = ite+1;
end

if nargout>1
 [v, qv] = PDSM_get_initial_velocity(alpha1(:,:,1),q1,alpha2(:,:,1),q2);
end 
    
if doplot==1
    figure;
    plot(linspace(0,1,T),gamma_opt);
    
    %% distance
    sprintf('The distance between the two curves is %0.3f',dis)
    
    %% get the geodesic between alpha1 and alpha2
    S = 10;
    alphab=PDSM_getGeodesicInSLn(alpha1(:,:,1),q1,alpha2(:,:,1),q2,S);
    
    %% project to PDSMs
    
    betab = zeros(m, m, T+1, S+1);
    for i = 1: T+1
        for k = 1: S+1
            betab(:,:,i,k) = sqrtm(alphab(:,:,i,k)*alphab(:,:,i,k)');
        end
    end
    
    % plot
    figure
    for i = 1: 20: T+1
        for k = 1: S+1
            C(:, k, i) = [i/2, 3*k, 0]';
            plotEllipse(betab(:,:,i,k),C(:,k,i)) %%%
            hold on
        end
    end
    axis equal
%     axis off
    hold off
    set(gcf,'color','w');
end


