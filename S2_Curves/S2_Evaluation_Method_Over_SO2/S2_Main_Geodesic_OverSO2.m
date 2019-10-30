% Calculate geodesic between curves beta1 and beta2
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

function [dis,v,qv,gamma_opt] = S2_Main_Geodesic_OverSO2(beta1,beta2,varargin)
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

[dimG, T] = size(beta1);
n=zeros(1,dimG);n(1,dimG)=1;n=n';
alpha1=S2_lift(S2_startPoint(beta1,n), beta1);
alpha2=S2_lift(S2_startPoint(beta2,n), beta2);

%% get q 

q1 = S2_curveToq(alpha1);
q2 = S2_curveToq(alpha2);

%% Get Optimal (alpha2(0), q2)

dis = inf;
Threshold_Value = inf;
ite = 0;
gamma_opt = linspace(0,1,T);

while (Threshold_Value > 0.001 && ite < max_iter)
    
Dis = dis;
if optRigid == 0 %no Rigid Motions
    [alpha2(:,:,1),q2]=...
    S2_optimization_over_SO2(alpha1(:,:,1),q1,alpha2(:,:,1),q2); % optimization over SO(2)
else %optimizing over rigid motions
    alpha2(:,:,1) = alpha1(:,:,1); %align the starting points
    [q2]=S2_optimization_over_SO2_modG(q1,q2);
end    

[q2, gamma]=S2_get_q_gamma_Rep(q1,q2); % reparametrizationSn 
dis = S2_distance_between_Orbits(alpha1(:,:,1),q1,alpha2(:,:,1),q2);
Threshold_Value = norm(Dis-dis);

gamma_opt = composeGamma(gamma_opt,gamma);
ite = ite+1;
end

if nargout>1
 [v, qv] = S2_get_initial_velocity(alpha1(:,:,1),q1,alpha2(:,:,1),q2);
end 

if doplot==1
    figure;
    plot(linspace(0,1,T),gamma_opt);

    %% distance
    sprintf('The distance between the two curves is %0.3f',dis)
    
    %% Get the geodesic between alpha1 and alpha2
    S=100; %number of point in the variation direction
    alphab= S2_get_GeodesicInSOn(alpha1(:,:,1),q1,alpha2(:,:,1),q2,S);
    
    %% Project to S^2 and Plot interpolating geodesics
    
    betab = zeros(dimG, T+1, S+1);
    figure
    for i = 1: 20: S+1
        betab(:,:,i) = S2_projection(alphab(:,:,:,i));
        plot3(betab(1,:,i),betab(2,:,i),betab(3,:,i),'k','LineWidth',1)
        hold on
    end
    plotS2
    plot3(betab(1,1,1),betab(2,1,1),betab(3,1,1),'*r',...
        betab(1,1,S+1),betab(2,1,S+1),betab(3,1,S+1),'*r','LineWidth',1);
    plot3(betab(1,:,1),betab(2,:,1),betab(3,:,1),...
        betab(1,:,S+1),betab(2,:,S+1),betab(3,:,S+1),'Color','k','LineWidth',2);
    axis equal;
    axis off
    set(gcf,'color','w');
    hold off
end
