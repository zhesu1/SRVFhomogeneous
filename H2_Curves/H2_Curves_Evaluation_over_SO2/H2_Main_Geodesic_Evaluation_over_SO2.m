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

function [dis, v, qv, gamma_opt]=...
    H2_Main_Geodesic_Evaluation_over_SO2(beta1,beta2,varargin)
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

[~, T] = size(beta1);
alpha1=H2_lift(beta1);
alpha2=H2_lift(beta2);

%% get q 
q1=H2_to_q(alpha1);
q2=H2_to_q(alpha2);

%% Get Optimal (alpha2(0),q2)

dis = inf;
Threshold_Value = inf;
ite = 0;
gamma_opt = linspace(0,1,T);

while (Threshold_Value > 0.001 && ite < max_iter)
    
Dis = dis;
if optRigid == 0 %no Rigid Motions
    [alpha2(:,:,1),q2]=...
        H2_optimization_over_SO2(alpha1(:,:,1),q1,alpha2(:,:,1),q2); % evaluation over SO2
else %optimizing over rigid motions
    alpha2(:,:,1) = alpha1(:,:,1); %align the starting points
    [q2] = H2_optimization_over_SO2_modG(q1,q2);
end    
       
[q2,gamma]=H2_get_q_gammaRep(q1,q2); % reparametrizationH2
dis = H2_distance_between_Orbits(alpha1(:,:,1),q1,alpha2(:,:,1),q2);
Threshold_Value = norm(Dis-dis);

gamma_opt = composeGamma(gamma_opt,gamma);
ite = ite+1;
end

if nargout>1
 [v, qv] = H2_get_initial_velocity(alpha1(:,:,1),q1,alpha2(:,:,1),q2);
end 
    
if doplot==1
    figure;
    plot(linspace(0,1,T),gamma_opt); 

    %% distance
    
    sprintf('The distance between the two curves is %0.3f',dis)
    
    %% get the geodesic between alpha1 and alpha2 and initial velocity (v, qv)
    
    S = 100; %number of point in the variation direction
    alphab = H2_get_GeodesicInSL2(alpha1(:,:,1),q1,alpha2(:,:,1),q2,S);
    
    %% project to H^2 and plot interpolating geodesics
   
    figure  
    x=zeros(T+1,S+1);
    y=zeros(T+1,S+1);
    for j=1:10:S+1
        x(:,j)=(alphab(1,2,:,j).*alphab(2,2,:,j)+alphab(1,1,:,j).*alphab(2,1,:,j))...
            ./(alphab(2,2,:,j).^2+alphab(2,1,:,j).^2);
        y(:,j)=(alphab(1,1,:,j).*alphab(2,2,:,j)-alphab(1,2,:,j).*alphab(2,1,:,j))...
            ./(alphab(2,2,:,j).^2+alphab(2,1,:,j).^2);
        plot(x(:,j),y(:,j),'Color','k','LineWidth',1)
        hold on
    end
    plot(x(:,1),y(:,1),x(:,S+1),y(:,S+1),'Color','k','LineWidth',2)
    plot(x(1,1),y(1,1),'*r',x(1,S+1),y(1,S+1),'*r','LineWidth',1)
%     plot(beta1(1,:), beta1(2,:), beta2(1,:), beta2(2,:),'Color','r','LineWidth',2)


%     x1=zeros(T+1,S+1);
%     y2=zeros(T+1,S+1);
%     for i=1:10:T+1
%         x1(i,:)=(alphab(1,2,i,:).*alphab(2,2,i,:)+alphab(1,1,i,:).*alphab(2,1,i,:))...
%             ./(alphab(2,2,i,:).^2+alphab(2,1,i,:).^2);
%         y2(i,:)=(alphab(1,1,i,:).*alphab(2,2,i,:)-alphab(1,2,i,:).*alphab(2,1,i,:))...
%             ./(alphab(2,2,i,:).^2+alphab(2,1,i,:).^2);
%         plot(x1(i,:),y2(i,:),'k--')
%     end
%     axis equal
    set(gcf,'color','w');
    box off
    hold off
end
