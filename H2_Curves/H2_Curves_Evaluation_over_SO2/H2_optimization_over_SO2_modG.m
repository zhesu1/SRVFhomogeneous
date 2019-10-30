function [q2new] = H2_optimization_over_SO2_modG(q1,q2)

m = 100;

[~,~,T]=size(q1);
% thetaa=linspace(0,2*pi,m+1);

S_int = q1(1,1,:).*q2(1,1,:)+q1(1,2,:).*q2(1,2,:);
T_int = q1(1,1,:).*q2(1,2,:)-q1(1,2,:).*q2(1,1,:);
SS=trapz(linspace(0,1,T),S_int,3);
TT=trapz(linspace(0,1,T),T_int,3);

dis = zeros(1,m+1);
for k = 1: m+1
    dis(k) = -4*(SS*cos(4*pi*(k-1)/m)+TT*sin(4*pi*(k-1)/m));
end
[min_value,ind]=min(dis);
theta_opt=2*pi*(ind-1)/m;

% figure;
% plot(thetaa,dis,theta_opt,min_value,'*');
% xlim([0, 2*pi])
% ylabel('d_s');

y_theta = [cos(theta_opt) -sin(theta_opt);...
                      sin(theta_opt) cos(theta_opt)];

q2new = zeros(size(q2));
for j=1:T
     q2new(:,:,j)=y_theta\q2(:,:,j)*y_theta;
 end