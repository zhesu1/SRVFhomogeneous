function [f2new,q2new] = H2_optimization_over_SO2(f1,q1,f2,q2)

m = 100;

[~,~,T]=size(q1);
% thetaa=linspace(0,2*pi,m+1);

theta = zeros(1, m+1);
y_theta = zeros(2,2,m+1);

for k = 1: m+1
    theta(k) = 2*pi*(k-1)/m;
    y_theta(:,:,k) = [cos(theta(k)) -sin(theta(k));...
                      sin(theta(k)) cos(theta(k))];
end

S_int = q1(1,1,:).*q2(1,1,:)+q1(1,2,:).*q2(1,2,:);
T_int = q1(1,1,:).*q2(1,2,:)-q1(1,2,:).*q2(1,1,:);
SS=trapz(linspace(0,1,T),S_int,3);
TT=trapz(linspace(0,1,T),T_int,3);

Norm__square1 = zeros(1, m+1);
dis = zeros(1,m+1);
for k = 1: m+1
    Norm__square1(k)=norm(H2_invRieExpOnSLn(f1\f2*y_theta(:,:,k)),'fro')^2;
    dis(k) = Norm__square1(k)-4*(SS*cos(4*pi*(k-1)/m)+TT*sin(4*pi*(k-1)/m));
end

[min_value,ind]=min(dis);
% theta_opt=2*pi*(ind-1)/m;
% 
% figure;
% plot(thetaa,dis,theta_opt,min_value,'*');
% xlim([0, 2*pi])
% ylabel('d_s');


f2new = f2*y_theta(:,:,ind);
q2new = zeros(size(q2));
for j=1:T
     q2new(:,:,j)=y_theta(:,:,ind)\q2(:,:,j)*y_theta(:,:,ind);
 end


