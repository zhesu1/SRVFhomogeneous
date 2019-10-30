function [q2n] = S2_optimization_over_SO2_modG(q1,q2)


m=200;
[~,~,T]=size(q1);
tvec=linspace(0,1,T);

Sint=q1(1,3,:).*q2(1,3,:)+q1(2,3,:).*q2(2,3,:);
Tint=q1(1,3,:).*q2(2,3,:)-q1(2,3,:).*q2(1,3,:);
SS=trapz(tvec(1:T),Sint,3);
TT=trapz(tvec(1:T),Tint,3);
 
d_s = zeros(1, m+1);
for i=1:(m+1)
    d_s(i)=-4*(SS*cos(2*pi*(i-1)/m)+TT*sin(2*pi*(i-1)/m));
end

[mindis,ii]=min(d_s);
args=2*pi*(ii-1)/m;
% figure;
% plot(thetaa,d_s,args,mindis,'*');
% xlim([0, 2*pi])
% ylabel('d_s');
 
 y_args=[cos(args) -sin(args) 0;sin(args) cos(args) 0;0 0 1];
 yinv_args=[cos(args) sin(args) 0;-sin(args) cos(args) 0;0 0 1];
 
 q2n = zeros(size(q2));
 for j=1:T
     q2n(:,:,j)=yinv_args*q2(:,:,j)*y_args;
 end
