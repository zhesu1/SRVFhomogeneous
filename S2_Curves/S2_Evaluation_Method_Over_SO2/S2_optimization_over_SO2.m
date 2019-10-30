function [f2n,q2n] = S2_optimization_over_SO2(f1,q1,f2,q2)

m=200;
[~,~,T]=size(q1);
tvec=linspace(0,1,T);

 A=f1'*f2;
 
 S_int=q1(1,3,:).*q2(1,3,:)+q1(2,3,:).*q2(2,3,:);
 T_int=q1(1,3,:).*q2(2,3,:)-q1(2,3,:).*q2(1,3,:);
 SS=trapz(tvec(1:T),S_int,3);
 TT=trapz(tvec(1:T),T_int,3);
 
 d_s = zeros(1,m+1);
for i=1:(m+1)
    d_s(i)=2*(acos(((A(1,1)+A(2,2))*cos(2*pi*(i-1)/m)+(A(1,2)-A(2,1)) ...
        *sin(2*pi*(i-1)/m)+A(3,3)-1)/2))^2-4*(SS*cos(2*pi*(i-1)/m) ...
        +TT*sin(2*pi*(i-1)/m));
end

% thetaa=linspace(0,2*pi,m+1);

[mindis,ii]=min(d_s);
args=2*pi*(ii-1)/m;

% figure;
% plot(thetaa,d_s,args,mindis,'*');
% xlim([0, 2*pi])
% ylabel('d_s');
 
 y_args=[cos(args) -sin(args) 0;sin(args) cos(args) 0;0 0 1];
 yinv_args=[cos(args) sin(args) 0;-sin(args) cos(args) 0;0 0 1];
 
 f2n=f2*y_args;
 q2n = zeros(size(q2));
 for j=1:T
     q2n(:,:,j)=yinv_args*q2(:,:,j)*y_args;
 end

