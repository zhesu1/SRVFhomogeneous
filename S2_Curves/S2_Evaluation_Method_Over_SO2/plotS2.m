function plotS2

[theta,rho]=meshgrid(0:0.1:pi,0:0.1:2*pi); 
RRR = 1;
xxx = RRR*sin(theta).*cos(rho); 
yyy =RRR*sin(theta).*sin(rho); 
zzz = RRR*cos(theta); 
plot3(xxx,yyy,zzz,'Color',[0.5,0.5,0.5])
