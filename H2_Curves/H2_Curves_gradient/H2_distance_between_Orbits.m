function d=H2_distance_between_Orbits(f1,q1,f2,q2)

[~,~,T]=size(q1);

y1=norm(H2_invRieExpOnSLn(inv(f1)*f2),'fro')^2;
inside = zeros(size(q1));
for i=1:T
    inside(:,:,i)=(q1(:,:,i)-q2(:,:,i))*(q1(:,:,i)-q2(:,:,i))';
end
y2=trace(trapz(linspace(0,1,T),inside,3));

d=sqrt(y1+y2);

