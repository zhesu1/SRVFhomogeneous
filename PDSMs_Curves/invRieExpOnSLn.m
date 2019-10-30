function v=invRieExpOnSLn(y)

%---------------------------------------------------
% using the polar decomposition Y=PV to give the responding symmetric
% matrix P=sqrtm(Y*Y')
z=sqrtm(y*y');
z1=z;
Norm=inf;
while (Norm>1e-4)
    z1=vfzTovzInvexpOnSLn(z1,y);
    f1=z1*expm(logm(z1')-logm(z1));
    vc=logm(f1\y);
    Norm=norm(vc,'fro');
    znew=z1;
end
v=logm(znew);
v=v';

