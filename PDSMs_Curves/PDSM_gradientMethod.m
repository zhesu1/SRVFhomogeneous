function [f2new,q2new] = PDSM_gradientMethod(f1,q1,f2,q2,epsilon)

% Gradient method

[m,~,T]=size(q1);
tvec=linspace(0,1,T+1);
Gnorm=inf;
ite = 0;

while (Gnorm>1e-5 && ite < 30)
    
    A=f2\f1;
    
    inside = zeros(m, m, T);
    for i=1:T
        inside(:,:,i)=q1(:,:,i)*q2(:,:,i)'-q2(:,:,i)'*q1(:,:,i);
    end
    
    GradientF=2*(-invRieExpOnSLn(A)+trapz(tvec(1:T),inside,3));
    GradientF=(GradientF-GradientF')/2;
    Gnorm = norm(GradientF,'fro');

    f2new=f2*RieExpOnSLn(-epsilon*GradientF);

    q2new = zeros(m, m, T);
    for j=1:T
        q2new(:,:,j)=inv(RieExpOnSLn(-epsilon*GradientF))*q2(:,:,j)*RieExpOnSLn(-epsilon*GradientF);
    end
    f2=f2new;
    q2=q2new;
    ite = ite+1;
end
end