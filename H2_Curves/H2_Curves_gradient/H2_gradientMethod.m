function [f2new,q2new] = H2_gradientMethod(f1,q1,f2,q2,epsilon)

% Gradient method

[~,~,T]=size(q1);
tvec=linspace(0,1,T+1);
Gnorm=inf;
ite = 0;

while (Gnorm>1e-5 && ite < 25)
    
    A=f2\f1;
    
    inside = zeros(size(q1));
    for i=1:T
        inside(:,:,i)=q1(:,:,i)*q2(:,:,i)'-q2(:,:,i)'*q1(:,:,i);
    end
    
    GradientF = 2*(-H2_invRieExpOnSLn(A)+trapz(tvec(1:T),inside,3));
    GradientF = (GradientF-GradientF')/2;
   
    Gnorm = norm(GradientF,'fro');
    
    f2new=f2*H2_RieExpOnSLn(-epsilon*GradientF);

    q2new = zeros(size(q2));
    for j=1:T
       q2new(:,:,j)=...
          inv(H2_RieExpOnSLn(-epsilon*GradientF))*q2(:,:,j)*H2_RieExpOnSLn(-epsilon*GradientF);
    end

    f2 = f2new;
    q2 = q2new;
    ite = ite + 1;
end