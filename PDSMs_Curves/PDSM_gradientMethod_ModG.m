function [q2new] = PDSM_gradientMethod_ModG(q1,q2,epsilon)

% Gradient method

[~,~,T]=size(q1);
tvec=linspace(0,1,T+1);
Gnorm=inf;
ite = 0;

while (Gnorm>1e-5 && ite <25)
    
    inside = zeros(size(q1));
    for i=1:T
        inside(:,:,i)=q1(:,:,i)*q2(:,:,i)'-q2(:,:,i)'*q1(:,:,i);
    end
    
    GradientF=2*(trapz(tvec(1:T),inside,3));
    GradientF=(GradientF-GradientF')/2;
    Gnorm = norm(GradientF,'fro');
    
    q2new = zeros(size(q1));
    for j=1:T
        q2new(:,:,j)=inv(RieExpOnSLn(-epsilon*GradientF))*q2(:,:,j)*RieExpOnSLn(-epsilon*GradientF);
    end
    q2=q2new;
    ite = ite+1;
end
end