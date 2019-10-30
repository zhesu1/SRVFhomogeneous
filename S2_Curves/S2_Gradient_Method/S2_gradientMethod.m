function [f2new,q2new] = S2_gradientMethod(f1,q1,f2,q2,epsilon)

% Gradient method

[dim,~,T]=size(q1);
tvec=linspace(0,1,T);
Gnorm=inf;
ite = 0;

while (Gnorm>1e-5 && ite<25) 
    
    A=f2\f1;
    logmA= logm(A);
    
    inside = zeros(dim,dim, T);
    for i=1:T
        inside(:,:,i)=q1(:,:,i)*q2(:,:,i)'-q2(:,:,i)'*q1(:,:,i);
    end
    
    GradientF=2*(-logmA+trapz(tvec,inside,3));
    GradientF=1/2*(GradientF-GradientF'); %project to so(n)
    GradientF(1:(dim-1),dim)=0;
    GradientF(dim,1:(dim-1))=0;
    Gnorm = norm(GradientF,'fro');
 

    f2new=f2*expm(-epsilon*GradientF);
    
    
    q2new = zeros(size(q2));
    for j=1:T
        q2new(:,:,j)=expm(epsilon*GradientF)*q2(:,:,j)*expm(-epsilon*GradientF);
    end
  
    f2=f2new;
    q2=q2new;
    ite = ite+1; 
end
end
