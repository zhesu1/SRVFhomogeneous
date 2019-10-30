function g = S2_lift(g_0,f)

[m,T] = size(f);

    g(:,:,1)=g_0;

for i=2:T
    g(:,:,i)=(eye(m)-2*(f(:,i-1)+f(:,i))*(f(:,i-1)+f(:,i))'/(norm(f(:,i-1)...
        +f(:,i)))^2)*(eye(m)-2*f(:,i-1)*f(:,i-1)')*g(:,:,i-1);
end