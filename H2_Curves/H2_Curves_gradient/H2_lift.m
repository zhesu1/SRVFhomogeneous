function y = H2_lift(x)

% lift
[~,T]=size(x);

%% first step

args=2*pi*rand;
y_args=[cos(args) -sin(args);sin(args) cos(args)];

yy = zeros(2,2,T);
for i=1:T
    yy(:,:,i)=[sqrt(x(2,i)), x(1,i)/sqrt(x(2,i));...
         0, 1/sqrt(x(2,i))]*y_args;
end

%% second step

y(:,:,1)=yy(:,:,1);


for i=2:T
    
    y(:,:,i)=y(:,:,i-1)*sqrtm(inv(y(:,:,i-1))*yy(:,:,i)...
        *yy(:,:,i)'*inv(y(:,:,i-1)'));
end