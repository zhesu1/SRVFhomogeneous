function beta_new = Resample_Curves_PDSM(beta, N)

[m, ~, T] = size(beta);


del = zeros(1, T);
del(1) = 0;
for r = 2:T
    del(r) = norm(beta(:,:,r) - beta(:,:,r-1),'fro');
end
cumdel = cumsum(del)/sum(del);

newdel = linspace(0, 1, N);

beta_new = zeros(m, m, N);
for i = 1: m
    for j = 1: m
        beta_new(i,j,:) = interp1(cumdel,squeeze(beta(i,j,:)),newdel,'linear');
    end
end
        
            