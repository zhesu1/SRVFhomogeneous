function gamma_new = composeGamma(gamma2,gamma1)

[~, T] = size(gamma1);

gamma_new = interp1(linspace(0,1,T) , gamma2 ,gamma1,'linear');


