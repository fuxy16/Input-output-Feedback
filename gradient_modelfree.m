function grad = gradient_modelfree(A,B,C,Q,R,K,p,r,max_grad)
% Calculate gradient with a simulator
M = 20; %100
dif = 0;
[m,n]=size(K);
for j = 1:M
    U = randn(m,n);
    U = U./norm(U,"fro");
    K1 = K + r*sqrt(m*n)*U;
    K2 = K - r*sqrt(m*n)*U;
    V1=cost_modelfree(A,B,C,K1,Q,R,p);
     V2=cost_modelfree(A,B,C,K2,Q,R,p);
    dif = dif + (V1-V2)*U;
end
grad = dif/(2*r*M);
if norm(grad)>max_grad
    grad=max_grad*grad/norm(grad);
end
end

