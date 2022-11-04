function grad = gradient_modelfree_SOF(A,B,C,Q,R,K,r,max_grad,sigma_0)
M = 20; %100
dif = 0;
[m,n]=size(K);

for j = 1:M
    U = randn(m,n);
    U = U./norm(U,"fro");
    K1 = K + r*sqrt(m*n)*U;
    K2 = K - r*sqrt(m*n)*U;
%     p_K_SOF1 = P_K_SOF(A,B,C,Q,R,K1);
%     p_K_SOF2 = P_K_SOF(A,B,C,Q,R,K2);
% V1 = trace(p_K_SOF1*sigma_0);
% V2 = trace(p_K_SOF2*sigma_0);
    V1=cost_modelfree_SOF(A,B,C,Q,R,K1,sigma_0);
     V2=cost_modelfree_SOF(A,B,C,Q,R,K2,sigma_0);
    dif = dif + (V1-V2)*U;
end
grad = dif/(2*r*M);
if norm(grad)>max_grad
    grad=max_grad*grad/norm(grad);
end
end


