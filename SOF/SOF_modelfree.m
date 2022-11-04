function K = SOF_modelfree(A,B,C,Q,R,sigma_0,K_ini,stepsize,num_iter)
K=K_ini;

for i=1:num_iter
% %p_K
% p_K=P_K_SOF(A,B,C,Q,R,K);
% %e_K
% e_K=(R+B'*p_K*B)*K*C-B'*p_K*A;
% %sigma_K
% sigma_K=Sigma_K_SOF(K,A,B,C,sigma_0);
% 
% delta_K0=2*e_K*sigma_K*C'
delta_K=gradient_modelfree_SOF(A,B,C,Q,R,K,0.01,10000,sigma_0);

K=K-stepsize*delta_K;

end
end

