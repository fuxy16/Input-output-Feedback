function K = SOF_modelbase(A,B,C,Q,R,sigma_0,K_ini,stepsize,num_iter)
K=K_ini;

for i=1:num_iter
%p_K
p_K=P_K_SOF(A,B,C,Q,R,K);
%e_K
e_K=(R+B'*p_K*B)*K*C-B'*p_K*A;
%sigma_K
sigma_K=Sigma_K_SOF(K,A,B,C,sigma_0);

delta_K=2*e_K*sigma_K*C';
K=K-stepsize*delta_K;

end
end

