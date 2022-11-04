function sigma_k = Sigma_K(K,A,B,M_inv,sigma_0)
% calculate Sigma_K with model
sigma_k=dlyap((A-B*K*M_inv)',sigma_0);
end

