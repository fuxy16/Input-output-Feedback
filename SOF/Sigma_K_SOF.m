function sigma_k = Sigma_K_SOF(K,A,B,C,sigma_0)
sigma_k=dlyap((A-B*K*C)',sigma_0);
end

