function p_K = P_K(A,B,M_inv,Q,R,K)
% calculate p_K 
p_K=dlyap((A-B*K*M_inv)',Q+M_inv'*K'*R*K*M_inv);
end

