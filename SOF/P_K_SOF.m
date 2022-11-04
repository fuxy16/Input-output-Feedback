function p_K = P_K_SOF(A,B,C,Q,R,K)
p_K=dlyap((A-B*K*C)',Q+C'*K'*R*K*C);
end

