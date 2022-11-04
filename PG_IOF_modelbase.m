function [K,costs] = PG_IOF_modelbase(A,B,M_inv,Qc,R,sigma_0,K_ini,stepsize,num_iter,num_rec)
% policy gradient with model

K=K_ini;
costs=[];

for i=1:num_iter
    %p_K
    p_K=P_K(A,B,M_inv,Qc,R,K);
    %e_K
    e_K=(R+B'*p_K*B)*K*M_inv-B'*p_K*A;
    %sigma_K
    sigma_K=Sigma_K(K,A,B,M_inv,sigma_0);
    
    delta_K=2*e_K*sigma_K*M_inv';
    K=K-stepsize*delta_K;
    
    if mod(i,num_rec)==1
        p_K_PG_output = P_K(A,B,M_inv,Qc,R,K);
        costs=[costs,cost(p_K_PG_output,sigma_0)];
    end
end
end

