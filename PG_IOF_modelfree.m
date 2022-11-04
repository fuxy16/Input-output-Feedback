function [K,costs] = PG_IOF_modelfree(A,B,C,M_inv,Q,R,p,K_ini,stepsize,num_iter,num_rec)
% policy gradient without model
K=K_ini;
sigma_0=20*(A*A)*(A*A)'+20*(A*B)*(A*B)'+20*B*B';
costs=[];

for i=1:num_iter
    delta_K=gradient_modelfree(A,B,C,Q,R,K,p,0.1,1000);%r=0.2
    K=K-stepsize*delta_K;
    if mod(i,num_rec)==1
        p_K_PG_output = P_K(A,B,M_inv,C'*Q*C,R,K);
        costs=[costs,cost(p_K_PG_output,sigma_0)];
    end
end
end

