% 不同方法的简单对比
%% init
clear
addpath(genpath('./SOF'))

num=10; % experiment iter
rou_A=0.8;
state_dim=4;
input_dim=2;
output_dim=4;
num_rec=100;

cost_optimal=zeros(num,1);
cost_PG_output=zeros(num,1);
cost_SOF=zeros(num,1);

for iter=1:num
A=randn(state_dim,state_dim);
A=A/max(abs(eig(A)))*rou_A;
B=randn(state_dim,input_dim);
C=randn(output_dim,state_dim);

Q=eye(output_dim);
R=0.01*eye(input_dim);
sigma_0=eye(state_dim);

%% optimal control with state feedback 1
[K_state_feedback,S,e]=dlqr(A,B,C'*Q*C,R);
cost_optimal(iter)=trace(S*sigma_0);

%% pg_output
stepsize=0.00001;%0.00001
num_iter=1e5;%5e4
Qc=C'*Q*C;

p=calculate_p(A,B,C,state_dim);
z_dim=p*(input_dim+output_dim);
M=calculate_M(A,B,C,p);
M_inv=M'*inv(M*M');
K_ini = zeros(input_dim,z_dim);

[K_PG_output,costs]=PG_IOF_modelbase(A,B,M_inv,Qc,R,sigma_0,K_ini,stepsize,num_iter,num_rec);
if max(abs(eig(A-B*K_PG_output*M_inv)))>1 
    "pg output wrong"
end
p_K_PG_output = P_K(A,B,M_inv,Qc,R,K_PG_output);
cost_PG_output(iter) = cost(p_K_PG_output,sigma_0);

%% SOF
stepsize=0.00005;
num_iter=1e5;
K_ini=zeros(input_dim,output_dim);

K_SOF=SOF_modelbase(A,B,C,C'*Q*C,R,sigma_0,K_ini,stepsize,num_iter);
p_K_SOF = P_K_SOF(A,B,C,C'*Q*C,R,K_SOF);
if max(abs(eig(A-B*K_SOF*C)))>1 
   "SOF wrong"
end
cost_SOF(iter) = trace(p_K_SOF*sigma_0);

end

c_optimal=mean(cost_optimal)
c_PG_output=mean(cost_PG_output)
c_SOF=mean(cost_SOF)