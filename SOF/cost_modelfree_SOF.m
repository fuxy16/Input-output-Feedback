function J = cost_modelfree_SOF(A,B,C,Q,R,K,sigma_0)
%tic
num=5; %num of Monte Carlo
step=20; 
[state_dim,input_dim]=size(B);
[output_dim,state_dim]=size(C);

% sigma=0;
%     M=calculate_M(A,B,C,p);
%     M_inv=M'*inv(M*M');
%     p_K = P_K(A,B,M_inv,C'*Q*C,R,K);

u=zeros(input_dim,step);
x=zeros(state_dim,step);
y=zeros(output_dim,step);

J=zeros(num,1);
for i=1:num
    x(:,1)=mvnrnd(zeros(state_dim,1),sigma_0)';
    y(:,1)=C*x(:,1);
    u(:,1)=-K*y(:,1);
    J(i)=J(i)+x(:,1)'*Q*x(:,1)+u(:,1)'*R*u(:,1);

    for t=2:step
        x(:,t)=A*x(:,t-1)+B*u(:,t-1);
        y(:,t)=C*x(:,t);
        u(:,t)=-K*y(:,t);
        J(i)=J(i)+x(:,t)'*Q*x(:,t)+u(:,t)'*R*u(:,t);
    end

end

J=mean(J);
end

