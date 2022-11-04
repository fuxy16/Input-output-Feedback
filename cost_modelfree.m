function J = cost_modelfree(A,B,C,K,Q,R,p)
% Calculate cost with a simulator
% Monte Carlo
%tic
num=1; %num of Monte Carlo
step=20; 
[state_dim,input_dim]=size(B);
[output_dim,state_dim]=size(C);

% sigma=0;
%     M=calculate_M(A,B,C,p);
%     M_inv=M'*inv(M*M');
%     p_K = P_K(A,B,M_inv,C'*Q*C,R,K);

u=zeros(input_dim,p+step);
x=zeros(state_dim,p+step);
y=zeros(output_dim,p+step);

J=zeros(num,1);
for i=1:num
    u(:,1)=mvnrnd(zeros(input_dim,1),20*eye(input_dim))';
    x(:,1)=mvnrnd(zeros(state_dim,1),20*eye(state_dim))';
    y(:,1)=C*x(:,1);
    %initial p steps without input
    for t=2:p
        u(:,t)=mvnrnd(zeros(input_dim,1),eye(input_dim))';
        x(:,t)=A*x(:,t-1)+B*u(:,t-1);
        y(:,t)=C*x(:,t);
    end
    
    % last steps
    for t=(p+1):(p+step)
        u_tp=vec(u(:,(t-1):-1:(t-p)));
        y_tp=vec(y(:,(t-1):-1:(t-p)));
        z_t=[u_tp;y_tp];
        u(:,t)=-K*z_t;
        x(:,t)=A*x(:,t-1)+B*u(:,t-1);
        y(:,t)=C*x(:,t);
        J(i)=J(i)+y(:,t)'*Q*y(:,t)+u(:,t)'*R*u(:,t);
    end

%     [p_K;C'*Q*C+M_inv'*K'*R*K*M_inv+(A-B*K*M_inv)'*p_K*(A-B*K*M_inv)];
%     J(i);
%     x(:,(p+1))'*p_K*x(:,(p+1));
    
%     sigma=sigma+x(:,(p+1))*x(:,(p+1))';
end
% sigma_00=A*A*A'*A'+A*B*B'*A'+B*B';
% sigma=sigma/num;
J=mean(J);
%t=toc
% cost(p_K,sigma_00);
end

