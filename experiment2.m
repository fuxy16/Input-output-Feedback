% model free experiment
%% init
clear
close all

num=20; % experiment iter
rou_A=0.8;
state_dim=4;
input_dim=2;
output_dim=2;
stepsize=0.00001;
maxiter=10000;
num_rec=100;

A=randn(state_dim,state_dim);
A=A/max(abs(eig(A)))*rou_A; 
B=randn(state_dim,input_dim);
C=randn(output_dim,state_dim);
Q=eye(output_dim);
R=0.01*eye(input_dim);
load('data_experiment.mat')
sigma_0=20*(A*A)*(A*A)'+20*(A*B)*(A*B)'+20*B*B';

cost_PG_output=zeros(num,1);
cost_PG_output_rec=[];
time_PG_output=zeros(num,1);

%% optimal control with state feedback 1
Qc=C'*Q*C;
[K_state_feedback,S,e]=dlqr(A,B,Qc,R);
cost_optimal=trace(S*sigma_0);

%% pg_output
for iter=1:num
tic
p=calculate_p(A,B,C,state_dim);
z_dim=p*(input_dim+output_dim);

M=calculate_M(A,B,C,p);
M_inv=M'*inv(M*M');
K_ini = zeros(input_dim,z_dim);

[K_PG_output,costs]=PG_IOF_modelfree(A,B,C,M_inv,Q,R,p,K_ini,stepsize,maxiter,num_rec);
p_K_PG_output = P_K(A,B,M_inv,Qc,R,K_PG_output);
cost_PG_output(iter) = cost(p_K_PG_output,sigma_0);

time_PG_output(iter)=toc;
cost_PG_output_rec=[cost_PG_output_rec;costs];
end

% %% plot
% mean_cost=mean(cost_PG_output_rec,1);
% max_cost=max(cost_PG_output_rec);
% min_cost=min(cost_PG_output_rec);
% 
% plot(cost_optimal*ones(size(mean_cost)),'--','LineWidth',1.5,'color',[65 105 225]/255)
% hold on
% plot(mean_cost,'LineWidth',1.5,'color',[255 153 18]/255)
% hold on
% temp=1:length(mean_cost);
% h = fill([temp fliplr(temp)], [min_cost, fliplr(max_cost)], [255 227 132]/255);
% set(h,'edgealpha',0,'facealpha',0.7) 
% 
% set(gca,'FontSize',14);
% xlabel('Iteration','FontSize',14)
% ylabel('cost','FontSize',16,'Interpreter','latex','rotation',0)
% legend('optimal cost','pg output','Location','northeast')
% grid on
%% plot
mean_cost=mean(cost_PG_output_rec,1);
max_cost=max(cost_PG_output_rec);
min_cost=min(cost_PG_output_rec);

%set(gca,'FontSize',14,'YScale','log','XLim',[0,4900],'YLim',[5e-2,3],'YtickMode','manual','YTick',10.^(-2:0),'Box','on');
set(gca,'FontSize',14,'YScale','log','XLim',[0,9900],'YLim',[1e-2,5],'YtickMode','manual','YTick',10.^(-2:0),'Box','on');
xlabel('Iteration','FontSize',14)
ylabel('Optimality gap','FontSize',14)
set(gcf,'unit','centimeters','position',[1,2,14,8])
hold on

%temp=0:100:4900;
temp=0:100:9900;
plot(temp,mean_cost/cost_optimal-1,'LineWidth',1.5,'color',[255 153 18]/255)
hold on
h = fill([temp fliplr(temp)], [min_cost, fliplr(max_cost)]/cost_optimal-1, [255 227 132]/255);
set(h,'edgealpha',0,'facealpha',0.3) 