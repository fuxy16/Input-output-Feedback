% model base experiment
%% init
clear
close all

num=20; % experiment iter num
rou_A=0.8; %
state_dim=4;
input_dim=2;
output_dim=2;
stepsize=0.001;%
num_iter=5000;
num_rec=100;

A=randn(state_dim,state_dim);
A=A/max(abs(eig(A)))*rou_A;
B=randn(state_dim,input_dim);
C=randn(output_dim,state_dim);
Q=eye(output_dim);
R=0.01*eye(input_dim);
sigma_0=eye(state_dim);
load('data_experiment.mat')

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
    K_ini=0.1*randn(input_dim,z_dim);
    while max(abs(eig(A-B*K_ini*M_inv)))>0.8
        K_ini=0.1*randn(input_dim,z_dim);
    end

[K_PG_output,costs]=PG_IOF_modelbase(A,B,M_inv,Qc,R,sigma_0,K_ini,stepsize,num_iter,num_rec);
p_K_PG_output = P_K(A,B,M_inv,Qc,R,K_PG_output);
cost_PG_output(iter) = cost(p_K_PG_output,sigma_0);

time_PG_output(iter)=toc;
cost_PG_output_rec=[cost_PG_output_rec;costs];

end

%% plot
mean_cost=mean(cost_PG_output_rec,1);
max_cost=max(cost_PG_output_rec);
min_cost=min(cost_PG_output_rec);

set(gca,'FontSize',14,'YScale','log','XLim',[0,4900],'YLim',[1e-6,1e1],'YtickMode','manual','YTick',10.^(-6:2:1),'Box','on');
%set(gca,'FontSize',14,'YScale','log','XLim',[0,4900],'YLim',[5e-2,3],'YtickMode','manual','YTick',10.^(-1:0),'Box','on');
xlabel('Iteration','FontSize',14)
ylabel('Optimality gap','FontSize',14)
set(gcf,'unit','centimeters','position',[1,2,14,8])
hold on

temp=0:100:4900;
plot(temp,mean_cost/cost_optimal-1,'LineWidth',1.5,'color',[255 153 18]/255)
hold on
h = fill([temp fliplr(temp)], [min_cost, fliplr(max_cost)]/cost_optimal-1, [255 227 132]/255);
set(h,'edgealpha',0,'facealpha',0.3) 
