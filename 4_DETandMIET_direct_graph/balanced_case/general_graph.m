close all
clear
clc
tic

%% initial matrices and parameters;
A = [0.1 1
      0  0];
B_u = [-1;1];
ADJ  = [0 0 0 0 1
        1 0 0 0 0
        0 1 0 0 0
        1 0 1 0 0
        0 0 0 1 0];
L = diag(sum(ADJ,2))-ADJ;

H = diag([2 1 1 1 2])/7;
h = [2 1 1 1 2]'/7;
phi = 0.618;
psi = norm(H^0.5*L*H^-0.5);

omega = 0.4;
c = 1/phi/(1-omega);
xishu = psi^2*(1+psi/phi)^2;

N = length(L);

nx = size(A,1);
nu = size(B_u,2);

% solov P and K;
Q = eye(nx);
P = are(A, B_u*B_u', Q); % are(A,B,C): A'*X + X*A - X*B*X + C = 0
K = -B_u'*P;


% event-triggered sigma setup
dt = 0.0001;
Ts = 30;
sigma = 0.0001*[1 1 1 1 1];
sigma_all = zeros(1+Ts/dt,N); % sigma information of all time
dsigma_dt = zeros(1,N); % the dot{sigma}
beta = 50;
delta = 0.01*[1 1 1 1 1];
Time = []; % triggering time
e = zeros(nx,N); % error information

%solve T_i
lambda = eig(K'*K/P);
lambda_k = max(lambda);
lambda = eig(Q/P);
lambda_q = min(lambda(lambda>10^-10));
for i = 1:(N-1)
    Deta(i) = 4*lambda_k^2*xishu*c^2/omega-(delta(i)+lambda_k-lambda_q)^2;
    L_T(i) = 2/sqrt(Deta(i))*(atan((2*beta+delta(i)+lambda_k-lambda_q)/sqrt(Deta(i)))-atan((delta(i)+lambda_k-lambda_q)/sqrt(Deta(i))));
end

for i = N:N
    Deta(i) = 4*lambda_k^2*xishu*c^2/omega-(delta(i)+lambda_k-lambda_q)^2
    L_T(i) = 2/sqrt(Deta(i))*(atan((2*beta+delta(i)+lambda_k-lambda_q)/sqrt(Deta(i)))-atan((delta(i)+lambda_k-lambda_q)/sqrt(Deta(i))))
end

% relevant matrix setup: the information of agent i at time t, and the
% information of all time
x = [];
x_all = zeros(1+Ts/dt,nx*N);
hatx = [];
hatx_all = zeros(1+Ts/dt,nx*N);

% 作图看具体的触发过程
ee_all = zeros(1+Ts/dt,nx*N);
eother_all = zeros(1+Ts/dt,nx*N);

% initial conditions of matrix
% x0 = 2*rand([1,nx*N])-1;
% xx = x0;
% save('x0_general','xx');
load x0_general.mat
x0 = xx*1;
hatx0 = x0;

%% simulation
% some calculation for control signals to update
syms s
% dt = 0.001;
% Ts = 10;
G = expm(A*dt);
H = int(expm(A*s),0,dt)*B_u;
H = eval(H);
count = 0;

% simulation for time
for t = 0:dt:Ts
    % The aim of matrix_i is to get information of every agent at every
    % instant, i.e., matrix_i is a tool, to get the whole information
    count = count+1;
    matrix_x = reshape(x0, nx, N);
    matrix_hat_x = reshape(hatx0, nx, N);

    % simulation for agents
    for i = 1:N
        if i == 1
            for j = 1:N
                matrix_hat_x(:,j) = G*matrix_hat_x(:,j); % update the estimator
            end
        end

        u = 0;
        delta_x = 0;
        for j = 1:N
            u = u+c*K*ADJ(i,j)*(matrix_hat_x(:,i)-matrix_hat_x(:,j));
            delta_x = delta_x+ADJ(i,j)*K*(matrix_hat_x(:,i)-matrix_hat_x(:,j));
        end
        matrix_x(:,i) = G*matrix_x(:,i)+H*u; %  update x;
        e(:,i) = matrix_hat_x(:,i)-matrix_x(:,i); % update error information

        fi = norm(K*e(:,i))^2-omega/xishu*norm(delta_x)^2-beta*sigma(i);

        ee_all(count,i) = norm(K*e(:,i))^2;
        eother_all(count,i) = omega/xishu*norm(delta_x)^2+beta*sigma(i);

        if fi >= 0 % event-triggered
            Time = [Time;t,i]; % triggering time, Time = [time1, agent i;time2, agent j]
            matrix_hat_x(:,i) = matrix_x(:,i); % update the information exchanged
        else
            dsigma_dt(i) = -delta(i)*sigma(i)+omega/xishu*norm(delta_x)^2-norm(K*e(:,i))^2;
            sigma(i) = sigma(i)+dsigma_dt(i)*dt;
        end
        sigma_all(count,i) = sigma(i);
    end

    % set new initial conditions for a new loop
    x0 = reshape(matrix_x, 1, nx*N);
    hatx0 = reshape(matrix_hat_x, 1, nx*N);

    % get the whole information of the system
    x_all(count,:) = x0;
    hatx_all(count,:) = hatx0;
end

%% plot figures
close all
set(0,'defaultfigurecolor','w')
color = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    0.87      0.45      0
    %0.07      0.62      1
    %1         0.41      0.16
    ];
color = kron(color,[1;1]);

% state x
figure(1)
for j = 1:nx
    subplot(nx,1,j)
    plot(0:dt:Ts,x_all(:,j:2:end),'-','linewidth',1.5)
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('${x}_{i',num2str(j),'}$'),'Interpreter','Latex');
    box on
end
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','FontName','Times','FontSize',8,'NumColumns',3,'FontWeight','bold','Location','southeast')

% % state hatx
% figure(2)
% for j = 1:nx
%     subplot(nx,1,j)
%     plot(0:dt:Ts,hatx_all(:,j:nx:end),'-','linewidth',1.5)
%     xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\hat{x}_{i',num2str(j),'}$'),'Interpreter','Latex');
%     box on
% end
% legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','FontName','Times','FontSize',8,'NumColumns',3,'FontWeight','bold','Location','southeast')
% 
% % trigger time
% % trigger记录了每个智能体触发的时间点，count记录了每个智能体触发的次数
% % interval记录了每个智能体的触发时间间隔
% figure(3)
% Trigger = zeros(1,N);
% count = ones(1,N);
% for i = 1:length(Time)
%     for j = 1:N
%         if Time(i,2) == j
%             Trigger(count(j),j) = Time(i,1);
%             count(j) = count(j)+1;
%         end
%     end
% end
% for i = 1:N
%     count(i) = count(i)-1;
% end
% interval = zeros(1,N);
% for i = 1:N
%     interval(1,i) = Trigger(1,i);
%     for j = 2:count(i)
%         interval(j,i) = Trigger(j,i)-Trigger(j-1,i);
%     end
% end
% 
% Min = 10*ones(1,N);
% for i = 1:N
%     for j = 1:count(i)
%         if interval(j,i) <= Min(i)
%             Min(i) = interval(j,i);
%         end
%     end
% end
% 
% for i = 1:N
%     subplot(N,1,i)
%     hold on
%     Mini = Min(i)*ones(1,Ts/dt+1);
%     plot(0:dt:Ts,Mini,'linestyle',':','color',[1 0 0],'LineWidth',1,'MarkerSize',2)% itself
%     plot(Trigger(:,i),interval(:,i),'o','color',[0 0 1],'LineWidth',0.5,...
%         'MarkerSize',2,'DisplayName',strcat('$\hat{T}^{',num2str(i),'}_{k}$'))
% 
%     % Triggeri = Trigger(count(i),i)/count(i)*ones(1,Ts/dt+1); %最后一次触发的时刻/触发次数
%     Triggeri = Ts/count(i)*ones(1,Ts/dt+1); %仿真时间/触发次数
%     plot(0:dt:Ts,Triggeri,'linestyle','--','LineWidth',1,'MarkerSize',2)
%     hold on
%     % text(Ts+0.2,Trigger(count(i),i)/count(i)+0.2,num2str(roundn(Trigger(count(i),i)/count(i),-4)),'FontSize',10,'FontWeight','bold')%最后一次触发的时刻/触发次数
%     text(Ts+0.2,Ts/count(i)+0.4,num2str(roundn(Ts/count(i),-4)),'FontSize',10,'FontWeight','bold') %仿真时间/触发次数
%     text(Ts+0.2,Min(i),num2str(roundn(Min(i),-4)),'FontSize',10,'FontWeight','bold')
%     xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\hat{T}^{',num2str(i),'}_{\hat{k}}$'),'Interpreter','Latex')
%     set(gca,'YScale','log');
%     ylim([10^-3,10])
%     box on
% end
% 
% % sigma 每个单独一张小图展示
% figure(4)
% for j = 1:N
%     subplot(N,1,j)
%     plot(0:dt:Ts,sigma_all(1:1:end,j),'-','linewidth',1.5)
%     xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\sigma_{i',num2str(j),'}$'),'Interpreter','Latex');
%     box on
% end
% 
% % sigma一起展示
% figure(5)
% plot(0:dt:Ts,sigma_all(1:1:Ts/dt+1,1:1:end),'-','linewidth',1.5)
% xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\sigma_{i}$'),'Interpreter','Latex');
% legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','FontName','Times','FontSize',8,'NumColumns',3,'FontWeight','bold','Location','southeast')
% 
% 
% figure(6)
% plot(0:dt:Ts,ee_all(1:1:end,1),'-','linewidth',1.5)
% hold on
% plot(0:dt:Ts,eother_all(1:1:end,1),'--','linewidth',1.5)
% xlabel('$t$','Interpreter','Latex');ylabel(strcat('$f_1(t)$'),'Interpreter','Latex');
% box on

toc
disp(['运行时间: ',num2str(toc)]);

