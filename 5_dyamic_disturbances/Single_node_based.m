close all
clear
clc
tic

%% initial conditions of graph
ADJ  = [0 1 1 1 0 0
    1 0 0 0 0 1
    1 0 0 1 0 0
    1 0 1 0 1 0
    0 0 0 1 0 1
    0 1 0 0 1 0];
L = diag(sum(ADJ,2))-ADJ;
N = length(L);
NI = [3 2 2 3 2 2];

nx = 1;
nu = 1;

% dynamic event-triggered sigma setup
dt = 0.001;
Ts = 30;
sigma0 = 1*[1 1 1 1 1 1 1];
sigma = sigma0;
sigma_all = zeros(1+Ts/dt,N); % sigma information of all time
sigma_all_tmp = zeros(1,N);
dsigma_dt = zeros(1,N); % the dot{sigma}
delta = 0.1*NI;
omega = 10;
rho = 0.5*ones(1,N);
o1 = 0.5*ones(1,N);
o2 = 0.5*ones(1,N);

Time = []; % triggering time
Trigger_count = [];
tcount = ones(nx,N);
e = zeros(nx,N); % error information

% Guaranteed MIETs
Deta = zeros(1,N);
T = zeros(1,N);
for i = 1:(N-1)
    Deta(i) = 16*L(i,i)*(L(i,i)/rho(i)+(0.1+0.1)^2/o1(i))-delta(i)^2;
    T(i) = 2/sqrt(Deta(i))*(atan((2*omega+delta(i))/sqrt(Deta(i)))-atan(delta(i)/sqrt(Deta(i))));
end

for i = N:N
    Deta(i) = 16*L(i,i)*(L(i,i)/rho(i)+(0.1+0.1)^2/o1(i))-delta(i)^2
    T(i) = 2/sqrt(Deta(i))*(atan((2*omega+delta(i))/sqrt(Deta(i)))-atan(delta(i)/sqrt(Deta(i))))
end

% relevant matrix setup: the information of agent i at time t, and the
% information of all time
x = [];
x_all = zeros(1+Ts/dt,nx*N);
hatx = [];
hatx_all = zeros(1+Ts/dt,nx*N);
tilde_x = [];
tilde_x_all = zeros(1+Ts/dt,nx*N);
hat_tilde_x = [];
hat_tilde_x_all = zeros(1+Ts/dt,nx*N);

% Noise
v_all = zeros(1+Ts/dt,N);
w_all = zeros(1+Ts/dt,N);
dotw_all = zeros(1+Ts/dt,N);

% initial conditions of matrix
% load x0.mat
x0 = 20*rand([1,nx*N])-10;%save('x0','x0');
savex0 = x0;
hatx0 = x0; 
tildex0 = x0;
hat_tildex0 = x0;

% 作图看具体的触发过程
ee_all = zeros(1+Ts/dt,nx*N);
eother_all = zeros(1+Ts/dt,nx*N);
hathat_all = zeros(1+Ts/dt,nx*N);

% consensus error
chi = zeros(1+Ts/dt,N);

%% simulation
% some calculation for control signals to update
syms s
% dt = 0.001;
% Ts = 10;
G = exp(1*dt);
count = 0;

% simulation for time
for t = 0:dt:Ts
    % The aim of matrix_i is to get information of every agent at every
    % instant, i.e., matrix_i is a tool, to get the whole information
    count = count+1;
    matrix_x = reshape(x0, nx, N);
    matrix_hat_x = reshape(hatx0, nx, N);
    matrix_tilde_x = reshape(tildex0, nx, N);
    matrix_hat_tilde_x = reshape(hat_tildex0, nx, N);

    % error information
    if count > 1
        for i = 1:N
            v_all(count,i) = 2*rand(1)-1;  % -1到1之间
            v_all(count,i) = v_all(count,i)/10;
            dotw_all(count,i) = 2*rand(1)-1;  % -1到1之间
            dotw_all(count,i) = dotw_all(count,i)/10;
            w_all(count,i) = w_all(count-1,i)+dotw_all(count,i);
            if w_all(count,i) < -1/10
                w_all(count,i) = -1/10;
            end
            if w_all(count,i) > 1/10
                w_all(count,i) = 1/10;
            end
        end
    end

    % simulation for agents
    for i = 1:N
        u = 0;
        delta_x = 0;
        for j = 1:N
            u = u-L(i,j)*matrix_hat_tilde_x(:,j);
            delta_x = delta_x+ADJ(i,j)*(matrix_hat_tilde_x(:,i)-matrix_hat_tilde_x(:,j))^2;
            chi(count,i) = (chi(count,i)+matrix_x(:,i)-matrix_x(:,j))/N; % update consensus error
        end
        if chi(count,i) < 0
            chi(count,i) = -chi(count,i);
        end
        matrix_x(:,i) = matrix_x(:,i)+u*dt+v_all(count,i)*dt; %  update x, u = c*K*matrix_hat_r(:,i);
        matrix_tilde_x(:,i) = matrix_x(:,i)+w_all(count,i);
        e(:,i) = matrix_hat_tilde_x(:,i)-matrix_tilde_x(:,i); % update error information

        fi = 4*L(i,i)*e(:,i)^2-rho(i)*delta_x-omega*sigma(i)-o2(i);

        ee_all(count,i) = 4*L(i,i)*e(:,i)^2;
        eother_all(count,i) = rho(i)*delta_x+omega*sigma(i)+o2(i);
        hathat_all(count,i) = rho(i)*delta_x;

        if fi >= 0 % event-triggered
            Time = [Time;t,i]; % triggering time, Time = [time1, agent i;time2, agent j]
            tcount(i) = tcount(i)+1;
            Trigger_count(tcount(i)-1,i) = count;
            matrix_hat_tilde_x(:,i) = matrix_tilde_x(:,i); % update the information exchanged
        else
            dsigma_dt(i) = -delta(i)*sigma(i)+rho(i)*delta_x-4*L(i,i)*e(:,i)^2+o1(i);
            sigma(i) = sigma(i)+dsigma_dt(i)*dt;
        end

        % get the whole information of eta, i.e. eta_all: Time×N
        if i < N
            sigma_all_tmp(i) = sigma(i);
        else
            sigma_all_tmp(i) = sigma(i);
            sigma_all(count,:) = sigma_all_tmp;
            sigma_all_tmp = [];
        end
    end

    % set new initial conditions for a new loop
    x0 = reshape(matrix_x, 1, nx*N);
    hatx0 = reshape(matrix_hat_x, 1, nx*N);
    tildex0 = reshape(matrix_tilde_x, 1, nx*N);
    hat_tildex0 = reshape(matrix_hat_tilde_x, 1, nx*N);

    % get the whole information of the system
    x_all(count,:) = x0;
    hatx_all(count,:) = hatx0;
    tilde_x_all(count,:) = tildex0;
    hat_tilde_x_all(count,:) = hat_tildex0;
end

%% plot figures
color = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    0.87      0.45      0
    0.07      0.62      1
    1         0.41      0.16
    0.39      0.83      0.07];
color = [color;color;color;color];
set(0,'defaultfigurecolor','w')
close all

Fig = 1;
% save('x_all_node_single','x_all_node')
% save('sigma_all_node_single','sigma_all')
% save('ee_all_node_single','ee_all')
% save('eother_all_node_single','eother_all')

%% consenuss check
% state x
Ts = 30;
figure(Fig)
Fig = Fig+1;
for j = 1:nx
    subplot(nx,1,j)
    plot(0:dt:Ts,x_all(:,j:1:end),'-','linewidth',1.5)
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$x_{i}$'),'Interpreter','Latex');
    box on
end
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','FontName','Times','FontSize',8,'NumColumns',3,'FontWeight','bold','Location','southeast')

% state tildex
figure(Fig)
Fig = Fig+1;
for j = 1:nx
    subplot(nx,1,j)
    plot(0:dt:Ts,tilde_x_all(:,j:1:end),'-','linewidth',1.5)
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\tilde{x}_{',num2str(j),'}$'),'Interpreter','Latex');
    box on
end
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','FontName','Times','FontSize',8,'NumColumns',3,'FontWeight','bold','Location','southeast')

% state hat_tildex
figure(Fig)
Fig = Fig+1;
for j = 1:nx
    subplot(nx,1,j)
    plot(0:dt:Ts,hat_tilde_x_all(:,j:1:end),'-','linewidth',1.5)
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\hat{\tilde{x}}_{',num2str(j),'}$'),'Interpreter','Latex');
    box on
end
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','FontName','Times','FontSize',8,'NumColumns',3,'FontWeight','bold','Location','southeast')

% consensus error 每个单独一张小图展示
figure(Fig)
Fig = Fig+1;
for j = 1:N
    subplot(N,1,j)
    plot(0:dt:Ts,chi(1:1:end,j),'-','linewidth',1.5)
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\chi_{',num2str(j),'}$'),'Interpreter','Latex');
    box on
end

% Ts = 15;
% sigma 合并一起显示
figure(Fig)
Fig = Fig+1;
plot(0:dt:Ts,sigma_all(1:1:Ts/dt+1,1:1:end),'-','linewidth',1.5)
xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\chi_{i}$'),'Interpreter','Latex');
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','FontName','Times','FontSize',8,'NumColumns',3,'FontWeight','bold','Location','southeast')

%%
% trigger time
% trigger记录了每个智能体触发的时间点，count记录了每个智能体触发的次数
% interval记录了每个智能体的触发时间间隔
figure(Fig)
Fig = Fig+1;
Trigger = zeros(1,N);
count = ones(1,N);
for i = 1:length(Time)
    for j = 1:N
        if Time(i,2) == j
            Trigger(count(j),j) = Time(i,1);
            count(j) = count(j)+1;
        end
    end
end
for i = 1:N
    count(i) = count(i)-1;
end
interval = zeros(1,N);
for i = 1:N
    interval(1,i) = Trigger(1,i);
    for j = 2:count(i)
        interval(j,i) = Trigger(j,i)-Trigger(j-1,i);
    end
end

Min = 10*ones(1,N);
for i = 1:N
    for j = 1:count(i)
        if interval(j,i) <= Min(i)
            Min(i) = interval(j,i);
        end
    end
    % TT(z,i) = Min(i);
    % AV(z,i) = Ts/count(i);
end

for i = 1:N
    subplot(N,1,i)
    hold on
    Mini = Min(i)*ones(1,Ts/dt+1);
    plot(0:dt:Ts,Mini,'linestyle',':','color',[1 0 0],'LineWidth',1,'MarkerSize',2)% itself
    plot(Trigger(:,i),interval(:,i),'o','color',[0 0 1],'LineWidth',0.5,...
        'MarkerSize',2,'DisplayName',strcat('$\hat{T}^{',num2str(i),'}_{k}$'))

    % Triggeri = Trigger(count(i),i)/count(i)*ones(1,Ts/dt+1); %最后一次触发的时刻/触发次数
    Triggeri = Ts/count(i)*ones(1,Ts/dt+1); %仿真时间/触发次数
    plot(0:dt:Ts,Triggeri,'linestyle','--','LineWidth',1,'MarkerSize',2)
    hold on
    % text(Ts+0.2,Trigger(count(i),i)/count(i)+0.2,num2str(roundn(Trigger(count(i),i)/count(i),-4)),'FontSize',10,'FontWeight','bold')%最后一次触发的时刻/触发次数
    text(Ts+0.2,Ts/count(i)+0.4,num2str(roundn(Ts/count(i),-4)),'FontSize',10,'FontWeight','bold') %仿真时间/触发次数
    text(Ts+0.2,Min(i),num2str(roundn(Min(i),-4)),'FontSize',10,'FontWeight','bold')
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$T^{',num2str(i),'}_{k}$'),'Interpreter','Latex')
    set(gca,'YScale','log');
    ylim([10^-3,10])
    box on
end

% sigma 每个单独一张小图展示
figure(Fig)
Fig = Fig+1;
for j = 1:N
    subplot(N,1,j)
    plot(0:dt:Ts,sigma_all(1:1:end,j),'-','linewidth',1.5)
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\sigma_{i',num2str(j),'}$'),'Interpreter','Latex');
    box on
end


% Ts = 15;
% sigma 合并一起显示
figure(Fig)
Fig = Fig+1;
plot(0:dt:Ts,sigma_all(1:1:Ts/dt+1,1:1:end),'-','linewidth',1.5)
xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\sigma_{i}$'),'Interpreter','Latex');
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','FontName','Times','FontSize',8,'NumColumns',3,'FontWeight','bold','Location','southeast')


figure(Fig)
Fig = Fig+1;
plot(0:dt:Ts,ee_all(1:1:Ts/dt+1,1),':','color',[0.8500 0.3250 0.0980],'linewidth',1)
hold on
plot(0:dt:Ts,eother_all(1:1:Ts/dt+1,1),'-','color',[0 0.4470 0.7410],'linewidth',1)
hold on
plot(Trigger(1:1:14,1),ee_all(Trigger_count(1:1:14,1),1),'o','color',[0 0.4470 0.7410],'LineWidth',1,...
    'LineWidth',1,'MarkerSize',3); %32是空格的代码
xlabel('$t$','Interpreter','Latex');%ylabel(strcat('$\sigma_{i',num2str(j),'}$'),'Interpreter','Latex');
box on
legend({strcat('$4 l_{11}{e}_{i}(t)^2$'), ...
    strcat('$\rho_i\sum_{j\in \mathcal{N}_{1}}a_{1j}(\hat{x}_{i}-\hat{x}_{j})^{2}+o_f^i+\omega_i\sigma_i(t)$'),...
    '$\hat{t}^i_k$'},...
    'Interpreter','Latex')

%%
% process disturbances v
Ts = 0.5;
figure(Fig)
Fig = Fig+1;
for j = 1:N
    subplot(N,1,j)
    % plot(0:dt:Ts,v_all(1:1:end,j),'-','linewidth',1.5)
    stem(0:dt:Ts,v_all(1:1:Ts/dt+1,j),'color',[1 0 0],'linewidth',1,'MarkerSize',0.5)
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$v_{',num2str(j),'}$'),'Interpreter','Latex');
    box on
end

% measurement noise w
figure(Fig)
Fig = Fig+1;
for j = 1:N
    subplot(N,1,j)
    % plot(0:dt:Ts,v_all(1:1:end,j),'-','linewidth',1.5)
    stem(0:dt:Ts,w_all(1:1:Ts/dt+1,j),'color',[1 0 0],'linewidth',1,'MarkerSize',0.5)
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$w_{',num2str(j),'}$'),'Interpreter','Latex');
    box on
end


% save('trigger_node_single','Trigger')
% save('trigger_count_node_single','Trigger_count')
% save('Trigger_analysis','T')
% save('Trigger_simulation','TT')
% save('Average','AV')

%
% figure(7)
% plot(0:dt:Ts,sigma_all(1:1:Ts/dt+1,1),'-','linewidth',1.5)
% hold on
% plot(0:dt:Ts,hathat_all(1:1:Ts/dt+1,1),':','linewidth',1.5)
% hold on
% plot(0:dt:Ts,eother_all(1:1:Ts/dt+1,1),'-','color',[0 0.4470 0.7410],'linewidth',1.5)
% hold on
% xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\sigma_{i}$'),'Interpreter','Latex');

toc
disp(['运行时间: ',num2str(toc)]);
