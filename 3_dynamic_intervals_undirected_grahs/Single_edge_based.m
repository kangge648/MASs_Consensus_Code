close all
clear
clc
tic

% initial conditions of graph
ADJ  = [0  0.6  0.8  1.2  0  0  0  0
    0.6  0  0  0  1  0.8  0  0
    0.8  0  0  0.9  0  0  0  0
    1.2  0  0.9  0  0.8  0  1.5  0
    0  1  0  0.8  0  0.9  0  1.2
    0  0.8  0  0  0.9  0  0  0
    0  0  0  1.5  0  0  0  0
    0  0  0  0  1.2  0  0  0]
L = diag(sum(ADJ,2))-ADJ;
N = length(L);
NI = [3 3 2 4 4 2 1 1];
Sum_NI = sum(NI);
index = [0  1  2  3  0  0  0  0
    4  0  0  0  5  6  0  0
    7  0  0  8  0  0  0  0
    9  0  10  0  11  0  12  0
    0  13  0  14  0  15  0  16
    0  17  0  0  18  0  0  0
    0  0  0  19  0  0  0  0
    0  0  0  0  20  0  0  0]; % 循环中的edge下标，方便写代码
End = [4 6 4 7 8 5 4 5]; % 每个智能体最后一个邻居的标号

nx = 1;
nu = 1;

%% initial conditions of matrix
% event-triggered sigma setup
dt = 0.00005;
Ts = 30;
sigma = 1*ones(Sum_NI);
sigma_all = zeros(1+Ts/dt,Sum_NI); % sigma information of all time
dsigma_dt = zeros(1,Sum_NI); % the dot{sigma}
omega = 50;
delta = 0*ones(Sum_NI);
Time = []; % triggering time
e = zeros(nx,Sum_NI); % error information

% for i = 1:(N-1)
%     T(i) = 2/sqrt(16*L(i,i)^2-delta(i)^2)*(atan((2*omega+delta(i))/sqrt(16*L(i,i)^2-delta(i)^2))-atan((delta(i))/sqrt(16*L(i,i)^2-delta(i)^2)));
% end
% 
% for i = N:N
%     T(i) = 2/sqrt(16*L(i,i)^2-delta(i)^2)*(atan((2*omega+delta(i))/sqrt(16*L(i,i)^2-delta(i)^2))-atan((delta(i))/sqrt(16*L(i,i)^2-delta(i)^2)))
% end

% relevant matrix setup: the information of agent i at time t, and the
% information of all time
x = [];
x_all_edge = zeros(1+Ts/dt,nx*N);
hatx = [];
hatx_all = zeros(1+Ts/dt,nx*Sum_NI);

% initial conditions of matrix
% x0 = 2*rand([1,nx*N])-1;save('x0','x0');
load x0.mat
hatx0 = [x0(1) x0(1) x0(1) x0(2) x0(2) x0(2) x0(3) x0(3) x0(4) x0(4) x0(4) x0(4) x0(5) x0(5) x0(5) x0(5) x0(6) x0(6) x0(7) x0(8)];

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
    matrix_hat_x = reshape(hatx0, nx, Sum_NI);

    % simulation for agents
    for i = 1:N % node
        u = 0;
        for j = 1:N
            if ADJ(i,j) > 0
                u = u-ADJ(i,j)*(matrix_hat_x(:,index(i,j))-matrix_hat_x(:,index(j,i)));
            end
        end
        matrix_x(:,i) = matrix_x(:,i)+u*dt; %  update x, u = c*K*matrix_hat_r(:,i);
        for j = 1:N % edge
            if ADJ(i,j) > 0 % edge exists
                e(:,index(i,j)) = matrix_hat_x(:,index(i,j))-matrix_x(:,i); % update error information

                delta_x = 0;
                for k = 1:N
                    if ADJ(i,k) > 0
                        delta_x = delta_x+ADJ(i,j)*(matrix_hat_x(:,index(i,k))-matrix_hat_x(:,index(k,i)))^2;
                    end
                end

                fij = 4*ADJ(i,j)*e(:,index(i,j))^2-delta_x/NI(i)-omega*sigma(index(i,j));

                if fij >= 0 % event-triggered
                    Time = [Time;t,i,j]; % triggering time, Time = [time1, agent i;time2, agent j]
                    matrix_hat_x(:,index(i,j)) = matrix_x(:,i); % update the information exchanged
                    sigma_all(count,index(i,j)) = sigma(index(i,j));
                else
                    dsigma_dt(index(i,j)) = -delta(index(i,j))*sigma(index(i,j))+delta_x/NI(i)-4*ADJ(i,j)*e(:,index(i,j))^2;
                    sigma(index(i,j)) = sigma(index(i,j))+dsigma_dt(index(i,j))*dt;
                    sigma_all(count,index(i,j)) = sigma(index(i,j));
                end
            end
        end
    end

    % set new initial conditions for a new loop
    x0 = reshape(matrix_x, 1, nx*N);
    hatx0 = reshape(matrix_hat_x, 1, nx*Sum_NI);

    % get the whole information of the system
    x_all_edge(count,:) = x0;
    hatx_all(count,:) = hatx0;
end

%% plot figures
save('x_all_edge_single','x_all_edge')
close all
set(0,'defaultfigurecolor','w')

% state x
figure(1)
for j = 1:nx
    subplot(nx,1,j)
    plot(0:dt:Ts,x_all_edge(:,j:1:end),'-','linewidth',1.5)
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$x_{i}$'),'Interpreter','Latex');
    box on
end
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','Agent 7','Agent 8','FontName','Times','FontSize',8,'NumColumns',3,'FontWeight','bold','Location','southeast')

% state hatx
figure(2)
for j = 1:nx
    subplot(nx,1,j)
    plot(0:dt:Ts,hatx_all(:,j:1:end),'-','linewidth',1.5)
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\hat{x}_{i',num2str(j),'}$'),'Interpreter','Latex');
    box on
end
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','Agent 7','Agent 8','FontName','Times','FontSize',8,'NumColumns',3,'FontWeight','bold','Location','southeast')

% trigger time
% trigger记录了每个智能体触发的时间点，count记录了每个智能体触发的次数
% interval记录了每个智能体的触发时间间隔
% 暂时只记录了
figure(3)
Trigger = zeros(1,Sum_NI);
count = ones(1,Sum_NI);
for k = 1:length(Time)
    Trigger(count(index(Time(k,2),Time(k,3))),index(Time(k,2),Time(k,3))) = Time(k,1);
    count(index(Time(k,2),Time(k,3))) = count(index(Time(k,2),Time(k,3)))+1;
end
for i = 1:Sum_NI
    count(i) = count(i)-1;
end
interval = zeros(Sum_NI);
for i = 1:N
    for j = 1:N
        if ADJ(i,j) > 0
            interval(1,index(i,j)) = Trigger(1,index(i,j));
            for k = 2:count(index(i,j))
                interval(k,index(i,j)) = Trigger(k,index(i,j))-Trigger(k-1,index(i,j));
            end
        end
    end
end

Min = 100*ones(1,Sum_NI);
for i = 1:N
    for j = 1:N
        if ADJ(i,j) > 0
            if interval(j,index(i,j)) <= Min(index(i,j))
                Min(index(i,j)) = interval(j,index(i,j));
            end
        end
    end
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
    text(Ts+0.2,Min(i)-0.2,num2str(roundn(Min(i),-4)),'FontSize',10,'FontWeight','bold')
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\hat{T}^{',num2str(i),'}_{\hat{k}}$'),'Interpreter','Latex')
    set(gca,'YScale','log');
    ylim([10^-3,10])
    box on
end

% sigma
figure(4)
for j = 1:N
    subplot(N,1,j)
    plot(0:dt:Ts,sigma_all(1:1:end,j),'-','linewidth',1.5)
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\sigma_{',num2str(j),'}$'),'Interpreter','Latex');
    box on
end

toc
disp(['运行时间: ',num2str(toc)]);
