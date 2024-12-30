close all
clear
clc
tic

%% initial conditions of graph
A = [0  1  0  0
    0  0  0  0
    0  0  0  1
    0  0  0  0];
B_u = [0 0;1 0;0 0;0 1];
C_y = [1 0 0 0;0 0 1 0];
ADJ  = [0  0.6  0.8  1.2  0  0  0  0
    0.6  0  0  0  1  0.8  0  0
    0.8  0  0  0.9  0  0  0  0
    1.2  0  0.9  0  0.8  0  1.5  0
    0  1  0  0.8  0  0.9  0  1.2
    0  0.8  0  0  0.9  0  0  0
    0  0  0  1.5  0  0  0  0
    0  0  0  0  1.2  0  0  0];
L = diag(sum(ADJ,2))-ADJ;
lambda = eig(L);
lambda2 = min(lambda(lambda>10^-10));
c = 1/lambda2;

N = length(L);
NI = [3 3 2 4 4 2 1 1];

nx = size(A,1);
nu = size(B_u,2);

%% initial conditions of matrix
Q = eye(nx);
P = are(A, B_u*B_u', Q); % are(A,B,C): A'*X + X*A - X*B*X + C = 0
K = -B_u'*P;

% event-triggered sigma setup
dt = 0.001;
Ts = 30;
sigma = 1*[1 1 1 1 1 1 1 1];
sigma_all = zeros(1+Ts/dt,N); % sigma information of all time
sigma_all_tmp = zeros(1,N);
dsigma_dt = zeros(1,N); % the dot{sigma}
omega = 0.1;
delta = 0*[3 3 2 4 4 2 1 1];


Time = []; % triggering time
e = zeros(nx,N); % error information

% relevant matrix setup: the information of agent i at time t, and the
% information of all time
x = [];
x_all = zeros(1+Ts/dt,nx*N);
hat_x = [];
hat_x_all = zeros(1+Ts/dt,nx*N);

% initial conditions of matrix
% x0 = 2*rand([1,nx*N])-1;
load x0_linear.mat
hat_x0 = x0;

%% simulation
% some calculation for control signals to update
syms s
G = expm(A*dt);
G_K = int(expm(A*s),0,dt)*B_u*K;
G_K = eval(G_K);
H = int(expm(A*s),0,dt)*B_u;
H = eval(H);
count = 0;

% simulation for time
for t = 0:dt:Ts
    % The aim of matrix_i is to get information of every agent at every
    % instant, i.e., matrix_i is a tool, to get the whole information
    count = count+1;
    matrix_x = reshape(x0, nx, N);
    matrix_hat_x = reshape(hat_x0, nx, N);

    % Leader_x = G*Leader_x; % update leader
    % Leader_x_all(count,:) = Leader_x;

    % simulation for agents
    for i = 1:N
        % event-triggered mechanism
        delta_x = [0;0;0;0];
        delta_hat_x = 0;
        for j = 1:N
            delta_x = delta_x+L(i,j)*c*matrix_hat_x(:,j);
            delta_hat_x = delta_hat_x+ADJ(i,j)*norm(K*(matrix_hat_x(:,i)-matrix_hat_x(:,j)))^2;
        end
        matrix_x(:,i) = G*matrix_x(:,i)+G_K*delta_x; %  update x, u = c*K*matrix_hat_r(:,i);
        e(:,i) = matrix_hat_x(:,i)-matrix_x(:,i); % update error information
        dsigma_dt(i) = -delta(i)*sigma(i)-4*L(i,i)*norm(K*e(:,i))^2+delta_hat_x;
        sigma(i) = sigma(i)+dsigma_dt(i)*dt;
        if sigma(i) < 0
            sigma(i) = 0;
        end

        if sigma(i) == 0 % event-triggered
            Time = [Time;t,i]; % triggering time, Time = [time1, agent i;time2, agent j]
            matrix_hat_x(:,i) = matrix_x(:,i); % update the information exchanged
        else
        end
    end

    % set new initial conditions for a new loop
    x0 = reshape(matrix_x, 1, nx*N);
    hat_x0 = reshape(matrix_hat_x, 1, nx*N);

    % get the whole information of the system
    x_all(count,:) = x0;
    hat_x_all(count,:) = hat_x0;
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

close all
set(0,'defaultfigurecolor','w')


% state x
figure(1)
for j = 1:nx
    subplot(nx,1,j)
    plot(0:dt:Ts,x_all(:,j:4:end),'-','linewidth',1.5)
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$x_{i',num2str(j),'}$'),'Interpreter','Latex');
    box on
end
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','FontName','Times','FontSize',8,'NumColumns',3,'FontWeight','bold','Location','southeast')


% trigger time
figure(2)
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
end

for i = 1:N
    subplot(N,1,i)
    hold on
    Mini = Min(i)*ones(1,Ts/dt+1);
    plot(0:dt:Ts,Mini,'linestyle',':','color',[1 0 0],'LineWidth',1,'MarkerSize',2)% itself
    hold on
    plot(Trigger(:,i),interval(:,i),'o','color',[0 0 1],'LineWidth',0.5,...
        'MarkerSize',2,'DisplayName',strcat('$\hat{T}^{',num2str(i),'}_{k}$'))

    Triggeri = Trigger(count(i),i)/count(i)*ones(1,Ts/dt+1);
    plot(0:dt:Ts,Triggeri,'linestyle','--','LineWidth',1,'MarkerSize',2)
    hold on
    text(Ts+0.2,Trigger(count(i),i)/count(i)+0.2,num2str(roundn(Trigger(count(i),i)/count(i),-4)),'FontSize',10,'FontWeight','bold')
    text(Ts+0.2,Min(i),num2str(roundn(Min(i),-4)),'FontSize',10,'FontWeight','bold')
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\hat{T}^{',num2str(i),'}_{\hat{k}}$'),'Interpreter','Latex')
    set(gca,'YScale','log');
    ylim([10^-3,10])
    box on
end
