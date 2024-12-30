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

nx = size(A,1);
nu = size(B_u,2);
ny = size(C_y,1);

%% initial conditions of matrix
Q = eye(nx);
P = are(A, B_u*B_u', Q); % are(A,B,C): A'*X + X*A - X*B*X + C = 0
K = -B_u'*P;

% event-triggered sigma setup
dt = 0.0001;
Ts = 30;
sigma = 1*ones(Sum_NI);
sigma_all = zeros(1+Ts/dt,Sum_NI); % sigma information of all time
sigma_all_tmp = zeros(1,Sum_NI);
dsigma_dt = zeros(1,Sum_NI); % the dot{sigma}
omega = 1;
delta = 0*[3 3 3 3 3 3 2 2 4 4 4 4 4 4 4 4 2 2 1 1];

% adaptive law setup
eta3 = 0.5;
eta4 = 1;
c0 = 0.0001*ones(1,Sum_NI);
c = c0;
hatc = c0;
alpha = 1*ones(1,Sum_NI);
c_all_edge_based = zeros(1+Ts/dt,Sum_NI); % c information of all time
hatc_all = zeros(1+Ts/dt,Sum_NI); % hatc information of all time

Time = []; % triggering time
e = zeros(nx,Sum_NI); % error information

% relevant matrix setup: the information of agent i at time t, and the
% information of all time
x = [];
x_all_based_ada = zeros(1+Ts/dt,nx*N);
hat_x = [];
hat_x_all = zeros(1+Ts/dt,nx*Sum_NI);
% initial conditions of matrix
% x0 = 2*rand([1,nx*N])-1;save('x0_adaptive_linear','x0')
load x0_adaptive_linear.mat
hat_x0 = [x0(1:4) x0(1:4) x0(1:4) x0(5:8) x0(5:8) x0(5:8) x0(9:12) x0(9:12) x0(13:16) x0(13:16) x0(13:16) x0(13:16) x0(17:20) x0(17:20) x0(17:20) x0(17:20) x0(21:24) x0(21:24) x0(25:28) x0(29:32)];

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
    matrix_hat_x = reshape(hat_x0, nx, Sum_NI);

    % simulation for agents
    for i = 1:N
        if i == 1
            for j = 1:N
                for k = 1:N
                    if ADJ(j,k) > 0
                        matrix_hat_x(:,index(j,k)) = G*matrix_hat_x(:,index(j,k)); % update the estimator
                    end
                end
            end
        end
        % compute u_i
        for j = 1:N
            if ADJ(i,j)>0
                dcij_dt = alpha(index(i,j))*norm(K*(matrix_hat_x(:,index(i,j))-matrix_hat_x(:,index(j,i))))^2;
                c(index(i,j)) = c(index(i,j))+dcij_dt*dt;
                c_all_edge_based(count,index(i,j)) = c(index(i,j));
            end
        end
        u = 0;
        for j = 1:N
            if ADJ(i,j)>0
                u = u+K*c(index(i,j))*ADJ(i,j)*(matrix_hat_x(:,index(i,j))-matrix_hat_x(:,index(j,i)));
            end
        end

        matrix_x(:,i) = G*matrix_x(:,i)+H*u; %  update x;

        % event-triggered mechanism
        for j = 1:N
            if ADJ(i,j)>0
                delta_hat_x = 0;
                for k = 1:N
                    if ADJ(i,k)>0
                        delta_hat_x = delta_hat_x+ADJ(i,k)*norm(K*(matrix_hat_x(:,index(i,k))-matrix_hat_x(:,index(k,i))))^2;
                    end
                end
                e(:,index(i,j)) = matrix_hat_x(:,index(i,j))-matrix_x(:,i); % update error information

                fij = 4*(1+eta3*c(index(i,j)))*ADJ(i,j)*norm(K*e(:,index(i,j)))^2-eta4/NI(i)*delta_hat_x-omega*sigma(index(i,j));

                if fij >= 0 % event-triggered
                    Time = [Time;t,i,j]; % triggering time, Time = [time1, agent i;time2, agent j]
                    matrix_hat_x(:,index(i,j)) = matrix_x(:,i); % update the information exchanged
                    sigma_all(count,index(i,j)) = sigma(index(i,j));  % 没触发也得等于，否则变成0了
                else
                    dsigma_dt(index(i,j)) = -delta(index(i,j))*sigma(index(i,j))+eta4/NI(i)*delta_hat_x-4*(1+eta3*c(index(i,j)))*ADJ(i,j)*norm(K*e(:,i))^2;
                    sigma(index(i,j)) = sigma(i)+dsigma_dt(index(i,j))*dt;
                    sigma_all(count,index(i,j)) = sigma(index(i,j));
                end
            end
        end
    end

    % set new initial conditions for a new loop
    x0 = reshape(matrix_x, 1, nx*N);
    hat_x0 = reshape(matrix_hat_x, 1, nx*Sum_NI);

    % get the whole information of the system
    x_all_based_ada(count,:) = x0;
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
save('cij_all_edge_linear','c_all_edge_based')
save('x_all_edge_ada','x_all_based_ada')

% state x
figure(1)
plot(x_all_based_ada(:,1:4:end),x_all_based_ada(:,3:4:end),'-','linewidth',1.5)
xlabel('$t$','Interpreter','Latex');ylabel(strcat('$x_{i',num2str(j),'}$'),'Interpreter','Latex');
box on
% legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','FontName','Times','FontSize',8,'NumColumns',3,'FontWeight','bold','Location','southeast')

% state \hat{x}
figure(2)
plot(0:dt:Ts,x_all_based_ada(:,1:1:end),'-','linewidth',1.5)
xlabel('$t$','Interpreter','Latex');ylabel(strcat('${x}_{i}$'),'Interpreter','Latex');
text(10,3,strcat('$x_{i,3}$'),'Interpreter','Latex','FontSize',12,'Color',[0 0 0],'FontWeight','bold');
text(15,2,strcat('$x_{i,1}$'),'Interpreter','Latex','FontSize',12,'Color',[0 0 0],'FontWeight','bold');
text(20,0.5,strcat('$x_{i,4}$'),'Interpreter','Latex','FontSize',12,'Color',[0 0 0],'FontWeight','bold');
text(25,-0.1,strcat('$x_{i,2}$'),'Interpreter','Latex','FontSize',12,'Color',[0 0 0],'FontWeight','bold');
%h = annotation('textarrow',[10/30 0.5],[11/30 0.5]);
box on
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','Agent 7','Agent 8','FontName','Times','FontSize',8,'NumColumns',4,'FontWeight','bold')

% trigger time
figure(3)
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
    text(Ts+0.2,Trigger(count(i),i)/count(i),num2str(roundn(Trigger(count(i),i)/count(i),-4)),'FontSize',10,'FontWeight','bold')
    text(Ts+0.2,Min(i),num2str(roundn(Min(i),-4)),'FontSize',10,'FontWeight','bold')
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\hat{T}^{',num2str(i),'}_{\hat{k}}$'),'Interpreter','Latex')
    set(gca,'YScale','log');
    ylim([10^-3,10])
    box on
end

% c
figure(4)
plot(0:dt:Ts,c_all_edge_based,'-','linewidth',1.5)
xlabel('$t$','Interpreter','Latex');ylabel(strcat('$c_{ij}$'),'Interpreter','Latex');
box on
% legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','FontName','Times','FontSize',8,'NumColumns',3,'FontWeight','bold')

% sigma
figure(5)
for j = 1:N
    subplot(N,1,j)
    plot(0:dt:Ts,sigma_all(1:1:end,j),'-','linewidth',1.5)
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\sigma_{',num2str(j),'}$'),'Interpreter','Latex');
    box on
end

toc
disp(['运行时间: ',num2str(toc)]);
