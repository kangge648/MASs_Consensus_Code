close all
clear
clc
tic

%% initial conditions of graph
A = [0  1  0
    -1  0.1  1
    0  0 0.1];
B = [0;0;1];
C = [1 0 0];
ADJ  = [0 1 0 0 1 0
    1 0 0 1 0 1
    0 0 0 0 1 0
    0 1 0 0 0 0
    1 0 1 0 0 1
    0 1 0 0 1 0];
L = diag(sum(ADJ,2))-ADJ;
N = length(L);
lambda = eig(L);
lambda_m = min(lambda(lambda>10^-10)); % minimum non-zero eigenvalues of L

nx = size(A,1);
nu = size(B,2);
ny = size(C,1);
nr = nx;

%% initial conditions of matrix
% control design of S1
Q = eye(nx);
P = are(A',C'*C,Q); % are(A,B,C): A'*X + X*A - X*B*X + C = 0
F = -P*C';

Tc = linspace(0.11,0.16,6); % the period of time-triggered mechanism

% control design of S2
Q = eye(nx);
P = are(A, B*B', Q); % are(A,B,C): A'*X + X*A - X*B*X + C = 0
K = -B'*P;
c = 1/lambda_m;

% event-triggered eta setup
dt = 0.001;
Ts = 10;
eta0 = [5 6 7 8 9 10];
eta = eta0;
eta_all = zeros(1+Ts/dt,N); % eta information of all time
eta_all_tmp = zeros(1,N);
deta_dt = zeros(1,N); % the dot{eta}
rho = 0.001; % Young
v_0 = 0; % exp
eta_0 = 10; % initial value

Time = []; % triggering time
e = zeros(nx,N); % error information

% relevant matrix setup: the information of agent i at time t, and the
% information of all time
x = [];
x_all = zeros(1+Ts/dt,nx*N);
r = [];
r_all = zeros(1+Ts/dt,nx*N);
hat_r_all = zeros(1+Ts/dt,nx*N);
zeta = [];
zeta_all = zeros(1+Ts/dt,nx*N);
check_zeta_all = zeros(1+Ts/dt,nx*N);
tilde_x = [0;0;0];

% initial conditions of matrix
x0 = [0.2, 0.3, 0.5, 0.7, 0.6, -0.9, -0.6, 0, 0.8, -0.1, -0.4, -0.2, 0.6, -0.9, -0.6, 0, 0.8, -0.1];
r0 = [0.2, 0.3, 0.5, 0.7, 0.6, -0.9, -0.6, 0, 0.8, -0.1, -0.4, -0.2, 0.6, -0.9, -0.6, 0, 0.8, -0.1];
hat_r0 = r0;
zeta0 = [1.1 -0.6 -2.0 -0.4 -1.4 -0.5 -0.5 -2.9 -0.2 0.8 -2.3 -1.9 -0.7 1.0 0.7 0.0 0.8 -0.1];
check_zeta0 = zeta0;

%% simulation
% some calculation for control signals to update
syms s
% dt = 0.001;
% Ts = 10;
G = expm(A*dt);
G_F = int(expm(A*s),0,dt)*F*C;
G_F = eval(G_F);
G_K = int(expm(A*s),0,dt)*c*B*K;
G_K = eval(G_K);
H = int(expm(A*s),0,dt)*B;
H = eval(H);
count = 0;

% simulation for time
for t = 0:dt:Ts
    % The aim of matrix_i is to get information of every agent at every
    % instant, i.e., matrix_i is a tool, to get the whole information
    count = count+1;
    matrix_x = reshape(x0, nx, N);
    matrix_r = reshape(r0, nx, N);
    matrix_hat_r = reshape(hat_r0, nx, N);
    matrix_zeta = reshape(zeta0, nx, N);
    matrix_check_zeta = reshape(check_zeta0, nx, N);

    % simulation for agents
    for i = 1:N

        % time-triggered mechanism
        % update zeta
        tilde_x = [0;0;0];
        for j = 1:N
            tilde_x = tilde_x+L(i,j)*matrix_x(:,j);
        end
        matrix_zeta(:,i) = matrix_r(:,i)-tilde_x;
        if (t>0) && (mod(t, Tc(i)) == 0) % time-triggered
            matrix_check_zeta(:,i) = matrix_zeta(:,i);
        end

        % event-triggered mechanism
        % upadate \hat{r} and r
        if i == 1
            for j = 1:N
                matrix_hat_r(:,j) = G*matrix_hat_r(:,j); % update the estimator
            end
        end
        delta_r = [0;0;0];
        for j = 1:N
            delta_r = delta_r+L(i,j)*matrix_hat_r(:,j);
        end
        matrix_r(:,i) = G*matrix_r(:,i)+G_F*matrix_check_zeta(:,i)+G_K*delta_r;
        matrix_x(:,i) = G*matrix_x(:,i)+H*c*K*matrix_hat_r(:,i); %  update x, u = c*K*matrix_hat_r(:,i);
        e(:,i) = matrix_hat_r(:,i)-matrix_r(:,i); % update error information

        if eta(i) <= 0 % event-triggered
            eta(i) = eta0(i);
            Time = [Time;t,i]; % triggering time, Time = [time1, agent i;time2, agent j]
            matrix_hat_r(:,i) = matrix_r(:,i); % update the information exchanged
        else
            delta_Kr1 = 0;
            delta_Kr2 = 0;
            for j = 1:N
                delta_Kr1 = delta_Kr1+ADJ(i,j)*norm(K*(matrix_hat_r(:,i)-matrix_hat_r(:,j)))^2;
                delta_Kr2 = delta_Kr2+ADJ(i,j)*(e(:,i))'*(K')*K*(matrix_hat_r(:,i)-matrix_hat_r(:,j));
            end
            deta_dt(i) = -rho*eta(i)^2-v_0*eta(i)-(2*c*L(i,i))*(e(:,i)')*(K')*K*e(:,i)/((e(:,i)')*P*e(:,i))...
                +2*((e(:,i))'*P*A*matrix_r(:,i))/((e(:,i))'*P*e(:,i))+c*delta_Kr1/(2*(e(:,i))'*P*e(:,i))-2*c*eta(i)*delta_Kr2/((e(:,i)')*P*e(:,i)); % triggering function
            if deta_dt(i) <= 1000 % limit the size, too samll is not for use
                eta(i) = eta(i)+deta_dt(i)*dt;
            else
                eta(i) = eta(i)+eta_0*dt;
            end
            if eta(i) < 0
                    eta(i) = 0;
            end
        end

        % get the whole information of eta, i.e. eta_all: Time×N
        if i < N
            eta_all_tmp(i) = eta(i);
        else
            eta_all_tmp(i) = eta(i);
            eta_all(count,:) = eta_all_tmp;
            eta_all_tmp = [];
        end
    end

    % set new initial conditions for a new loop
    x0 = reshape(matrix_x, 1, nx*N);
    r0 = reshape(matrix_r, 1, nx*N);
    hat_r0 = reshape(matrix_hat_r, 1, nx*N);
    zeta0 = reshape(matrix_zeta, 1, nx*N);
    check_zeta0 = reshape(matrix_check_zeta, 1, nx*N);

    % get the whole information of the system
    x_all(count,:) = x0;
    r_all(count,:) = r0;
    hat_r_all(count,:) = hat_r0;
    zeta_all(count,:) = zeta0;
    check_zeta_all(count,:) = check_zeta0;
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

% event interval
figure(1)
Trigger1 = zeros(1,2);
for i = 1:length(Time)
    if Time(i,2) == 4
        Trigger1 = [Trigger1;Time(i,:)];
    end
end
interval1 = [];
for j = 2:(length(Trigger1))
    interval1 = [interval1;Trigger1(j,1),(Trigger1(j,1)-Trigger1(j-1,1))];
end
stem(interval1(:,1),interval1(:,2),'linewidth',1.5)

% state x
figure(2)
for j = 1:nx
    subplot(nx,1,j)
    plot(0:dt:Ts,x_all(:,j:3:end),'-','linewidth',1.5)
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$x_{i',num2str(j),'}$'),'Interpreter','Latex');
    box on
end
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','FontName','Times','FontSize',8,'NumColumns',3,'FontWeight','bold','Location','southeast')

% protocol state r
figure(3)
for j = 1:nr
    subplot(nr,1,j)
    plot(0:dt:Ts,r_all(:,j:3:end),'-','linewidth',1.5)
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$r_{i',num2str(j),'}$'),'Interpreter','Latex');
    box on
end
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','FontName','Times','FontSize',8,'NumColumns',3,'FontWeight','bold')

% protocol state \hat{r}
figure(4)
for j = 1:nr
    subplot(nr,1,j)
    plot(0:dt:Ts,hat_r_all(:,j:3:end),'-','linewidth',1.5)
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\hat{r}_{i',num2str(j),'}$'),'Interpreter','Latex');
    box on
end
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','FontName','Times','FontSize',8,'NumColumns',3,'FontWeight','bold')

% zeta
figure(5)
for j = 1:nx
    subplot(nx,1,j)
    plot(0:dt:Ts,zeta_all(:,j:3:end),'-','linewidth',1.5)
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\zeta_{i',num2str(j),'}$'),'Interpreter','Latex');
    box on
end
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','FontName','Times','FontSize',8,'NumColumns',3,'FontWeight','bold')

% \check{zeta}
figure(6)
for j = 1:nx
    subplot(nx,1,j)
    plot(0:dt:Ts,check_zeta_all(:,j:3:end),'-','linewidth',1.5)
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\check{\zeta}_{i',num2str(j),'}$'),'Interpreter','Latex');
    box on
end
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','FontName','Times','FontSize',8,'NumColumns',3,'FontWeight','bold')

% eta
figure(7)
subplot(2,1,1)
plot(0:dt:Ts,eta_all(1:1:end,4),'color',[0 0 1],'linewidth',1.5) 
hold on
eta4 = eta0(4)*ones(1,Ts/dt+1);
plot(0:dt:Ts,eta4,'linestyle',':','color',[0 0 1],'LineWidth',1,'MarkerSize',2)
hold on
xlim([0,Ts])
text(1,eta0(4)+1,'$\eta_{0}^4=8$','Interpreter','Latex','FontSize',10,'FontWeight','bold')
xlabel('$t$','Interpreter','Latex');ylabel('$\eta_{4}$','Interpreter','Latex')
subplot(2,1,2)
stem(interval1(:,1),interval1(:,2),'color',[1 0 0],'linewidth',1,'MarkerSize',5)
xlim([0,Ts])
xlabel('$t$','Interpreter','Latex');ylabel('$\check{T}_{k}^{4}$','Interpreter','Latex')

% trigger time
figure(8)
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

compare = [0.0095, 0.0053, 0.0612, 0.0821, 0.0221, 0.0352];

for i = 1:N
    subplot(N,1,i)
    hold on
    Mini = Min(i)*ones(1,Ts/dt+1);
    Comparei = compare(i)*ones(1,Ts/dt+1);
    plot(0:dt:Ts,Mini,'linestyle',':','color',[1 0 0],'LineWidth',1,'MarkerSize',2)% itself
    hold on
    plot(0:dt:Ts,Comparei,'linestyle',':','color',[0 1 0],'LineWidth',1,'MarkerSize',2)% compare
    text(Ts+0.2,compare(i),num2str(roundn(compare(i),-4)),'FontSize',10,'FontWeight','bold')
    hold on
    plot(Trigger(:,i),interval(:,i),'o','color',[0 0 1],'LineWidth',0.5,...
        'MarkerSize',2,'DisplayName',strcat('$\hat{T}^{',num2str(i),'}_{k}$'))

    Triggeri = Trigger(count(i),i)/count(i)*ones(1,Ts/dt+1);
    plot(0:dt:Ts,Triggeri,'linestyle','--','LineWidth',1,'MarkerSize',2)
    hold on
    text(Ts+0.2,Trigger(count(i),i)/count(i)+0.2,num2str(roundn(Trigger(count(i),i)/count(i),-4)),'FontSize',10,'FontWeight','bold')
%     if i == 1
%         text(Ts+1.5,Min(i),'00','FontSize',10,'FontWeight','bold')
%     end
%     if i == 3
%         text(Ts+1.9,Min(i),'0','FontSize',10,'FontWeight','bold')
%     end
%     if i == 6
%         text(Ts+1.9,Min(i),'0','FontSize',10,'FontWeight','bold')
%     end
    text(Ts+0.2,Min(i),num2str(roundn(Min(i),-4)),'FontSize',10,'FontWeight','bold')
    xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\hat{T}^{',num2str(i),'}_{\hat{k}}$'),'Interpreter','Latex')
    set(gca,'YScale','log');
    ylim([10^-3,10])
    box on
end

toc
disp(['运行时间: ',num2str(toc)]);
