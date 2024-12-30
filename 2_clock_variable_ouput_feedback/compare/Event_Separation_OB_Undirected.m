%% OF,静态,协议状态ET,相对输出TT,无向图,无扰动,Leader-Follower & Leaderless

%% 
clc; close all

A = [0  1  0
    -1  0.1  1
     0  0  0.1];
Bu = [0;0;1];
Cy = [1 0 0];

ADJ  = [0 1 0 0 1 0
        1 0 1 1 0 1
        0 1 0 0 1 0
        0 1 0 0 1 0
        1 0 1 1 0 1
        0 1 0 0 1 0];
BB = [0  0  1  0  1  0]*0 %系数都为0表示leaderless,系数大于0表示leader-follower 
LL = diag(sum(ADJ,2))-ADJ
N = length(LL)
lambda = eig(LL+diag(BB));
lambda1 = min(lambda(lambda>10^-10));
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = size(A,1);
nu = size(Bu,2);
ny = size(Cy,1);
nr = nx;

%% 计算控制器增益
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%C1
Q = eye(nx); 
P = are(A',Cy'*Cy,Q) % are(A,B,C): A'*X + X*A - X*B*X + C = 0
F = -P*Cy'

syms t s
Ad = expm(A*t);
Bd = int(expm(A*(t-s))*F*Cy,s,[0,t]);
G = Ad+Bd
rhoG = abs(eig(G)) %
Tc = 10^20;
% for i = 1:nx %数值求解，可能求不出来
%     Tc_temp = eval(solve(rhoG(i)==1,t,'Real',true));
%     Tc_temp = min(Tc_temp(Tc_temp>10^-15));
%     if ~isempty(Tc_temp)
%         Tc = min(Tc,Tc_temp)
%     end
% end
% for dt = 0.89:0.0001:Tc+0.1 %作图求解
%    [Ad1,Bd1,Cd,Dd]=ssdata(c2d(ss(A,F*Cy,[1 0 0],0),dt,'zoh'));
%     G = Ad1+Bd1;
%     hold on;
%     plot(dt,max(abs(eig(G))),'x');grid on; drawnow;
%     G = Ad+Bd;G = double(subs(G,t,dt));
%     plot(dt,max(abs(eig(G))),'o');grid on; drawnow;
%     if max(abs(eig(G)))>=1
%         break
%     end
% end
% dt = dt-0.0001 %0.8966
% 人为选取满足条件的采样间隔
% !!!之前选择linspace(0.2,0.8,7)
% 这个采样间隔太大，而且不同智能体之间以及与S2之间差异太大，
% 导致仿真会出现非常chaos的结果，要仿真时间非常小才行
Tc = linspace(0.11,0.16,6);%Tc的数值精度不能比仿真步长小，否则simulink中的mod函数不工作
for i = 1:N
   [Ad1,Bd1,Cd,Dd]=ssdata(c2d(ss(A,F*Cy,[1 0 0],0),Tc(i),'zoh'));
    G = Ad1+Bd1;
    alpha(i) = 2/Tc(i)*log(1/max(abs(eig(G))));
end

%C2
Q = eye(nx); 
P = are(A,Bu*Bu',Q) % are(A,B,C): A'*X + X*A - X*B*X + C = 0
K = -Bu'*P

c = 1/lambda1

delta = 0.9*min(min(abs(eig(Q)))/max(abs(eig(P))),min(alpha))
theta = 1

Th = 10; %\hat{T}

eig(A+F*Cy)
eig(A+Bu*K)

%% 初始值
%x0 = 10*rand(N*nx,1)-5
% load x0.mat
x0 = 10*rand(N*nx,1)-5;
x0_leader = [1 1 1];

%% 仿真
I = eye(N);
for i = 1:N
    eval(strcat('ty0_',int2str(i),...
        ' = kron(LL(i,:),Cy)*x0;')); %sim文件中时间触发条件中的Memory模块需要
    eval(strcat('Kri_m_Krj_',int2str(i),...
        ' = kron(ones(N,1)*I(i,:),K)-kron(I,K);')); %sim文件中计算触发条件需要用到，避免仿真中进行运算
end
configSet = getActiveConfigSet('Example_Event_Separation_OB_Undirected');
set_param(configSet, 'StopTime', '30');
set_param(configSet, 'FixedStep', '0.0001');
sim('Example_Event_Separation_OB_Undirected')  

%% 画图
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

%% agent state x
t = x_r_u_tk.Time;
Data = x_r_u_tk.Data;
x = [];
for i = 1:N
    x = [x Data(:,(i-1)*(nx+nr+nu+1+nr+ny)+1:(i-1)*(nx+nr+nu+1+nr+ny)+nx)];
end
indx = 1:floor(size(x,1)/10000):size(x,1);
figure
if ~isempty(find(BB>0,1))
    x0 = x_leader.Data;
    for j = 1:nx
        subplot(nx,1,j)
        hold on
        plot(t(indx),x0(indx,j),'-','color',color(1,:),...
            'LineWidth',2,'DisplayName','Leader'); 
    end
end
for i = 1:N
    for j = 1:nx
        subplot(nx,1,j)
        hold on
        plot(t(indx),x(indx,(i-1)*nx+j),'-','color',color(i,:),...
            'LineWidth',1,'DisplayName',strcat('Agent',32,num2str(i))); %32是空格的代码
        box on
        xlabel('$t$','Interpreter','Latex')
        ylabel(strcat('$x_{i',num2str(j),'}$'),'Interpreter','Latex')
        drawnow;
    end
end

%% controller state r
t = x_r_u_tk.Time;
Data = x_r_u_tk.Data;
r = [];
for i = 1:N
    r = [r Data(:,(i-1)*(nx+nr+nu+1+nr+ny)+nx+1:(i-1)*(nx+nr+nu+1+nr+ny)+nx+nr)];
end
indx = 1:floor(size(r,1)/10000):size(r,1);
figure
for i = 1:N
    for j = 1:nr
        subplot(nr,1,j)
        hold on
        plot(t(indx),r(indx,(i-1)*nr+j),'-','color',color(i,:),'LineWidth',1,...
            'LineWidth',1,'DisplayName',strcat('Agent',32,num2str(i))); %32是空格的代码
        box on
        xlabel('$t$','Interpreter','Latex')
        ylabel(strcat('$r_{i',num2str(j),'}$'),'Interpreter','Latex')
        drawnow;
    end
end

%% controller state r(tk)
t = x_r_u_tk.Time;
Data = x_r_u_tk.Data;
rhat = [];
for i = 1:N
    rhat = [rhat Data(:,(i-1)*(nx+nr+nu+1+nr+ny)+nx+nr+nu+1+1:(i-1)*(nx+nr+nu+1+nr+ny)+nx+nr+nu+1+nr)];
end
indx = 1:floor(size(rhat,1)/10000):size(rhat,1);
figure
for i = 1:N
    for j = 1:nr
        subplot(nr,1,j)
        hold on
        plot(t(indx),rhat(indx,(i-1)*nr+j),'-','color',color(i,:),'LineWidth',1,...
            'LineWidth',1,'DisplayName',strcat('Agent',32,num2str(i))); %32是空格的代码
        box on
        xlabel('$t$','Interpreter','Latex')
        ylabel(strcat('$\hat{r}_{i',num2str(j),'}$'),'Interpreter','Latex')
        drawnow;
    end
end

%% intermediate variable zeta(tk) or xi(tk)
t = x_r_u_tk.Time;
Data = x_r_u_tk.Data;
czeta = [];
for i = 1:N
    czeta = [czeta Data(:,(i-1)*(nx+nr+nu+1+nr+ny)+nx+nr+nu+1+nx+1:(i-1)*(nx+nr+nu+1+nr+ny)+nx+nr+nu+1+nx+ny)];
end
indx = 1:floor(size(czeta,1)/10000):size(czeta,1);
figure
for i = 1:N
    for j = 1:ny
        subplot(ny,1,j)
        hold on
        plot(t(indx),czeta(indx,(i-1)*ny+j),'-','color',color(i,:),'LineWidth',1,...
            'LineWidth',1,'DisplayName',strcat('Agent',32,num2str(i))); %32是空格的代码
        box on
        xlabel('$t$','Interpreter','Latex')
        if ~isempty(find(BB>0,1))
            ylabel(strcat('$\check{\xi}_{i',num2str(j),'}$'),'Interpreter','Latex')
            drawnow;
        else
            ylabel(strcat('$C\check{\zeta}_{i',num2str(j),'}$'),'Interpreter','Latex')
            drawnow;
        end
            
    end
end

%% Sampling instants tki
Data = x_r_u_tk.Data;
t_k = [];
for i = 1:N
    t_k = [t_k Data(:,(i-1)*(nx+nr+nu+1+nr+ny)+nx+nr+nu+1)];
end
figure
for i = 1:N
    tki = t_k(:,i);
    tki = tki(find(diff(tki,1)>0)+1);
    subplot(N,1,i)
    hold on
    plot(tki,diff([0;tki],1),'o','color',color(1,:),'LineWidth',0.5,...
        'MarkerSize',2,'DisplayName',strcat('$\hat{T}^{',num2str(i),'}_{k}$'))
    plot([0,30],mean(diff([0;tki],1))*[1 1],'--','color',color(2,:),...
        'LineWidth',1,'MarkerSize',2,'DisplayName',strcat('Mean of $\hat{T}^{',num2str(i),'}_{k}$'))
    text(31,mean(diff([0;tki],1)),num2str(mean(diff([0;tki],1)),'%.4f'),'FontSize',10,'FontWeight','bold')
    plot([0,30],min(diff([0;tki],1))*[1 1],':','color',color(2,:),...
        'LineWidth',1,'MarkerSize',2,'DisplayName',strcat('Minimum of $\hat{T}^{',num2str(i),'}_{k}$'))
    text(31,min(diff([0;tki],1)),num2str(min(diff([0;tki],1)),'%.4f'),'FontSize',10,'FontWeight','bold')
    ylabel(strcat('$\hat{T}^{',num2str(i),'}_{\hat{k}}$'),'Interpreter','Latex')
    xlabel('$t$','Interpreter','Latex')
    set(gca,'YScale','log');
    ylim([10^-3,10])
    box on
end
set(findall(gcf,'-property','FontSize'),'FontSize',10)


%% Triggering functions f_i (only for leaderless)
t = x_r_u_tk.Time;
Data = x_r_u_tk.Data;
r = [];
rhat = [];
t_k = [];
for i = 1:N
    r = [r Data(:,(i-1)*(nx+nr+nu+1+nr+ny)+nx+1:(i-1)*(nx+nr+nu+1+nr+ny)+nx+nr)];
    rhat = [rhat Data(:,(i-1)*(nx+nr+nu+1+nr+ny)+nx+nr+nu+1+1:(i-1)*(nx+nr+nu+1+nr+ny)+nx+nr+nu+1+nr)];
    t_k = [t_k Data(:,(i-1)*(nx+nr+nu+1+nr+ny)+nx+nr+nu+1)];
end

f_e = [];
f_up = [];
for i = 1:N
    r_i = r(:,(i-1)*nr+1:i*nr);
    rhat_i = rhat(:,(i-1)*nr+1:i*nr);
    ehat_i = rhat_i-r_i;
    f_e = [f_e 4*LL(i,i)*vecnorm(ehat_i*K',2,2).^2]; %vecnorm(A,p,dim)
    
    f_up_i = theta*exp(-delta*t);
    for j = 1:N
        rhat_j = rhat(:,(j-1)*nr+1:j*nr);
        f_up_i = f_up_i+ADJ(i,j)*vecnorm(rhat_i*K'-rhat_j*K',2,2).^2;
    end
    f_up = [f_up f_up_i];
end

indx = 1:floor(size(f_e,1)/50000):size(f_e,1);
figure
for i = 1:1
    subplot(i,1,1)
    hold on
    plot(t(indx),f_up(indx,i),'-','color',color(i,:),'LineWidth',1,...
        'LineWidth',1,'DisplayName',strcat('Agent',32,num2str(i))); %32是空格的代码
    plot(t(indx),f_e(indx,i),':','color',color(i,:),'LineWidth',1,...
        'LineWidth',1,'DisplayName',strcat('Agent',32,num2str(i))); %32是空格的代码
    
    tki = t_k(:,i);
    tki_idx = find(diff(tki,1)>0)+1;
    plot(t(tki_idx),f_e(tki_idx,i),'o','color',color(i,:),'LineWidth',1,...
        'LineWidth',1,'MarkerSize',3,'DisplayName',strcat('Agent',32,num2str(i))); %32是空格的代码
    
    box on
    xlabel('$t$','Interpreter','Latex')
    legend({strcat('$\sum_{j\in \mathcal{N}_{i}}a_{ij}\Vert K\hat{r}_{i}+K\hat{r}_{j}\Vert ^{2}-\theta _{i}\mathrm{e}^{-\delta t}$'),...
        strcat('$4 l_{ii} \|K\hat{e}_{i}\|^2$'),...
        '$\hat{t}^i_k$'},...
        'Interpreter','Latex')
    drawnow;
end





