clear
clc

load ee_all.mat
load eother_all.mat
load interval.mat
load sigma_all.mat
load Time.mat
load Trigger.mat
load Trigger_count.mat
load x_all.mat
nx = 2;
N = 5;

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

box on
hold on
grid on

Ts = 30;
dt = 0.0001;
% state x
pos1 = [0.05 0.78125 0.23 0.16875];
subplot('Position',pos1)
plot(0:dt:Ts,x_all(1:1:end,1:2:end),'-','linewidth',1.5)
set(gca,'xticklabel',[])
set(gca,'XTick',[10,20,30]) 
set(gca,'YTick',[-8,-4,0]) 
% xlabel('$t$','Interpreter','Latex');
ylabel(strcat('${x}_{i1}$'),'Interpreter','Latex');
%box on
%hold on
pos1 = [0.05 0.6 0.23 0.16875];
subplot('Position',pos1)
plot(0:dt:Ts,x_all(1:1:end,2:2:end),'-','linewidth',1.5)
set(gca,'XTick',[10,20,30]) 
set(gca,'YTick',[-1,0,1]) 
xlabel('$t$','Interpreter','Latex');ylabel(strcat('${x}_{i2}$'),'Interpreter','Latex');
%box on
%hold on
%legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','FontName','Times','FontSize',8,'NumColumns',3,'FontWeight','bold','Location','southeast')


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

% for i = 1:N
%     pos1 = [0.66 0.95-i*0.16-(i-1)*0.0125 0.28 0.16];
%     subplot('Position',pos1)
%     Mini = Min(i)*ones(1,Ts/dt+1);
%     plot(0:dt:Ts,Mini,'linestyle',':','color',[1 0 0],'LineWidth',1,'MarkerSize',2)% itself
%     hold on
%     plot(Trigger(:,i),interval(:,i),'o','color',[0 0 1],'LineWidth',0.5,...
%         'MarkerSize',2,'DisplayName',strcat('$\hat{T}^{',num2str(i),'}_{k}$'))
% 
%     % Triggeri = Trigger(count(i),i)/count(i)*ones(1,Ts/dt+1); %最后一次触发的时刻/触发次数
%     Triggeri = Ts/count(i)*ones(1,Ts/dt+1); %仿真时间/触发次数
%     plot(0:dt:Ts,Triggeri,'linestyle','--','LineWidth',1,'MarkerSize',2)
%     hold on
%     % text(Ts+0.2,Trigger(count(i),i)/count(i)+0.2,num2str(roundn(Trigger(count(i),i)/count(i),-4)),'FontSize',10,'FontWeight','bold')%最后一次触发的时刻/触发次数
%     text(Ts+0.2,Ts/count(i),num2str(roundn(Ts/count(i),-4)),'FontSize',12,'FontWeight','bold') %仿真时间/触发次数
%     text(Ts+0.2,Min(i),num2str(roundn(Min(i),-4)),'FontSize',12,'FontWeight','bold')
%     if i == 5
%         xlabel('$t$','Interpreter','Latex');
%     else
%         set(gca,'xticklabel',[])
%     end
%     ylabel(strcat('${T}^{',num2str(i),'}_{{k}}$'),'Interpreter','Latex')
%     %set(gca,'XTick',[10,20,30]) 
%     set(gca,'YTick',[0.1,1]) 
%     set(gca,'YScale','log');
%     ylim([0.02,2])
%     box on
% end

ensure_min = [0.3657 0.3657 0.3657 0.3657 0.3657]*10^-3;
for i = 1:N
    pos1 = [0.65 0.9-i*0.15-(i-1)*0.0125 0.28 0.15];
    subplot('Position',pos1)
    plot(Trigger(:,i),interval(:,i),'o','color',[0 0 1],'LineWidth',0.5,...
        'MarkerSize',2,'DisplayName',strcat('$\hat{T}^{',num2str(i),'}_{k}$'))
    hold on
    % ensure_mini = ensure_min(i)*ones(1,Ts/dt+1);
    % stri = strcat(num2str(ensure_min(i)*10^4), '*10^{-4}');
    % plot(0:dt:Ts,ensure_mini,'linestyle','-.','LineWidth',1,'MarkerSize',2)
    % text(Ts+0.2,ensure_min(i),stri,'FontSize',12,'FontWeight','bold') %仿真时间/触发次数
    % hold on
    Mini = Min(i)*ones(1,Ts/dt+1);
    plot(0:dt:Ts,Mini,'linestyle',':','color',[1 0 0],'LineWidth',1,'MarkerSize',2)% itself
    text(Ts+0.2,Min(i),num2str(roundn(Min(i),-4)),'FontSize',12,'FontWeight','bold')
    hold on
    Triggeri = Ts/count(i)*ones(1,Ts/dt+1); %仿真时间/触发次数
    plot(0:dt:Ts,Triggeri,'linestyle','--','LineWidth',1,'MarkerSize',2)
    text(Ts+0.2,Ts/count(i),num2str(roundn(Ts/count(i),-4)),'FontSize',12,'FontWeight','bold') %仿真时间/触发次数
    hold on
    if i == 5
        xlabel('$t$','Interpreter','Latex');
    else
        set(gca,'xticklabel',[])
    end
    ylabel(strcat('${T}^{',num2str(i),'}_{{k}}$'),'Interpreter','Latex')
    set(gca,'XTick',[10,20,30]) 
    set(gca,'YTick',[0.1,1]) 
    set(gca,'YScale','log');
    ylim([0.01,2])
    box on
end

legend({strcat('$T_k^i$'), ...
    '$\min_k\{T_k^i\}$', '$\mathrm{mean}_k\{T_k^i\}$'},'NumColumns',3,'Interpreter','Latex','FontSize',10', 'Location','northeast')


Ts = 10;
% sigma 每个单独一张小图展示
pos1 = [0.35 0.6 0.24 0.35];
subplot('Position',pos1)
plot(0:dt:Ts,sigma_all(1:1:Ts/dt+1,1:1:end),'-','linewidth',1.5)
xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\sigma_{i}^\mathrm{st}$'),'Interpreter','Latex');
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','FontName','Times','FontSize',8,'NumColumns',1,'FontWeight','bold','Location','southeast')
%hold on

pos1 = [0.05 0.1 0.54 0.42];
subplot('Position',pos1)
plot(0:dt:Ts,eother_all(1:1:end,3),'-','linewidth',1.5)
hold on
plot(0:dt:Ts,ee_all(1:1:end,3),':','linewidth',1.5)
hold on
plot(Trigger(1:1:24,3),ee_all(Trigger_count(1:1:24,3),3),'o','color',[0 0.4470 0.7410],'LineWidth',1,...
    'LineWidth',1,'MarkerSize',4);
xlabel('$t$','Interpreter','Latex');%ylabel(strcat('$f_3$'),'Interpreter','Latex');
set(gca,'XTick',[2,4,6,8,10]) 
%ylim([0,1.2]);
legend({strcat('$\|K{e}_{3}\|^2$'), ...
    strcat(['$\frac{\omega_0\omega}{\left(\psi_{\mathrm{F}}\sqrt{\omega}+\frac{\xi_{\mathrm{F}}}{\phi_{\mathrm{F}}}\right)^2}\left\|\sum_{j\in\mathcal{N}_3}a_{3j}K(\hat{x}_3-\hat{x}_j)\right\|^2+\beta_3 \sigma^\mathrm{st}_3(t)$']),'${t}^3_k$'},'NumColumns',3,'Interpreter','Latex','FontSize',10')

% 
% 
% % % zp = BaseZoom();
% % % zp.run;
% % 
% % 
