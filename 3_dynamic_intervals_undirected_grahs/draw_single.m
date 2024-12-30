close all
clear
clc

color = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    0.87      0.45      0];

set(0,'defaultfigurecolor','w')

load x_all_node_single.mat
load x_all_edge_single.mat
load sigma_all_node_single.mat
load ee_all_node_single.mat
load eother_all_node_single.mat
load trigger_node_single.mat
load trigger_count_node_single.mat

dt = 0.00005;
Ts = 30;
N = 8;

figure(1)
subplot(3,2,1);
for i=1:N
    plot(0:dt:Ts,x_all_node(1:1:end,i),'-','color',color(i,:),'linewidth',1)
    hold on
end
ylim([-1,1])
xlabel('$t$','Interpreter','Latex');ylabel(strcat('$x_{i}$'),'Interpreter','Latex');
title('A: States in the Node-based Case','FontName','Times','FontSize',10,'FontWeight','bold')
box on
pbaspect([1.2 1 1]) %控制图画比例

subplot(3,2,2);
for i=1:N
    plot(0:dt:Ts,x_all_edge(1:1:end,i),'-','color',color(i,:),'linewidth',1)
    hold on
end
xlabel('$t$','Interpreter','Latex');ylabel(strcat('$x_{i}$'),'Interpreter','Latex');
title('B: States in the Edge-based Case','FontName','Times','FontSize',10','FontWeight','bold')
box on
pbaspect([1.2 1 1])

subplot(3,2,[3,4]);
for i=1:N
    plot(0:dt:15,sigma_all(1:1:15/dt+1,i),'-','color',color(i,:),'linewidth',1)
    hold on
end
%ylim([0,1.7])
xlabel('$t$','Interpreter','Latex');ylabel(strcat('$\sigma_{i}$'),'Interpreter','Latex');
title('C: Dynamic Variables in the Node-based Case','FontName','Times','FontSize',10','FontWeight','bold')
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','Agent 7','Agent 8','FontName','Times','FontSize',8,'NumColumns',3,'FontWeight','bold','Location','southeast')

subplot(3,2,[5,6]);
plot(0:dt:15,ee_all(1:1:15/dt+1,1),':','color',[0.8500 0.3250 0.0980],'linewidth',1)
hold on
plot(0:dt:15,eother_all(1:1:15/dt+1,1),'-','color',[0 0.4470 0.7410],'linewidth',1)
hold on
plot(Trigger(1:1:16,1),ee_all(Trigger_count(1:1:16,1),1),'o','color',[0 0.4470 0.7410],'LineWidth',1,...
    'LineWidth',1,'MarkerSize',3); %32是空格的代码
xlabel('$t$','Interpreter','Latex');%ylabel(strcat('$\sigma_{i',num2str(j),'}$'),'Interpreter','Latex');
title('D: Triggering Function in the Node-based Case','FontName','Times','FontSize',10','FontWeight','bold')
box on
legend({strcat('$4 l_{11}{e}_{1}(t)^2$'), ...
    strcat('$\sum_{j\in \mathcal{N}_{1}}a_{1j}(\hat{x}_{1}-\hat{x}_{j})^{2}+\omega\sigma_1(t)$'),'${t}^1_k$'},'Interpreter','Latex','FontSize',10')

%%
close all
figure(2)
load Trigger_analysis.mat
load Trigger_simulation.mat
load Average.mat
omega = [0.000001 0.000003 0.00001 0.00003 0.0001 0.0003 0.001 0.003 0.01 0.03 0.1 0.3 1 3 10 30 100];
N = 8;
TTT = zeros(17,N); %比值
for i = 1:17
    for j = 1:N
        TTT(i,j) = T(i,j)/TT(i,j);
    end
end
xtrick = [0.000001 0.00001 0.0001 0.001 0.01 0.1 1 10 100];%x轴上显示的坐标标签
%ytrick = [0.0001 0.001 0.01 0.1 1];
for i = 1:N
    subplot(4,2,i)
    %plot(omega,T(:,i),'^-','linewidth',1)
    %hold on
    plot(omega,TT(:,i),'v-','linewidth',1)
    hold on
    plot(omega,AV(:,i),'*-','linewidth',1)
    hold on
    % plot(omega,TTT(:,i),'o-','linewidth',1)%比值
    set(gca,'XScale','log')%改为log坐标格式
    %set(gca,'YScale','log')
    %axis([0.0005,130,0.000005,7]);%限定范围
    axis([0.0000005,150,-1,4]);
    xlabel('$\omega$','Interpreter','Latex');ylabel(strcat('Agent ${',num2str(i),'}$'),'Interpreter','Latex');
    xticks(xtrick);%yticks(ytrick);%yticklabels({'-4\pi','-2\pi','0','2\pi','4\pi'})
    grid on
end
% hh = legend('Guaranteed $T_i$','$\min_k\{T^i_k\}$','$\mathrm{min}_k\{T^i_k\}$','${T_i}/{\min_k\{T^i_k\}}$','Interpreter','Latex','Orientation','horizontal', 'location',[0.13,0.05,0.74,0.05],'FontName','Times New Roman','FontSize',10);
hh = legend('$\min_k\{T^i_k\}$','$\mathrm{mean}_k\{T^i_k\}$','Interpreter','Latex','Orientation','horizontal','FontName','Times New Roman','FontSize',10);
% set(hh,'Interpreter','latex','Location','SouthOutside') %设置legend为latex解释器显示分式
