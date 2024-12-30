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

close all
figure(1)
load Linear_Trigger_analysis.mat
load Linear_Trigger_simulation.mat
load Linear_Average.mat
omega = [0.000001 0.000003 0.00001 0.00003 0.0001 0.0003 0.001 0.003 0.01 0.03 0.1 0.3 1 3 10 30 100];
N = 8;
TTT = zeros(17,N); %比值
for i = 1:17
    for j = 1:N
        TTT(i,j) = L_T(i,j)/L_TT(i,j);
    end
end
xtrick = [0.000001 0.00001 0.0001 0.001 0.01 0.1 1 10 100];%x轴上显示的坐标标签
%ytrick = [0.0001 0.001 0.01 0.1 1];
miny = -0.03;
xianzhi = [0.0000005,150,miny,0.4
    0.0000005,150,miny,0.4
    0.0000005,150,miny,0.4
    0.0000005,150,miny,0.2
    0.0000005,150,miny,0.3
    0.0000005,150,miny,0.4
    0.0000005,150,miny,0.6
    0.0000005,150,miny,0.6];
for i = 1:N
    subplot(2,4,i)
    plot(omega,L_T(:,i),'^-','linewidth',1)
    hold on
    plot(omega,L_TT(:,i),'v-','linewidth',1)
    hold on
    plot(omega,L_AV(:,i),'*-','linewidth',1)
    hold on
    %plot(omega,TTT(:,i),'o-','linewidth',1)%比值
    set(gca,'XScale','log')%改为log坐标格式
    %set(gca,'YScale','log')
    % axis([0.0000005,150,miny,1]);%限定范围
    % xlim([0.0000005,150])
    axis(xianzhi(i,:));
    xlabel('$\omega$','Interpreter','Latex');ylabel(strcat('Agent ${',num2str(i),'}$'),'Interpreter','Latex');
    xticks(xtrick);%yticks(ytrick);%yticklabels({'-4\pi','-2\pi','0','2\pi','4\pi'})
    grid on
end
% hh = legend('Guaranteed $T_i$','$\min_k\{T^i_k\}$','$\mathrm{min}_k\{T^i_k\}$','${T_i}/{\min_k\{T^i_k\}}$','Interpreter','Latex','Orientation','horizontal', 'location',[0.13,0.05,0.74,0.05],'FontName','Times New Roman','FontSize',10);
hh = legend('Guaranteed $T_i$','$\min_k\{T^i_k\}$','$\mathrm{mean}_k\{T^i_k\}$','Interpreter','Latex','Orientation','horizontal','FontName','Times New Roman','FontSize',10);
% set(hh,'Interpreter','latex','Location','SouthOutside') %设置legend为latex解释器显示分式
