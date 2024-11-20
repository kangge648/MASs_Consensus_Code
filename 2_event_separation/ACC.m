clear 
close all

A = [0, 1;
     -1, 0];
B = [1;
     0];
K = [-1.35, 0.41];
L = [3 -1 0 -1 -1 0;
     -1 1 0 0 0 0;
     0 0 1 0 0 -1;
    -1 0 0 1 0 0;
    -1 0 0 0 2 -1;
     0 0 -1 0 -1 2];  % laplace矩阵
P=[1.3522 -0.4142;
   -0.4142 1.9132];
lM_PhiP=1.352;

syms s
dt = 0.001;
G = expm(A*dt);
H = int(expm(A*s),0,dt)*B;
H = eval(H);

x0 = [0.2, 0.3, 0.5, 0.7, 0.6, -0.9, -0.6, 0, 0.8, -0.1, -0.4, -0.2];  % 初始状态
xtri0 = x0;
c= [0,1.1,0,0.9,0.8,0;
    1.1,0,0,0,0,0;
    0,0,0,0,0,1.2;
    0.9,0,0,0,0,0;
    0.8,0,0,0,0,1;
    0,0,1.2,0,1,0];%自适应参数初始值
%rho0=5+5*rand([1,6]);%记录参数rho的初始值
rho0=[1.05 0.9 0.8 1.2 1.1 1];
rho=rho0;%参数rho初始化
gamma=5;
xi = [gamma,gamma,gamma,gamma,gamma,gamma];%时钟变量初始值
xi_e = [gamma,gamma,gamma,gamma,gamma,gamma];%时钟变量2初始值
dcdt  = zeros(6,6);%自适应参数导数初始化
drhodt = zeros(1,6);%参数\rho导数初始化
dxidt = zeros(1,6); %时钟变量导数初始化
%此处还缺k，epsilon值的设置
k=ones(6,6);
epsilon=ones(1,6);
a = [0 2 0 2 2 0;
     2 0 0 0 0 0;
     0 0 0 0 0 2;
     2 0 0 0 0 0;
     2 0 0 0 0 2;
     0 0 2 0 2 0];%邻接矩阵
cL=zeros(6,6);%自适应参数和拓扑的混合矩阵
e = zeros(2,6);%误差值
M = []; % 记录时钟变量值
M_e=[]; % 记录时钟变量2初始值;
M_c=[]; % 记录自适应参数值：
E = []; % 记录误差
T_1 = []; %记录时间和触发agent的索引（条件1）
flag_T_1=zeros(1,6);%标志变量，置1说明本次xi(i)<=0的事件已经被记录过
flag_T_g=zeros(1,6);%0代表本个周期内事件for agent i没有触发
T_2 = []; % 记录时间和触发agent的索引（条件2）
x_his = []; % 记录实际状态值
xtri_his = []; % 记录触发状态值
Rho =[]; %记录自适应参数rho_i
flag_t=1;%全局参数，用来记录时间是否已经到10s;
M_tmp = [];

Ts = 10;
for t = 0:dt:Ts
    matrix_x = reshape(x0, 2, 6);
    matrix_xtri = reshape(xtri0, 2, 6);
    cL=zeros(6,6); %初始化混合矩阵cL
    Gamma=zeros(1,6);%初始化triggering function矩阵Gamma
    dcdt  = zeros(6,6);%自适应参数导数初始化
    drhodt = zeros(1,6);%参数\rho导数初始化
    dxidt = zeros(1,6); %时钟变量导数初始化
    dxi_edt = zeros(1,6);%时钟变量2导数初始化
    for i = 1:length(L)
        if xi(i) <= 0
            T_1 = [T_1;t,i]; % triggering time
            matrix_xtri(:,i) = matrix_x(:,i);   % 事件发生，预测值更新为实际值.
            T_2 = [T_2; t, i]; %记录第二段条件被达成的事件
            xi(i) = gamma; %拨回时钟变量到初始值
        else
            %正常更新时钟变量
            sumcijaij=0;%\sum_{j=1}^N(a_{ij}c_{ij})
            ihat_jhat=0;%为新版时钟变量的计算而准备的
            sum_1=0;%新版时钟变量计算中所使用的第一项
            sum_2=0;%新版时钟变量计算中所使用的第二项
            for j=1:length(L)
                    sumcijaij=sumcijaij+a(i,j)*c(i,j);
                    ihat_jhat=matrix_xtri(:,i)-matrix_xtri(:,j);%为新版时钟变量的计算准备的
                    sum_1=sum_1+c(i,j)*a(i,j)*norm(K*ihat_jhat)^2;%为新版时钟变量的计算而准备的
                    sum_2=sum_2+c(i,j)*a(i,j)*(e(:,i)')*(K'*K)*(ihat_jhat);%为新版时钟变量的计算而准备的
            end
%                 dxidt(i)=-lM_PhiP*(rho(i)+xi(i)+2*sumcijaij*(xi(i)+1)^2);%更新时钟变量1（旧版）
                 dxidt(i)=sum_1/(2*(e(:,i))'*P*e(:,i)) - ((xi(i)+rho(i))*(e(:,i)')*(K')*K*e(:,i))/((e(:,i)')*P*e(:,i)) - 2*(xi(i)+1)*sum_2/((e(:,i)')*P*e(:,i));%更新时钟变量1（新版）
                 if dxidt(i)<=1e3%如果上升率太大，那么不要起作用了
                 xi(i)=xi(i)+dt*dxidt(i);
                 end
                
%                 dxi_edt(i)=-lM_PhiP*(rho(i)+xi_e(i)+4*sumcijaij*(xi_e(i)^2+1));%更新时钟变量2
%                 xi_e(i)=xi_e(i)+dt*dxi_edt(i);
        end

        if i < length(L)
            M_tmp = [M_tmp, xi(i)];
        else
            M_tmp = [M_tmp, xi(i)];
            M = [M; M_tmp];%记录时钟变量1
            M_tmp = [];
        end

%         M_e=[M_e;xi_e];%记录时钟变量2
        
        x = matrix_x(:,i);
        xtri = matrix_xtri(:,i);
        %此处需要做一些处理，把自适应参数加进去
        for j=1:length(L)
            cL(i,i) = cL(i,i)+c(i,j)*a(i,j)/L(i,i);%这里的L(i,i)是N_i
            if j~=i
            cL(i,j) = c(i,j)*a(i,j);
            end
        end
        %处理完毕
        u = K(1) * cL(i,:).* L(i,:) * matrix_xtri(1,:)' + ...
            K(2) * cL(i,:).* L(i,:) * matrix_xtri(2,:)';
        x = G * x + H * u;  % 实际值更新
      
        matrix_x(:,i) = x; %将更新完的实际值传入原矩阵
       
        
        %更新自适应参数\rho
        drhodt(i) = epsilon(i)*norm(e(:,i)'*K')^2;
%           if t>=10 && flag_t==1 && i==length(L)
%               rho10=rho;%记录下rho(i)在t=10时的值
%               rho(i)= rho(i)+dt*drhodt(i)-exp(-1/rho10(i)*(t-10))*dt;
%               flag_t=0;
%           elseif  t>=10 && flag_t==1
%               rho10=rho;
%               rho(i)= rho(i)+dt*drhodt(i)-exp(-1/rho10(i)*(t-10))*dt;
%           elseif  t>=10
%               rho(i)= rho(i)+dt*drhodt(i)-exp(-1/rho10(i)*(t-10))*dt;
%           else
              rho(i)=rho(i)+dt*drhodt(i);
%           end
        
           for j=1:length(L)
                dcdt(i,j)=k(i,j)*a(i,j)*norm((matrix_xtri(:,i)-matrix_xtri(:,j))'*K')^2;%旧版自适应参数更新
                c(i,j)=c(i,j)+dcdt(i,j)*dt;
           end    
    end
    M_c =[M_c;c];%记录自适应参数c
    
     %在更新的时间先后上有一定问题，如果希望cij=cji，那么必须使用同一个matrix_xtri来更新cij和cji.
    %更新自适应参数cij(新)
%         if t>=10 && flag_t==1 
%             cij10=c;
%             for j=1:length(L)
%                 if c(i,j)>=1e-4
%                 dcdt(i,j)=k(i,j)*a(i,j)*norm((matrix_xtri(:,i)-matrix_xtri(:,j))'*K')^2;%旧版自适应参数更新
%                 c(i,j)=c(i,j)+dcdt(i,j)*dt-exp(-1/cij10(i,j)*(t-10))*dt;
%                 end
%             end
%         elseif t>=10
%             for j=1:length(L)
%                 if c(i,j)>=1e-4
%                 dcdt(i,j)=k(i,j)*a(i,j)*norm((matrix_xtri(:,i)-matrix_xtri(:,j))'*K')^2;%旧版自适应参数更新
%                 c(i,j)=c(i,j)+dcdt(i,j)*dt-exp(-1/cij10(i,j)*(t-10))*dt;
%                 end
%             end
%         else
            
 
    
    
    
    
%    更新预测值，放在agent循环之后进行更新，确保每个agent所使用的matrix_xtri都是相同的
%     for i=1:length(L)
%         if flag_T_g(i) ==1 %标志着agent i的事件在本周期内触发了
%             matrix_xtri(:,i)=matrix_x(:,i);%按照实际值更新
%             flag_T_g(i) = 0;%更新完毕，把标志变量拨回去
%         else%未触发agent i的事件
            matrix_xtri = G * matrix_xtri;
%         end
%     end
    
    %更新预测误差e，放在agent循环之后进行更新
    e = matrix_xtri-matrix_x;
    
    
    
    
    x0 = reshape(matrix_x, 1, 12);
    xtri0 = reshape(matrix_xtri, 1, 12);
    E = [E; e];
    x_his = [x_his; x0];
    xtri_his = [xtri_his; xtri0];
    Rho = [Rho;rho];
end

% 画图
figure(1)
plot(0:dt:Ts,(x_his(:,1)-x_his(:,7)),'--');
hold on
plot(0:dt:Ts,xtri_his(:,1:2:2*length(L)),'k');
xlabel('time');ylabel('x_1');

figure(2)
plot(0:dt:Ts,x_his(:,1:2:end),'-')
hold on
plot(0:dt:Ts,xtri_his(:,1:2:end),'k')
xlabel('time');ylabel('x_i^1');

%用来画event interval
figure(3)
T=zeros(1,2);
for i=1:length(T_1)
    if T_1(i,2)==1
        T=[T;T_1(i,:)];
    end
end
interval=[];
for j=2:(length(T))
    interval=[interval;T(j,1),(T(j,1)-T(j-1,1))];
end
stem(interval(:,1),interval(:,2),'linewidth',1.5)

figure(4)
% for i = 1:length(L)
%     subplot(3,2,i)
%     plot(0:dt:Ts,M(i:6:end,i),'r--')
%     hold on
% %     plot(0:dt:Ts,M_e(i:6:end,i),'b--')
% %     hold on
% %     plot(0:dt:10,E(i:6:end,i),'k')
% end
 plot(0:dt:Ts,M(1:1:end,1),'r--')
 axis([0 10 0 50])
%      hold on
%     plot(0:dt:Ts,M_e(2:6:end,i),'b--')
    
figure(5)
plot(0:dt:Ts,Rho(:,1:end));
hold on
xlabel('time');ylabel('\mu_i');

%用来画cij
figure(6)
subplot(2,1,1)
plot(0:dt:Ts,Rho(:,1:end));
hold on
xlabel('t');ylabel('\mu_i');
grid
subplot(2,1,2)
plot(0:dt:Ts,M_c(1:6:end,2));
hold on
plot(0:dt:Ts,M_c(1:6:end,4));
hold on
plot(0:dt:Ts,M_c(1:6:end,5));
hold on
plot(0:dt:Ts,M_c(3:6:end,6));
hold on
plot(0:dt:Ts,M_c(5:6:end,6));
grid
xlabel('t');ylabel('c_{ij}');

figure(7)
subplot(2,1,1)
plot(0:dt:Ts,x_his(:,1:2:end),'-','linewidth',1.5)
xlabel('t');ylabel('x_{i}^1');
grid
subplot(2,1,2)
plot(0:dt:Ts,x_his(:,2:2:end),'-','linewidth',1.5)
xlabel('t');ylabel('x_{i}^2');
grid

figure(8)
subplot(2,1,1)
plot(0:dt:Ts,M(1:1:end,1),'r-','linewidth',1.5)
hold on
scatter(interval(:,1),zeros(length(interval),1),25,'filled','g')
xlabel('t');ylabel('\xi_1')
grid
axis([0 4 0 50])
subplot(2,1,2)
stem(interval(:,1),interval(:,2),'linewidth',1.5)
xlabel('t');ylabel('Event intervals of agent 1')
axis([0 4 0 1])
grid
