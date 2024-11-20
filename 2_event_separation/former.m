close all
clear
clc

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

%% 计算控制器增益
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%C1
Q = eye(nx); 
P = are(A',C'*C,Q) % are(A,B,C): A'*X + X*A - X*B*X + C = 0
F = -P*C'

% syms t s
% Ad = expm(A*t);
% Bd = int(expm(A*(t-s))*F*Cy,s,[0,t]);
% G = Ad+Bd
% rhoG = abs(eig(G)) %
% Tc = 0.89;
% for i = 1:nx %数值求解，可能求不出来
%     Tc_temp = eval(solve(rhoG(i)==1,t,'Real',true));
%     Tc_temp = min(Tc_temp(Tc_temp>10^-15));
%     if ~isempty(Tc_temp)
%         Tc = min(Tc,Tc_temp)
%     end
% end
% for dt = 0.89:0.001:Tc+0.1 %作图求解
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


Tc = linspace(0.11,0.16,6); % the period of time-triggered mechanism
for i = 1:N
   [Ad1,Bd1,Cd,Dd]=ssdata(c2d(ss(A,F*C,[1 0 0],0),Tc(i),'zoh'));
    G = Ad1+Bd1;
    alpha(i) = 2/Tc(i)*log(1/max(abs(eig(G))))
end

%C2
Q = eye(nx); 
P = are(A,B*B',Q) % are(A,B,C): A'*X + X*A - X*B*X + C = 0
K = -B'*P

c = 1/lambda_m

delta = 0.9*min(min(abs(eig(Q)))/max(abs(eig(P))),min(alpha));
theta = 1;

Th = 10; %\hat{T}

eig(A+F*C);
eig(A+B*K);