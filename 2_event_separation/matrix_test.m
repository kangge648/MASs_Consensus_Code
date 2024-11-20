clear
clc

% L = [3 -1 0 -1 -1 0;
%      -1 1 0 0 0 0;
%      0 0 1 0 0 -1;
%     -1 0 0 1 0 0;
%     -1 0 0 0 2 -1;
%      0 0 -1 0 -1 2];  % laplace矩阵
%
% L(end-1, 1) % end代表最后一行

% check_zeta_tmp = [];
% a1 = [1;2;3];
% a2 = [4;5;6];
% check_zeta_tmp = [a1]
% check_zeta_tmp = [a1,a2] % 只要维度满足就可以串联矩阵，不过得注意变量的维度
% check_zeta_tmp(4,1) = 1;
% check_zeta_tmp(5,1) = 1; % 可以直接给不存在的地方添加矩阵元素，其他地方默认补0
% check_zeta_tmp(6,5) = 1;
%
% check_zeta_tmp(end+1,5) = 1

% 计算zeta初值
x0 = [0.2, 0.3, 0.5, 0.7, 0.6, -0.9, -0.6, 0, 0.8, -0.1, -0.4, -0.2, 0.6, -0.9, -0.6, 0, 0.8, -0.1];
r0 = [0.2, 0.3, 0.5, 0.7, 0.6, -0.9, -0.6, 0, 0.8, -0.1, -0.4, -0.2, 0.6, -0.9, -0.6, 0, 0.8, -0.1];
matrix_x0 = reshape(x0,3,6)
matrix_r0 = reshape(r0,3,6)

ADJ  = [0 1 0 0 1 0
        1 0 0 1 0 1
        0 0 0 0 1 0
        0 1 0 0 0 0
        1 0 1 0 0 1
        0 1 0 0 1 0];
L = diag(sum(ADJ,2))-ADJ;
sum = 0;
for i = 1:6
    for j = 1:6
        sum = sum+ADJ(i,j)*(matrix_x0(:,i)-matrix_x0(:,j));
    end
    answer(:,i) = matrix_r0(:,i)-sum
end

% A = [0  1  0
%     -1  0.1  1
%     0  0  0.1];
% dt = 0.0001;
% syms s
% C = [1 0 0];
% Q = eye(3);
% P = are(A',C'*C,Q); % are(A,B,C): A'*X + X*A - X*B*X + C = 0
% F = -P*C';
% G_F = int(expm(A*s),0,dt)*F*C;

% Ts = 10;
% dt = 0.01;
% 
% figure(1)
% test = ones(1,Ts/dt+1);
% 
% for i = 1:3
%     subplot(3,1,i)
%     hold on
%     plot(0:dt:Ts,test,'linestyle','--'); % plot的x和y需要数组一一对应
%     text(Ts+0.2,1,num2str(1))
% end

% A = [0  1  0
%     -1  0.1  1
%     0  0  0.1];
% eig(A)

% % 优化
% x = zeros(10,6);
% x0 = ones(1,6);
% x(2,:) = x0

