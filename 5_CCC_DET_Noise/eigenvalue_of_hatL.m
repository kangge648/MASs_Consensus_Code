% 外部给定的n
% n = input('请输入矩阵的维度 n: ');
n = 100;
 
% 创建n阶全1列向量
ones_vector = ones(n, 1);
 
% 创建n阶单位矩阵
I = eye(n);
 
% 计算矩阵 I - (1/n) * ones_vector * ones_vector'
result_matrix = I - (1/n) * (ones_vector * ones_vector');
 
% 显示结果矩阵
disp(result_matrix);

% 显示特征值
eig(result_matrix)