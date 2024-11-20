% 生成一些示例数据
t = linspace(0, 30, 100);
x1 = sin(t) + randn(6, 100) * 0.1;
x2 = cos(t) + randn(6, 100) * 0.1;
x3 = sin(t/2) + randn(6, 100) * 0.1;
r = randn(6, 100) * 0.1;

figure;

% 图 (a)
subplot(4, 1, 1);
hold on;
for i = 1:6
    plot(t, x1(i, :), 'DisplayName', ['Agent ', num2str(i)]);
end
legend;
title('图 (a)');
hold off;

subplot(4, 1, 2);
hold on;
for i = 1:6
    plot(t, x2(i, :), 'DisplayName', ['Agent ', num2str(i)]);
end
legend;
title('图 (a)');
hold off;

subplot(4, 1, 3);
hold on;
for i = 1:6
    plot(t, x3(i, :), 'DisplayName', ['Agent ', num2str(i)]);
end
legend;
title('图 (a)');
hold off;

% 图 (b)
subplot(4, 1, 4);
hold on;
for i = 1:6
    plot(t, r(i, :), 'DisplayName', ['Agent ', num2str(i)]);
end
legend;
title('图 (b)');
hold off;

% 创建新figure用于图 (c) 和 (d)
figure;

% 图 (c)
subplot(2, 1, 1);
hold on;
for i = 1:6
    scatter(t, i * ones(size(t)), 'DisplayName', ['Agent ', num2str(i)]);
end
legend;
title('图 (c)');
hold off;

% 图 (d)
for i = 1:6
    subplot(2, 3, i + 3); % 这里使用 2x3 的布局
    hold on;
    plot(t, abs(x1(i, :)), 'o-', 'DisplayName', '||e||^2');
    plot(t, abs(x2(i, :)), 'r--', 'DisplayName', 'ω_i |y_i(t_k) + B_rC_yr_i(t_k)|^2 + θ_i');
    legend;
    title(['Agent ', num2str(i)]);
    hold off;
end