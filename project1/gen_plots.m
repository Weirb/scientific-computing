function [] = gen_plots()

laplace_2d_banded_normal_time_iter_comparison_plot();
laplace_1d_banded_normal_time_comparison_plot();
laplace_1d_banded_iter_time_plot();
laplace_1d_iter_time_plot();
laplace_2d_plot();
laplace_1d_sol_plot();

matrix2_time_iter_plot();
matrix2_sol_plot();
matrix2_create_exact_solution_plot();

close all;
end

function [] = matrix2_create_exact_solution_plot()

n = 1000;
N = 0:1/(n-1):1;
m = [0, 1, 10];

A = zeros(n);
b = ones(n,1).*2.5;
figure();
hold on;
for k = 1:3
    for i = 1:n
        for j = 1:n
            if i == j
                A(i,j) = 2* (i)^2 + m(k);
            elseif abs(i-j) == 1
                A(i,j) = -(i)^2;
            end
        end
    end
    
    x = A\b;
    plot(N, x);
end

set(findall(gcf,'-property','FontSize'),'FontSize',16);
title_s = sprintf('Plot of MATLAB solution to second matrix\nfor different values of m.');
title(title_s);
xlabel('x');
ylabel('y');
legend('m=0', 'm=1', 'm=10');
grid on;

saveas(gcf, sprintf('figures/mat2-exact-sol'), 'epsc');


end

function [] = matrix2_time_iter_plot()

data = load('task3.2.4.8-2.dat');

for i = 0:3
    tmp = data(1+i*20:(i+1)*20, :);
    
    n = tmp(:,1);
    m = tmp(:,2);
    iter = tmp(:,2);
    time = tmp(:,3);
    
    % Iteration plot
    figure();
    plot(n, iter, 'ro');
    set(findall(gcf,'-property','FontSize'),'FontSize',16);
    title_s = sprintf('Plot of size of problem against iterations to converge\nfor the second matrix, with m=%.1f',m(1));
    title(title_s);
    xlabel('Size of matrix');
    ylabel('Iterations to converge');
    grid on;

    saveas(gcf, sprintf('figures/mat2-iter-%d',i), 'epsc');
    
    % Time plot
    figure();
    loglog(n, time, 'bx');
    set(findall(gcf,'-property','FontSize'),'FontSize',16);
    title_s = sprintf('Log-log plot of size of problem against time to converge\nfor the second matrix, with m=%.1f',m(1));
    title(title_s);
    xlabel('Log size of matrix');
    ylabel('Log time to converge');
    grid on;

    saveas(gcf, sprintf('figures/mat2-time-%d',i), 'epsc');
end

end

function [] = matrix2_sol_plot()

data = load('task3.2.4.8-1.dat');

n = 0;
k = 0;
for j = 1:12
    k = k + 1 + n;
    n = data(k,1);
    m = data(k,2);
    tmp = data(k+1 : k+n, :);
    
    i = tmp(:,1);
    ip1overnp1 = tmp(:,2);
    x = tmp(:,3);

    % Plot with the standard axis
    figure();
    plot(i, x);
    set(findall(gcf,'-property','FontSize'),'FontSize',16);
    title_s = sprintf('Solution plot for second matrix, n=%d, m=%.1f.', n,m);
    title(title_s);
    xlabel('i');
    ylabel('x_i');
    grid on;

    saveas(gcf, sprintf('figures/mat2-sol-1-%d-%.1f.eps', n, m), 'epsc');

    % Plot with alternate x axis
    figure();
    plot(ip1overnp1, x);
    set(findall(gcf,'-property','FontSize'),'FontSize',16);
    title_s = sprintf('Solution plot for second matrix, n=%d, m=%.1f.', n,m);
    title(title_s);
    xlabel('(i+1)/(n+1)');
    ylabel('x_i');
    grid on;

    saveas(gcf, sprintf('figures/mat2-sol-2-%d-%.1f.eps', n, m), 'epsc');
end

end

function [] = laplace_2d_banded_normal_time_iter_comparison_plot()

data = load('task3.3.4.dat');
n = data(:,1);

% Plot iterations to converge
figure();
hold on;

set(gca,'XScale','linear','YScale','linear');
plot(n, data(:,3), 'bx'); % MMatrix iterations
plot(n, data(:,5), 'ro'); % MBandedMatrix iterations

set(findall(gcf,'-property','FontSize'),'FontSize',16);
title('Comparison of execution time using MMatrix and MBandedMatrix classes.');
xlabel('Size of matrix');
ylabel('Time to converge');
legend('MMatrix', 'MBandedMatrix');
legend('Location', 'northwest');
grid on;

saveas(gcf, 'figures/2d-matrix-type-compare-iter', 'epsc');

% Plot time to converge
figure();
hold on;

set(gca,'XScale','log','YScale','log');
plot(n, data(:,2), 'bx'); % MMatrix timing
plot(n, data(:,4), 'rx'); % MBandedMatrix timing

set(findall(gcf,'-property','FontSize'),'FontSize',16);
title('Comparison of execution time for 2D Poisson problem using MMatrix and MBandedMatrix classes.');
xlabel('Size of matrix');
ylabel('Time to converge');
legend('MMatrix', 'MBandedMatrix');
legend('Location', 'northwest');
grid on;

saveas(gcf, 'figures/2d-matrix-type-compare-time', 'epsc');

end

function [] = laplace_1d_banded_normal_time_comparison_plot()

data_banded = load('task3.2.5.6.dat');
data_matrix = load('task3.2.4.7.dat');

figure();
hold on;

n1 = data_banded(:,1);
n2 = data_matrix(:,1);

k = min(numel(n1), numel(n2));
set(gca,'XScale','log','YScale','log');
plot(n1(1:k), data_matrix(1:k,3), 'bx');
plot(n1(1:k), data_banded(1:k,3), 'rx');

set(findall(gcf,'-property','FontSize'),'FontSize',16);
title('Comparison of execution time for 1D Poisson problem using MMatrix and MBandedMatrix classes.');
xlabel('Size of matrix');
ylabel('Time to converge');
legend('MMatrix', 'MBandedMatrix');
legend('Location', 'northwest');
grid on;

saveas(gcf, 'figures/1d-matrix-type-compare-time', 'epsc');

end

function [] = laplace_1d_banded_iter_time_plot()
data = load('task3.2.5.6.dat');

n = data(:,1);
iter = data(:,2);
time = data(:,3);

figure();
plot(n, iter, 'ro');
set(findall(gcf,'-property','FontSize'),'FontSize',16);
title('Plot of size of problem against iterations to converge for 1D Poisson problem using MBandedMatrix.');
xlabel('Size of matrix');
ylabel('Iterations to converge');
grid on;

saveas(gcf, 'figures/1d-banded-iter', 'epsc');

figure();
loglog(n, time, 'bx');
set(findall(gcf,'-property','FontSize'),'FontSize',16);
title('Log-log plot of size of problem against time to converge for 1D Poisson problem using MBandedMatrix.');
xlabel('Log size of matrix');
ylabel('Log time to converge');
grid on;

saveas(gcf, 'figures/1d-banded-time', 'epsc');

end

function [] = laplace_1d_iter_time_plot()

data = load('task3.2.4.7.dat');

n = data(:,1);
iter = data(:,2);
time = data(:,3);

% Plot iterations
figure();
plot(n, iter, 'ro');
set(findall(gcf,'-property','FontSize'),'FontSize',16);
title_s = sprintf('Plot of size of problem against iterations\nto converge for 1D Poisson problem.');
title(title_s);
xlabel('Size of matrix');
ylabel('Iterations to converge');
grid on;

saveas(gcf, 'figures/1d-matrix-iter', 'epsc');

% Plot time
figure();
loglog(n, time, 'bx');
set(findall(gcf,'-property','FontSize'),'FontSize',16);
title_s = sprintf('Log-log plot of size of problem against time\nto converge for 1D Poisson problem.');
title(title_s);
xlabel('Log size of matrix');
ylabel('Log time to converge');
grid on;

saveas(gcf, 'figures/1d-matrix-time', 'epsc');

end

function [] = laplace_2d_plot()

data = load('task3.3.5.dat');

figure();
surf(data);
set(findall(gcf,'-property','FontSize'),'FontSize',16);
title('3D surface plot of the solution to the 2D Poisson problem, n=25^2.')
xlabel('x');
ylabel('y');

saveas(gcf, 'figures/2d-surf', 'epsc');

end

function [] = laplace_1d_sol_plot()
    clear;
    data = load('task3.2.4.6.dat');

    n = 0;
    k = 0;
    for j = 1:3
        k = k + 1 + n;
        n = data(k,1);
        tmp = data(k+1 : k+n, :);

        i = tmp(:,1);
        ip1overnp1 = tmp(:,2);
        x = tmp(:,3);
        
        % Plot with the standard axis
        figure();
        plot(i, x);
        set(findall(gcf,'-property','FontSize'),'FontSize',16);
        title_s = sprintf('Solution plot for 1D Poisson problem, n=%d.', n);
        title(title_s);
        xlabel('i');
        ylabel('x_i');
        grid on;

        saveas(gcf, sprintf('figures/1d-matrix-sol-1-%d', n), 'epsc');
        
        % Plot with alternate x axis
        figure();
        plot(ip1overnp1, x);
        set(findall(gcf,'-property','FontSize'),'FontSize',16);
        title_s = sprintf('Solution plot for 1D Poisson problem, n=%d.', n);
        title(title_s);
        xlabel('(i+1)/(n+1)');
        ylabel('x_i');
        grid on;

        saveas(gcf, sprintf('figures/1d-matrix-sol-2-%d', n), 'epsc');
    end

end