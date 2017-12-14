function [] = gen_plots()

plot_burgers_all();

end

function [x] = square(x)
x = (0 <= x).*(x <= 1);
end

function [] = plot_burgers_square(filename, t)
    figure();
    hold on;
    grid on;
    
    F = @(x,u) square(x-u*t) - u;
    fimplicit(F, [0 2*pi],'-', 'LineWidth',2,'Color',[1,0.4,0],'MeshDensity',2000);
    
    data = load(filename);
%     plot(data(:,1),data(:,2), 'ko--', 'markers',8,'LineWidth',1);
    plot(data(:,1),data(:,2), 'ko', 'markers',8,'LineWidth',1);
    
    xlabel('x');
    ylabel('u(x,t)');
    legend('Exact solution', 'Approximate solution');
    axis([ 0 2*pi -0.2 1.4]);
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9, 6]);
    [~,name,~] = fileparts(filename);
    saveas(gcf, sprintf('figures/%s', name), 'epsc');
end

function [] = plot_burgers_sine(filename, t)
    figure();
    hold on;
    grid on;
    
    F = @(x,u) 1.5 + sin(x-u*t) - u;
    fimplicit(F, [0 2*pi],'-', 'LineWidth',2,'Color',[1,0.4,0]);
    
    data = load(filename);
    plot(data(:,1),data(:,2), 'ko', 'markers',8,'LineWidth',1);
    
    xlabel('x');
    ylabel('u(x,t)');
    legend('Exact solution', 'Approximate solution');
    axis([ 0 2*pi 0 3]);
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9, 6]);
    [~,name,~] = fileparts(filename);
    saveas(gcf, sprintf('figures/%s', name), 'epsc');
end

function [] = plot_burgers_all()

T = [0., 0.25, 0.5, 1., 1.25, 1.5, 1.75, 2.];
for i=1:numel(T)
    filename = sprintf('data/burgers_sin_t_%d.dat', i-1);
    plot_burgers_sine(filename, T(i));

    filename = sprintf('data/burgers_sq_t_%d.dat', i-1);
    plot_burgers_square(filename, T(i));
end
end

function [] = plot_advection_square(filename, t)
    figure();
    hold on;
    grid on;
    
    x = linspace(0,2*pi,1000);
    plot(x,square(x-t),'-', 'LineWidth',2,'Color',[1,0.4,0]);
    
    data = load(filename);
    plot(data(:,1),data(:,2), 'ko', 'markers',8,'LineWidth',1);
    
    xlabel('x');
    ylabel('u(x,t)');
    legend('Exact solution', 'Approximate solution');
    axis([ 0 2*pi -0.2 1.4]);
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9, 6]);
    [~,name,~] = fileparts(filename);
    saveas(gcf, sprintf('figures/%s', name), 'epsc');
end

function [] = plot_advection_sine(filename, t)
    figure();
    hold on;
    grid on;
    
    x = linspace(0,2*pi,1000);
    plot(x,1.5+sin(x-t),'-', 'LineWidth',2,'Color',[1,0.4,0]);
    
    data = load(filename);
    plot(data(:,1),data(:,2), 'ko', 'markers',8,'LineWidth',1);
    
    xlabel('x');
    ylabel('u(x,t)');
    legend('Exact solution', 'Approximate solution');
    axis([ 0 2*pi 0 3]);
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9, 6]);
    [~,name,~] = fileparts(filename);
    saveas(gcf, sprintf('figures/%s', name), 'epsc');
end

function [] = plot_advection_all()

T = [0., 0.25, 0.5, 1.];
for i=1:numel(T)
    filename = sprintf('data/advec_sin_t_%d.dat', i-1);
    plot_advection_sine(filename, T(i));
    
    filename = sprintf('data/advec_sq_t_%d.dat', i-1);
    plot_advection_square(filename, T(i));
end
end


function [] = plot_3d()
    N = 28;
    data = load('plot_data.dat');
    y = data(1:N, :)';
    x = linspace(0,2*pi,50);
    t = 0:0.1:0.1*(N-1);

    [T,X] = meshgrid(t,x);

    % Generate a surface plot of the solution
    surf(T,X,y);

    % Generate a contour plot of the solution
    % contour(T,X,y);

    xlabel('t');
    ylabel('x');
    zlabel('u(x,t)');
end

function [] = animate_plot_file(filename)
    data = load(filename);
    [N, M] = size(data);
    k = N;
    y = data(1:k, :);
    x = linspace(0,2*pi,M);
    
    h = plot(NaN,NaN);
    axis([min(x) max(x) min(min(y)) max(max(y))]);
    for i = 1:15:k
        pause(0.1);
        set(h, 'XData', x);
        set(h, 'YData', y(i,:));
        drawnow;
    end
end
