function [] = gen_plots()


hold on;
data = load('data/burgers_sin_t_1.dat');
plot(data(:,1),data(:,2), 'ko', 'markers',8);

x = linspace(0,2*pi,1000); 

t=0.5; 
u=(2*sinh(x))./(cosh(x)-exp(-t)); 
plot(x,u,'--'); 
% for t=0.2:0.5:1.2 
%     u=(-2*sinh(x))./(cosh(x)-exp(-t)); 
%     plot(x,u); 
% end

end

function [x] = square(x)
x = (0 <= x).*(x <= 1);
end

function [] = plot_burgers_sine(filename, t)
    figure();
    hold on;
    grid on;
    
%     x = linspace(0,2*pi,1000);
%     plot(x,1.5+sin(x-t),'-', 'LineWidth',2,'Color',[1,0.4,0]);
    
    data = load(filename);
    plot(data(:,1),data(:,2), 'ko', 'markers',8);
    
    xlabel('x');
    ylabel('u(x,t)');
    legend('Exact solution', 'Approximate solution');
    axis([ 0 2*pi 0 3]);
    set(findall(gcf,'-property','FontSize'),'FontSize',20);
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9, 6]);
    [~,name,~] = fileparts(filename);
%     saveas(gcf, sprintf('figures/%s', name), 'epsc');
end

function [] = plot_advection_square(filename, t)
    figure();
    hold on;
    grid on;
    
    x = linspace(0,2*pi,1000);
    plot(x,square(x-t),'-', 'LineWidth',2,'Color',[1,0.4,0]);
    
    data = load(filename);
    plot(data(:,1),data(:,2), 'ko', 'markers',8);
    
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
    plot(data(:,1),data(:,2), 'ko', 'markers',8);
    
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
% SINE WAVES
plot_advection_sine('data/advec_sin_t_0.dat', 0);
plot_advection_sine('data/advec_sin_t_1.dat', 0.25);
plot_advection_sine('data/advec_sin_t_2.dat', 0.5);
plot_advection_sine('data/advec_sin_t_3.dat', 1);

% SQUARE WAVES
plot_advection_square('data/advec_sq_t_0.dat', 0);
plot_advection_square('data/advec_sq_t_1.dat', 0.25);
plot_advection_square('data/advec_sq_t_2.dat', 0.5);
plot_advection_square('data/advec_sq_t_3.dat', 1);
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
