function [] = gen_plots()

animate_plot();
end

function [] = plot1()
    files = ["plot_data_10.dat";"plot_data_100.dat";"plot_data_200.dat"];
    hold on;
    for i = 1:3
        data = load(files(i,:));
        plot(data(:,1),data(:,2));
    end
    legend('10','100','200');
end

function [] = plot2(c, dt)
    data = load('plot_data.dat');
    x = data(:,1);
    hold on;
    plot(x, data(:,2), 'rx');
    plot(x, 1.5+sin(x-c*dt));
    legend('Approximation', 'Exact');
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

function [z] = square_wave(x)
if (0<= x) && (x <= 1)
    z = 1.5;
else
    z = 0;
end
end

function [] = animate_plot()
    N = 10;
    data = load('animation_data.dat');
    y = data(1:N, :);
    x = linspace(0,2*pi,50);
    
    figure();
    hold on;
    h = plot(NaN,NaN);
    k = plot(NaN,NaN);
    axis([min(x) max(x) min(min(y)) max(max(y))]);
    legend('Approximate', 'Exact');
    for i = 1:N
        pause(0.1)
%         ex = 1.5+sin(x-0.05*(i-1));
%         ex = square_wave(x-0.05*(i-1));
%         ex = (0 <= x) .* (x <= 1)
        ex = square(x);
%         set(h, 'XData', x, 'YData', y(i,:));
        set(k, 'XData', x, 'YData', ex);
        drawnow
    end
end