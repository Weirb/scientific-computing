function [] = gen_plots()

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

function [] = animate_plot()
    N = 32;
    data = load('plot_data.dat');
    y = data(1:N, :);
    x = linspace(0,2*pi,50);

    h = plot(NaN,NaN);
    axis([min(x) max(x) min(min(y)) max(max(y))]);

    for k = 1:N
        pause(0.1)
        set(h, 'XData', x, 'YData', y(k,:));
        drawnow
    end
end