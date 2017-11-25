function [] = gen_plots()

plot1();

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