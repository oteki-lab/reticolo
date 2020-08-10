load('.\results\back grating\period_2400_diam_27wav400_1200_nbpoints101_Fourier10.mat')

layers_name = layers(:,1);
layer_numss = layers(:,2);
legends = cell(1,length(layers));
Abs_array = cell(1,length(layers));
for index=1:length(layers)
    layer_nums = layer_numss{index};
    if length(layer_nums)>1
        legends(index) = append(layers_name(index),'(',int2str(layer_nums(1)),'-',int2str(layer_nums(length(layer_nums))),')');
    else
        legends(index) = append(layers_name(index),'(',int2str(layer_nums(1)),')');
    end
    Abs_temp = zeros(1,length(wavelength));
    for index2=layer_nums(1):layer_nums(length(layer_nums))
        Abs_temp = Abs_temp + Abs(index2,:);
    end
    Abs_array{index} = Abs_temp;
end
    
figure
hold on
for index=1:length(layers)
    plot(wavelength,cell2mat(Abs_array(:,index)), 'Linewidth',3);
end
hold off
legend(legends)
xlabel('\lambda (um)')
ylabel('Absorption')
xlim([min(wavelength) max(wavelength)])
ylim([0 1])
set(gca,'Fontsize',12)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w');
box on