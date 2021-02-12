function mapping_field(in, zou, x, y, z, E)
    text=append(in.map_dir, '\I_map_', int2str(zou), '.mat');
    save(text, 'x', 'y', 'z', 'E');

    if false
        M1 = mean(abs(E).^2, [1 2]);
        figure
        plot(z, M1(1,:), 'Linewidth',3);
        xlabel('Height (um)')
        ylabel('I')
        xlim([0 max(z)])
        ylim([0 max(abs(E).^2,[],'all')])
        set(gca,'Fontsize',12)
        set(gca,'XMinorTick','on','YMinorTick','on')
        set(gcf,'color','w');
        box on

        text = append(in.prefix,"Mean normalized field_", int2str(zou), ".png");
        saveas(gcf, text);
    end
end