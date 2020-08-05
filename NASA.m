wavelength=0.27:0.1:4;
I=[];
exp_a=exp((-(wavelength-0.26503)/0.15994));
I=1000*sqrt(2)*(0.06977+7.0625*(1-exp_a).^2.28411.*exp_a);
clf
hold on
plot(wavelength,I,'Linewidth',3,'Linestyle','-.')
xlabel('\lambda (É m)')
ylabel('Intensity')
xlim([min(wavelength) max(wavelength)])
ylim([0 1000])
set(gca,'Fontsize',12)
legend('NASA')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w');
box on
