fig = openfig('Results\20210213044231506\no_1\I_mean.fig');
a = get(gca,'Children');
xdata = get(a, 'XData');
ydata = get(a, 'YData');
zdata = get(a, 'ZData');
size(zdata)