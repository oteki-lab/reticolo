function filename = save_cross_section(in, horizontal, vertical, contour, cs, E, horizontal_label, vertical_label, diagonal_label, horizontal_lim, vertical_lim, diagonal_lim)
%%%% The separation between the layers is in white (indice contains the position of all indices)
figure
hP=surf(horizontal,vertical,repmat(cs,size(E)),abs(E).^2);   %, 'AlphaData',gradient(abs(E).^2),'FaceAlpha','interp');
shading interp;
colormap(jet);
colorbar;
caxis([0 Inf]);
view(2)
hP.DataTipTemplate.DataTipRows(1).Label = horizontal_label;
hP.DataTipTemplate.DataTipRows(2).Label = vertical_label;
hP.DataTipTemplate.DataTipRows(3).Label = diagonal_label;
hP.DataTipTemplate.DataTipRows(4).Label = '|E|^2';
hP.DataTipTemplate.DataTipRows(4).Value = hP.CData;
hold on

contour3(real(horizontal),real(vertical),contour,[cs cs],'black','linewidth',0.1);
xlabel(horizontal_label);   ylabel(vertical_label); zlabel(diagonal_label);
xlim(horizontal_lim);       ylim(vertical_lim);     zlim(diagonal_lim);
hold off

filename = append(in.prefix,"cross-section_",diagonal_label,"_.png");
saveas(gcf, filename);