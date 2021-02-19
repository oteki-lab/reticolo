function filename = save_cross_section(in, horizontal, vertical, indice, E, horizontal_label, vertical_label, diagonal_label, horizontal_lim, vertical_lim, zou)
%%%% The separation between the layers is in white (indice contains the position of all indices)
figure
pcolor(horizontal, vertical, abs(E).^2);
shading interp;
colormap(jet);
hold on
contour(real(horizontal), real(vertical), real(indice), 'black','linewidth',1)

colorbar;
caxis([0 Inf]);
xlabel(horizontal_label);   ylabel(vertical_label);
xlim(horizontal_lim);       ylim(vertical_lim);
hold off

filename = append(in.cs_dir, "\cross-section_",diagonal_label, "_", string(zou),"_.png");
saveas(gcf, filename);