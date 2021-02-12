load("Results\20210212172832465\no_1\I_map\I_map_2.mat")

figure
isosurface(x,y,z,abs(E).^2)
colorbar
colormap jet
[fe, ve, ce] = isocaps(x,y,z,abs(E).^2,10);
p2 = patch('Faces',fe, 'Vertices',ve, 'FaceVertexCData',ce);
p2.FaceColor = 'interp';
p2.EdgeColor = 'none';
%set(gca, 'clim',[0 14])
grid on
xlabel('x');
ylabel('y');
zlabel('z');
set(gca,'ZDir','reverse')

%saveas(gcf, "I_map.fig");