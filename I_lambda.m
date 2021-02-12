load('Results\20210212174505328\no_1\no_1_period_700_diam_0wav900_1200_npoints51_Fourier0.mat', 'wavelength_A');

I = dir('Results\20210212174505328\no_1\I_map\*.mat');
[~, reindex] = sort( str2double( regexp( {I.name}, '\d+', 'match', 'once' )));
I = I(reindex);

I_mean = []; 
for i=1:length(I)
    load([I(i).folder,'\',I(i).name], 'E');
    M1 = mean(abs(E).^2, [1 2]);
    I_mean = [I_mean, M1];
end

size(I_mean)