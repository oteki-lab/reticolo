function [x,y]=layers(in)
wavelength=linspace(in.lambdamin, in.lambdamax, in.npoints);    % range of wavelength

% diameter of each layer
% Thicknesses, from top to bottom   (0 si if no layer)
% Refraction indices (from top to bottom), can be a function of the wavelength
% params = [height, ni, nim, in.period_x, in.period_y]
structure_params = {
%    0.08,   retindice_chen(wavelength,23.21),   retindice_chen(wavelength,23.21),   in.period_x,           in.period_y;
	0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802),   in.diam_x*1/9,         in.diam_y*1/9;
	0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802),   in.diam_x*2/9,         in.diam_y*2/9;
	0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802),   in.diam_x*3/9,         in.diam_y*3/9;
	0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802),   in.diam_x*4/9,         in.diam_y*4/9;
	0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802),   in.diam_x*5/9,         in.diam_y*5/9;
	0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802),   in.diam_x*6/9,         in.diam_y*6/9;
	0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802),   in.diam_x*7/9,         in.diam_y*7/9;
	0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802),   in.diam_x*8/9,         in.diam_y*8/9;
    0.04,   retindice_chen(wavelength,4.802),   0.00*ones(size(wavelength)),        in.period_x,           in.period_y;
    0.16,   retindice_chen(wavelength,4.707),   0.00*ones(size(wavelength)),        in.period_x,           in.period_y;
    0.14,   retindice_chen(wavelength,4.708),   0.00*ones(size(wavelength)),        in.period_x,           in.period_y;
    1.7,    retindice_chen(wavelength,4.707),   0.00*ones(size(wavelength)),        in.period_x,           in.period_y;
    0.04,   retindice_chen(wavelength,4.802),   0.00*ones(size(wavelength)),        in.period_x,           in.period_y;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   in.period_x*11/12,     in.period_y*11/12;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   in.period_x*10/12,     in.period_y*10/12;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   in.period_x*9/12,      in.period_y*9/12;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   in.period_x*8/12,      in.period_y*8/12;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   in.period_x*7/12,      in.period_y*7/12;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   in.period_x*6/12,      in.period_y*6/12;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   in.period_x*5/12,      in.period_y*5/12;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   in.period_x*4/12,      in.period_y*4/12;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   in.period_x*3/12,      in.period_y*3/12;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   in.period_x*2/12,      in.period_y*2/12;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   in.period_x*1/12,      in.period_y*1/12;
    0.5,    retindice_chen(wavelength,1.72),    0.00*ones(size(wavelength)),        in.period_x,           in.period_y;
};
x = structure_params;

layers = {
    'SiNx ARC',         [1,2,3,4,5,6,7];
    'AlInP window',     [8];
    'GaAs emitter',     [9];
    'QD',               [10];
    'GaAs base',        [11];
    'AlInP BSF',        [12];
    'Ag mirror',        [14];
    'GaAs',             [9,11]
    'Active region',    [9,10,11]
    'Total',            [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
};
y = layers;
