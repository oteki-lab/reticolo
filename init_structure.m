function [x,y]=init_structure(in)
npoints=in.npoints;                                 % 1 for only structure
lambdamin=in.lambdamin;                             % min wavelength
lambdamax=in.lambdamax;                             % max wavelength
wavelength=linspace(lambdamin,lambdamax,npoints);   % range of wavelength
period_x=in.period_x;                               % period in x
period_y=in.period_y;                               % period in y
diam_x=in.diam_x;
diam_y=in.diam_y;
Mx=in.Mx;

% diameter of each layer
% Thicknesses, from top to bottom   (0 si if no layer)
% Refraction indices (from top to bottom), can be a function of the wavelength
% params = [height, ni, nim, period_x, period_y]
structure_params = {
    0.08,   retindice_chen(wavelength,23.21),   retindice_chen(wavelength,23.21),   period_x,           period_y;
    0.04,   retindice_chen(wavelength,4.802),   0.00*ones(size(wavelength)),        period_x,           period_y;
    0.16,   retindice_chen(wavelength,4.707),   0.00*ones(size(wavelength)),        period_x,           period_y;
    0.14,   retindice_chen(wavelength,4.708),   0.00*ones(size(wavelength)),        period_x,           period_y;
    1.7,    retindice_chen(wavelength,4.707),   0.00*ones(size(wavelength)),        period_x,           period_y;
    0.04,   retindice_chen(wavelength,4.802),   0.00*ones(size(wavelength)),        period_x,           period_y;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   period_x*11/12,     period_y*11/12;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   period_x*10/12,     period_y*10/12;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   period_x*9/12,      period_y*9/12;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   period_x*8/12,      period_y*8/12;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   period_x*7/12,      period_y*7/12;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   period_x*6/12,      period_y*6/12;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   period_x*5/12,      period_y*5/12;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   period_x*4/12,      period_y*4/12;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   period_x*3/12,      period_y*3/12;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   period_x*2/12,      period_y*2/12;
	0.9/12, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802),   period_x*1/12,      period_y*1/12;
    0.5,    retindice_chen(wavelength,1.72),    0.00*ones(size(wavelength)),        period_x,           period_y;
};
x = structure_params;

layers = {
    'SiNx ARC',         [1];
    'AlInP window',     [2];
    'GaAs emitter',     [3];
    'QD',               [4];
    'GaAs base',        [5];
    'AlInP BSF',        [6,7,8,9,10,11,12,13,14,15,16,17];
    'Ag mirror',        [18];
    'GaAs',             [3,5]
    'Active region',    [3,4,5]
    'Total',            [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
};
y = layers;
