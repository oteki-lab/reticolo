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
    0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802),   diam_x*1/9,         diam_y*1/9;
	0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802),   diam_x*2/9,         diam_y*2/9;
	0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802),   diam_x*3/9,         diam_y*3/9;
	0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802),   diam_x*4/9,         diam_y*4/9;
	0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802),   diam_x*5/9,         diam_y*5/9;
	0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802),   diam_x*6/9,         diam_y*6/9;
	0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802),   diam_x*7/9,         diam_y*7/9;
	0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802),   diam_x*8/9,         diam_y*8/9;
    0.04,   retindice_chen(wavelength,4.802),   0.00*ones(size(wavelength)),        period_x,           period_y;
    0.16,   retindice_chen(wavelength,4.707),   0.00*ones(size(wavelength)),        period_x,           period_y;
    0.14,   retindice_chen(wavelength,4.708),   0.00*ones(size(wavelength)),        period_x,           period_y;
    1.7,    retindice_chen(wavelength,4.707),   0.00*ones(size(wavelength)),        period_x,           period_y;
    0.04,   retindice_chen(wavelength,4.802),   0.00*ones(size(wavelength)),        period_x,           period_y;
	0.9/12, retindice_chen(wavelength,1.73),    retindice_chen(wavelength,4.802),   period_x*11/12,     period_x*11/12;
	0.9/12, retindice_chen(wavelength,1.73),    retindice_chen(wavelength,4.802),   period_x*10/12,     period_x*10/12;
	0.9/12, retindice_chen(wavelength,1.73),    retindice_chen(wavelength,4.802),   period_x*9/12,      period_x*9/12;
	0.9/12, retindice_chen(wavelength,1.73),    retindice_chen(wavelength,4.802),   period_x*8/12,      period_x*8/12;
	0.9/12, retindice_chen(wavelength,1.73),    retindice_chen(wavelength,4.802),   period_x*7/12,      period_x*7/12;
	0.9/12, retindice_chen(wavelength,1.73),    retindice_chen(wavelength,4.802),   period_x*6/12,      period_x*6/12;
	0.9/12, retindice_chen(wavelength,1.73),    retindice_chen(wavelength,4.802),   period_x*5/12,      period_x*5/12;
	0.9/12, retindice_chen(wavelength,1.73),    retindice_chen(wavelength,4.802),   period_x*4/12,      period_x*4/12;
	0.9/12, retindice_chen(wavelength,1.73),    retindice_chen(wavelength,4.802),   period_x*3/12,      period_x*3/12;
	0.9/12, retindice_chen(wavelength,1.73),    retindice_chen(wavelength,4.802),   period_x*2/12,      period_x*2/12;
	0.9/12, retindice_chen(wavelength,1.73),    retindice_chen(wavelength,4.802),   period_x*1/12,      period_x*1/12;
	0.9/12, retindice_chen(wavelength,1.73),    retindice_chen(wavelength,4.802),   period_x*0.5/12,    period_x*0.5/12;
    0.5,    retindice_chen(wavelength,1.72),    0.00*ones(size(wavelength)),        period_x,           period_y;
};
x = structure_params;

layers = {
    'AlInP window',     [1,2,3,4,5,6,7,8,9];
    'GaAs emitter',     [10];
    'QD',               [11];
    'GaAs base',        [12];
    'AlInP BSF',        [13,14,15,16,17,18,19,20,21,22,23,24,25];
    'Ag mirror',        [26];
    'active region',    [10,11,12]
    'GaAs',             [10,12]
    'Total',            [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26]
};
y = layers;
