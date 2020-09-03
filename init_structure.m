function [x,y]=init_structure(in)
npoints=in.npoints;                    % 1 for only structure
lambdamin=in.lambdamin;                % min wavelength
lambdamax=in.lambdamax;                % max wavelength
wavelength=linspace(lambdamin,lambdamax,npoints);   % range of wavelength
period_x=in.period_x;                  % period in x
diam_perc=in.diam_perc;
height=in.height;
diam_x=in.diam_x;
Mx=in.Mx;

% diameter of each layer
% Thicknesses, from top to bottom   (0 si if no layer)
% Refraction indices (from top to bottom), can be a function of the wavelength
% params = [diameter_x, height, ni, nim]
structure_params = {
    period_x,          0.08,   retindice_chen(wavelength,23.21),   retindice_chen(wavelength,23.21);   % 
    period_x,          0.04,   retindice_chen(wavelength,4.802),   0.00*ones(size(wavelength));        % 
    period_x,          0.16,   retindice_chen(wavelength,4.707),   0.00*ones(size(wavelength));        % 
    period_x,          0.14,   retindice_chen(wavelength,4.708),   0.00*ones(size(wavelength));        % 
    period_x,          1.70,   retindice_chen(wavelength,4.707),   0.00*ones(size(wavelength));        % 
    period_x,          0.04,   retindice_chen(wavelength,4.802),   0.00*ones(size(wavelength));        % 
    period_x*diam_perc,height, retindice_chen(wavelength,1.73),    retindice_chen(wavelength,4.802);   % 
    period_x,          0.01,   retindice_chen(wavelength,1.72),    0.00*ones(size(wavelength));        % 
};
x = structure_params;

layers = {
    'SiNx ARC',         [1];
    'AlInP window',     [2];
    'GaAs emitter',     [3];
    'QD',               [4];
    'GaAs base',        [5];
    'AlInP BSF',        [6,7];
    'Ag mirror',        [8];
    'active region',    [3,4,5]
};
y = layers;
