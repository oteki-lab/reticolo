function [x,y,z]=init_structure(in)
npoints=in.npoints;                    % 1 for only structure
lambdamin=in.lambdamin;                % min wavelength
lambdamax=in.lambdamax;                % max wavelength
wavelength=linspace(lambdamin,lambdamax,npoints);   % range of wavelength
period_x=in.period_x;                  % period in x
diam_x=in.diam_x;
Mx=in.Mx;

% diameter of each layer
% Thicknesses, from top to bottom   (0 si if no layer)
% Refraction indices (from top to bottom), can be a function of the wavelength
% params = [diameter_x, height, ni, nim]
structure_params = {
    period_x,          0.08,   retindice_chen(wavelength,23.21),   retindice_chen(wavelength,23.21);   % 
%    diam_x*1/9,             0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802);
%    diam_x*2/9,             0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802);
%    diam_x*3/9,             0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802);
%    diam_x*4/9,             0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802);
%    diam_x*5/9,             0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802);
%    diam_x*6/9,             0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802);
%    diam_x*7/9,             0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802);
%    diam_x*8/9,             0.4/8,  ones(size(wavelength)),             retindice_chen(wavelength,4.802);
    period_x,          0.04,   retindice_chen(wavelength,4.802),   0.00*ones(size(wavelength));        % 
    period_x,          0.16,   retindice_chen(wavelength,4.707),   0.00*ones(size(wavelength));        % 
    period_x,          0.14,   retindice_chen(wavelength,4.708),   0.00*ones(size(wavelength));        % 
    period_x,          1.7,    retindice_chen(wavelength,4.707),   0.00*ones(size(wavelength));        % 
    period_x,          0.04,   retindice_chen(wavelength,4.802),   0.00*ones(size(wavelength));        % 
%    period_x*24/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*23/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*22/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*21/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*20/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*19/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*18/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*17/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*16/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*15/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*14/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*13/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*12/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*11/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*9/25,     0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*8/25,     0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*7/25,     0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*6/25,     0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*5/25,     0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*4/25,     0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*3/25,     0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*2/25,     0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    period_x*1/25,     0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
    period_x,          0.5,   retindice_chen(wavelength,1.72),    0.00*ones(size(wavelength));        % 
};
x = structure_params;

layers = {
    'SiNx ARC',         [1];
    'AlInP window',     [2];
    'GaAs emitter',     [3];
    'QD',               [4];
    'GaAs base',        [5];
    'AlInP BSF',        [6];%,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30];
    'Ag mirror',        [7];
    'active region',    [3,4,5]
};
y = layers;

z = ['period_',int2str(period_x*1000),'_diam_',int2str(diam_x*1000),'wav',int2str(lambdamin*1000),'_',int2str(lambdamax*1000),'_npoints',int2str(npoints),'_Fourier',int2str(Mx)];
