function y=structure(wavelength);
periodicity_x=2.4;                  % period in x
periodicity_y=periodicity_x;        % period in y
diam=0.215/8;
structure_params = {
    periodicity_x,          0.075,   retindice_chen(wavelength,23.21),   retindice_chen(wavelength,23.21);   % 
%    diam*1,                 0.4/7,  ones(size(wavelength)),             retindice_chen(wavelength,4.802);
%    diam*2,                 0.4/7,  ones(size(wavelength)),             retindice_chen(wavelength,4.802);
%    diam*3,                 0.4/7,  ones(size(wavelength)),             retindice_chen(wavelength,4.802);
%    diam*4,                 0.4/7,  ones(size(wavelength)),             retindice_chen(wavelength,4.802);
%    diam*5,                 0.4/7,  ones(size(wavelength)),             retindice_chen(wavelength,4.802);
%    diam*6,                 0.4/7,  ones(size(wavelength)),             retindice_chen(wavelength,4.802);
%    diam*7,                 0.4/7,  ones(size(wavelength)),             retindice_chen(wavelength,4.802);
    periodicity_x,          0.04,   retindice_chen(wavelength,4.802),   0.00*ones(size(wavelength));        % 
    periodicity_x,          0.16,   retindice_chen(wavelength,4.707),   0.00*ones(size(wavelength));        % 
    periodicity_x,          0.14,   retindice_chen(wavelength,4.708),   0.00*ones(size(wavelength));        % 
    periodicity_x,          1.7,    retindice_chen(wavelength,4.707),   0.00*ones(size(wavelength));        % 
    periodicity_x,          0.04,   retindice_chen(wavelength,4.802),   0.00*ones(size(wavelength));        % 
%    periodicity_x*24/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*23/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*22/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*21/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*20/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*19/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*18/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*17/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*16/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*15/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*14/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*13/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*12/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*11/25,    0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*9/25,     0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*8/25,     0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*7/25,     0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*6/25,     0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*5/25,     0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*4/25,     0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*3/25,     0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*2/25,     0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
%    periodicity_x*1/25,     0.9/24, 1.58*ones(size(wavelength)),        retindice_chen(wavelength,4.802);   % 
    periodicity_x,          0.5,   retindice_chen(wavelength,1.72),    0.00*ones(size(wavelength));        % 
};
y = structure_params;