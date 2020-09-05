clearvars
close all

E_G=1.424;          % bandgap
h=6.63e-34;         % planck constant
c=3e8;              % speed of light
q=1.60e-19;         % elementary charge
k_B = 1.381e-23;    % Boltzmann's constant
T_s=5800;           % Sun Temperature

%sequence 5
seq_num=5;
period=[2.40];
diam_perc=[0.215];
height=[2.0];
ref_index=1.3;
npoints=201;
Mx=15;
wavmin=400;
wavmax=1200;
h_plan=0.100;

load('Solarspectrum.mat') % spectra loaded in [W.m^-2.nm^-1]
% wavelength
wwavelength = wavelength(find(wavelength==wavmin):find(wavelength==wavmax));    % wavelength from wavmin to wavmax [nm]
E = flipud(h*c./(wwavelength*1e-9))/q;                                          % from wavelength [nm] to energy [eV]
% solar spectrum
AM1_5G = AM1_5G(find(wavelength==wavmin):find(wavelength==wavmax));             % solar spectrum from wavmin to wavmax [W.m^-2.nm^-1]
irr_nm = 1e-4*AM1_5G;                                                           % solar spectrum from [W.m^-2.nm^-1] to [W.cm^-2.nm^-1]
irr_eV = flipud(irr_nm*1e9)*h*c./(q*E.^2);                                      % solar spectrum from [W.cm-2.nm^-1] to [W.cm^-2.eV^-1]
irr_kin_eV = max((1-E_G./E).*irr_eV,0);                                         % the kinetic power in [W.cm^-2.eV^-1]
J_tot = nansum((irr_eV(1:end-1)./E(1:end-1)).*(E(2:end)-E(1:end-1)));           % Total current that can be converted from sunlight
% back in nm
irr_kin_nm = flipud(irr_kin_eV*q.*E.^2/(1e9*h*c));                              % the kinetic power from [W.cm-2.eV^-1] to [W.cm^-2.nm^-1]
P_kin_tot = nansum(irr_kin_eV(1:end-1).*(E(2:end)-E(1:end-1)));                 % Heat irradiance, or ideal thermalization intensity, assuming A=1, in W.cm^-2

number=1;

dir = 'results';
file = 'period_2400_diam_27wav400_1200_nbpoints101_Fourier0.mat';
text = [dir,'\',file];

% read aborption rate
load(text,'A_tot','Abs','Abs_array','layers','wavelength');

% from wavelength [nm] to energy [eV]
E_A=fliplr(h*c./(q*wavelength*1e-6));
for index = 1:length(Abs_array)
    Abs_array_E{index} = fliplr(Abs_array{index});
end

% 1-D data interpolation
E=fliplr(E);
for index = 1:length(Abs_array)
    Abs_array_eV{index} = interp1(E_A,Abs_array_E{index},E,'makima');
end

% absorption table
header = ["wavelength","energy"];
data = [flipud(wwavelength),E];
for index = 1:length(Abs_array)
    header= horzcat(header, layers{index});
    data = horzcat(data, Abs_array_eV{index});
end
Abs_table = vertcat(header,data);
Abs_layers = containers.Map(layers(:,1), Abs_array_eV);
%disp(Abs_layers(char(layers(1,1))));

% calculate J
header = ["No"];
J_array = [number];
for index = 1:length(Abs_array)
    header= horzcat(header, layers{index});
    J_layer=nansum(Abs_array_eV{index}(1:end-1).*(irr_eV(1:end-1)./E(1:end-1)).*(E(2:end)-E(1:end-1)));
    J_array = horzcat(J_array, J_layer);
end
J_table = vertcat(header,J_array);

A_active = Abs_layers(char('Active region'));
P_kin(number)=nansum(A_active(1:end-1).*irr_kin_eV(1:end-1).*(E(2:end)-E(1:end-1)));
J(number)=nansum(A_active(1:end-1).*(irr_eV(1:end-1)./E(1:end-1)).*(E(2:end)-E(1:end-1)));

save([dir,'\J_',file]);
