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
period=[0.48:0.01:0.52];%:0.05:0.55];
diam_perc=[0.71:0.01:0.755];
height=[0.04];
ref_index=1.3;
npoints=144;
Mx=10;
wavmin=400;
wavmax=900;
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


load('GaAs_20nm_wav400_850_flat_mirror.mat') % absorbance loaded in [-]
E_A_flat = fliplr(h*c./(q*lambda));                                             % from wavelength [nm] to energy [eV]
% absorbance
A_flat = A;                                                                     % total absorbance
A_GaAs_flat = Abs_layer(4,:);                                                   % absorbance in layer 4 (GaAs)
A_GaAs_flat_eV = interp1(E_A_flat, fliplr(A_GaAs_flat), E, 'makima');           % order from [nm] to [eV]
J_flat = nansum(A_GaAs_flat_eV(1:end-1).*(irr_eV(1:end-1)./E(1:end-1)).*(E(2:end)-E(1:end-1)));
% ratio
ratio_to_GaAs = Abs_layer(4,:)./(Abs_layer(3,:)+Abs_layer(4,:)+Abs_layer(5,:)); % A rough idea to know which fraction of the absorbed power will end up in the GaAs
ratio_to_GaAs_E = fliplr(ratio_to_GaAs);
ratio_to_GaAs_eV = interp1(E_A_flat, ratio_to_GaAs_E,E,'makima');
% back in nm
P_kin_flat = nansum(A_GaAs_flat_eV(1:end-1).*irr_kin_eV(1:end-1).*(E(2:end)-E(1:end-1)));
clear A


%load('Reflection_ARC_20GaAs_wav400_850.mat')
%A_max=1-R;
%A_max_E=fliplr(A_max);
%A_max_eV=interp1(E_A_flat,A_max_E,E,'makima');
%P_kin_max=nansum(A_max_eV(1:end-1).*ratio_to_GaAs_eV(1:end-1).*irr_kin_eV(1:end-1).*(E(2:end)-E(1:end-1)));


ii=1;
jj=1;
kk=1;
ll=1;
for number=1:length(period)*length(diam_perc)*length(height)*length(ref_index)
    sequence{number}=[period(ii),diam_perc(jj),height(kk),ref_index(ll)];
    if ll<length(ref_index)
        ll=ll+1;
    else
        ll=1;
        if kk<length(height)
            kk=kk+1;
        else
            kk=1;
            if jj<length(diam_perc)
                jj=jj+1;
            else
                jj=1;
                ii=ii+1;
            end
        end
    end
end

for number=1:length(sequence)
    text=['Optim muc-Si 20nm GaAs ' int2str(seq_num) '\AbsRCWA_20nm_GaAswav' num2str(wavmin) '_' num2str(wavmax) '_period_' int2str(sequence{number}(1)*1000) '_diam_' int2str(sequence{number}(1)*sequence{number}(2)*1000) 'h_struc' int2str(sequence{number}(3)*1000) '_h_plan_' num2str(h_plan*1000) '_n_plan_' num2str(sequence{number}(4)) '_Fourier' num2str(Mx) '_nbpoints' num2str(npoints) '_muc_Si_grid.mat'];
    load(text,'R0_TE_TE','A_tot','Abs','Abs_plots','wavelength');
    
%     figure
%     plot(wavelength,1-R0_TE_TE,'linewidth',2)
%     hold on
%     plot(wavelength,A_tot,'linewidth',2)
%     plot(wavelength,Abs(3,:),'linewidth',2)
%     plot(wavelength,Abs(4,:),'linewidth',2)
%     plot(wavelength,Abs(5,:),'linewidth',2)
%     plot(wavelength,Abs_plots(7,:)+Abs(8,:),'linewidth',2)
%     plot(lambda*1e6,A_GaAs_flat,'linewidth',2)
%     legend({'1-R','Total abs','Front AlGaAs abs','GaAs abs','Back AlGaAs abs','Ag abs'})
% %     xlim([-inf h*c*1e9./(E_G*q)])
% %     xlabel('$\lambda \: \mathrm{(µm)}$','Interpreter','Latex')
%     ylabel('$A$','Interpreter','Latex')
%     box on
%     set(gca,'Fontsize',12)
%     set(gca,'XMinorTick','on','YMinorTick','on')
%     set(gcf,'color','w');

    A{number}=1-R0_TE_TE;
    A_GaAs{number}=Abs(4,:);
    A_front_AlGaAs{number}=Abs(3,:);
    A_back_AlGaAs{number}=Abs(5,:);
    A_muc_Si{number}=Abs(6,:);
    A_Ag{number}=Abs(8,:);
    A_total{number}=A_tot;
    
    E_A=fliplr(h*c./(q*wavelength*1e-6));
    A_GaAs_E=fliplr(A_GaAs{number});
    A_front_AlGaAs_E=fliplr(A_front_AlGaAs{number});
    A_back_AlGaAs_E=fliplr(A_back_AlGaAs{number});
    A_Ag_E=fliplr(A_Ag{number});
    A_total_E=fliplr(A_total{number});
    
    A_GaAs_eV{number}=interp1(E_A,A_GaAs_E,E,'makima');
    A_front_AlGaAs_eV{number}=interp1(E_A,A_front_AlGaAs_E,E,'makima');
    A_back_AlGaAs_eV{number}=interp1(E_A,A_back_AlGaAs_E,E,'makima');
    A_Ag_eV{number}=interp1(E_A,A_Ag_E,E,'makima');
    A_total_eV{number}=interp1(E_A,A_total_E,E,'makima');
    
    P_kin(number)=nansum(A_GaAs_eV{number}(1:end-1).*irr_kin_eV(1:end-1).*(E(2:end)-E(1:end-1)));
    J(number)=nansum(A_GaAs_eV{number}(1:end-1).*(irr_eV(1:end-1)./E(1:end-1)).*(E(2:end)-E(1:end-1)));

end

% % for number=1:length(sequence)
% %     FOM(number)=nansum(A_GaAs_eV{number}(1:end-1).*irr_eV(1:end-1).*(E(2:end)-E(1:end-1)));
% % end
[val,idx]=sort(P_kin,'descend');
for i=1:length(sequence)
    seq_ordered(i,1)=sequence{idx(i)}(1);
    seq_ordered(i,2)=sequence{idx(i)}(2);
    seq_ordered(i,3)=sequence{idx(i)}(3);
    seq_ordered(i,4)=sequence{idx(i)}(4);
    seq_ordered(i,5)=P_kin(idx(i));
    seq_ordered(i,6)=J(idx(i));
end
% % sequence{idx}
% % 
% 
for i=1:length(sequence)
    figure
    
    plot(wavelength,A{idx(i)},'linewidth',2)
    hold on
    plot(wavelength,A_total{idx(i)},'linewidth',2)
    plot(wavelength,A_front_AlGaAs{idx(i)},'linewidth',2)
    plot(wavelength,A_GaAs{idx(i)},'linewidth',2)
    plot(wavelength,A_back_AlGaAs{idx(i)},'linewidth',2)
    plot(wavelength,A_muc_Si{idx(i)},'linewidth',2)
    plot(wavelength,A_Ag{idx(i)},'linewidth',2)
    plot(lambda*1e6,A_GaAs_flat,'linewidth',2)
    legend({'1-R','Total','Front AlGaAs','GaAs','Back AlGaAs','µc-Si','Ag','Flat mirror GaAs'})
    ylabel('$A$','Interpreter','Latex')
    box on
    set(gca,'Fontsize',12)
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gcf,'color','w');
    ylim([0 1])
end

% figure
% plot(wwavelength,flipud(A_GaAs_eV{idx(1)}))
% hold on
% plot(wwavelength,flipud(A_flat_eV))
% % plot(E,A_max_eV)
% % plot(E,A_flat_eV)
% ylim([0 1])