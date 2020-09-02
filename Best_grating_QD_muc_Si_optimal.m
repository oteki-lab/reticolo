clearvars
close all

E_G=1.424;
h=6.63e-34;
c=3e8;
q=1.60e-19;
k_B = 1.381e-23; %Boltzmann's constant
T_s=5800; %Sun Temperature

% sequence 1
% seq_num=1;
% period=[0.4:0.05:0.55];
% diam_perc=[0.4:0.1:0.6];
% height=[0.04:0.02:0.1];
% ref_index=1.3;
% npoints=101;
% Mx=10;
% wavmin=400;
% wavmax=900;
% h_plan=0.100;
%
% %sequence 2
% seq_num=2;
% period=[0.45];%:0.05:0.55];
% diam_perc=[0.6:0.05:0.80];
% height=[0.03:0.02:0.05];
% ref_index=1.3;
% npoints=72;
% Mx=10;
% wavmin=400;
% wavmax=900;
% h_plan=0.100;


% %sequence 3
% seq_num=3;
% period=[0.45:0.05:0.55];%:0.05:0.55];
% diam_perc=[0.70:0.05:0.75];
% height=[0.03:0.01:0.05];
% ref_index=1.3;
% npoints=72;
% Mx=10;
% wavmin=400;
% wavmax=900;
% h_plan=0.100;

% %sequence 4
% seq_num=4;
% period=[0.50];%:0.05:0.55];
% diam_perc=[0.75];
% height=[0.03:0.005:0.05];
% ref_index=1.3;
% npoints=72;
% Mx=10;
% wavmin=400;
% wavmax=900;
% h_plan=0.100;

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
thickness_GaAs=2000e-9;



load('solarspectrum') %spectra loaded in W.m^-2.nm^-1
wwavelength=wavelength(find(wavelength==wavmin):find(wavelength==wavmax));
AM1_5G=AM1_5G(find(wavelength==wavmin):find(wavelength==wavmax));
AAM1_5G=AM1_5G(1:4:801,1);
E=h*c./(wwavelength*1e-9)/q; %in eV flipud;�s��̏㉺�𔽓]
irr_nm=1e-4*AM1_5G; % in W.cm^-2.nm^-1 100�~100
irr_eV=irr_nm*1e9*h*c./(q.*E.^2); %from W.cm-2.nm^-1 to W.cm^-2.eV^-1
irr_kin_eV=max((1-E_G./E).*irr_eV,0); %the kinetic power in W.cm^-2.eV^-1
%back in nm
irr_kin_nm=irr_kin_eV*q.*E.^2/(1e9*h*c); %from W.cm-2.eV^-1 to W.cm^-2.nm^-1
P_kin_tot=nansum(irr_kin_eV(1:end-1).*(E(2:end)-E(1:end-1))); %Heat irradiance, or ideal thermalization intensity, assuming A=1, in W.cm^-2
% J_tot=nansum((irr_eV(find(E>=E_G,1,'first'):end-1)./E(find(E>=E_G,1,'first'):end-1)).*-(E(find(E>=E_G,1,'first')+1:end)-E(find(E>=E_G,1,'first'):end-1)));

J_tot=nansum((irr_eV(1:find(E<=E_G,1,'first')-2)./E(1:find(E<=E_G,1,'first')-2).*-(E(2:find(E<=E_G,1,'first')-1)-E(1:find(E<=E_G,1,'first')-2))));
J_total=nansum((irr_eV(1:end-1)./E(1:end-1).*-(E(2:end)-E(1:end-1))));

% E=flipud(h*c./(wwavelength*1e-9))/q; %in eV 
% irr_nm=1e-4*AM1_5G; % in W.cm^-2.nm^-1 100�~100
% irr_eV=flipud(irr_nm*1e9)*h*c./(q.*E.^2); %from W.cm-2.nm^-1 to W.cm^-2.eV^-1
% irr_kin_eV=max((1-E_G./E).*irr_eV,0); %the kinetic power in W.cm^-2.eV^-1
% %back in nm
% irr_kin_nm=flipud(irr_kin_eV*q.*E.^2/(1e9*h*c)); %from W.cm-2.eV^-1 to W.cm^-2.nm^-1
% P_kin_tot=nansum(irr_kin_eV(1:end-1).*(E(2:end)-E(1:end-1))); %Heat irradiance, or ideal thermalization intensity, assuming A=1, in W.cm^-2
% J_tot=nansum((irr_eV(find(E>=E_G,1,'first'):end-1)./E(find(E>=E_G,1,'first'):end-1)).*(E(find(E>=E_G,1,'first')+1:end)-E(find(E>=E_G,1,'first'):end-1)));
% J_tot=nansum(
%(irr_eV(find(E>=E_G,1,'first'):end-1)�@AM1.5[W.cm^-2.eV^-1]��E>E_G����Ō�̈��O�܂�
%./E(find(E>=E_G,1,'first'):end-1)
% )
%.*(
% E(find(E>=E_G,1,'first')+1:end)
% -E(find(E>=E_G,1,'first'):end-1)) (E2-E1,E3-E2,,
% );
% AM1.5/E*(E-E)

% load('GaAs_20nm_wav300_1000_flat_mirror.mat')
% A_flat=A;
% A_GaAs_flat=Abs_layer(4,:);
% ratio_to_GaAs=Abs_layer(4,:)./(Abs_layer(3,:)+Abs_layer(4,:)+Abs_layer(5,:)); % A rough idea to know which fraction of the absorbed power will end up in the GaAs
% 
% E_A_flat=fliplr(h*c./(q*lambda));
% ratio_to_GaAs_E=fliplr(ratio_to_GaAs);
% ratio_to_GaAs_eV=interp1(E_A_flat,ratio_to_GaAs_E,E,'makima');
% A_GaAs_flat_E=fliplr(A_GaAs_flat);
% A_GaAs_flat_eV=interp1(E_A_flat,A_GaAs_flat_E,E,'makima');
% P_kin_flat=nansum(A_GaAs_flat_eV(1:end-1).*irr_kin_eV(1:end-1).*(E(2:end)-E(1:end-1)));
% J_flat=nansum(A_GaAs_flat_eV(1:end-1).*(irr_eV(1:end-1)./E(1:end-1)).*-(E(2:end)-E(1:end-1)));
% clear A
%
% load('Reflection_ARC_20GaAs_wav400_850.mat')
% A_max=1-R;
% A_max_E=fliplr(A_max);
% A_max_eV=interp1(E_A_flat,A_max_E,E,'makima');
% P_kin_max=nansum(A_max_eV(1:end-1).*ratio_to_GaAs_eV(1:end-1).*irr_kin_eV(1:end-1).*(E(2:end)-E(1:end-1)));
%



ii=1;
jj=1;
kk=1;
ll=1;
number=1;
sequence{number}=[period(ii),diam_perc(jj),height(kk),ref_index(ll)];
% for number=1:length(period)*length(diam_perc)*length(height)*length(ref_index)
%     sequence{number}=[period(ii),diam_perc(jj),height(kk),ref_index(ll)];
%     if ll<length(ref_index)
%         ll=ll+1;
%     else
%         ll=1;
%         if kk<length(height)
%             kk=kk+1;
%         else
%             kk=1;
%             if jj<length(diam_perc)
%                 jj=jj+1;
%             else
%                 jj=1;
%                 ii=ii+1;
%             end
%         end
%     end
% end
wavelength=[];
text=['a.mat'];%�T���v���Ƃ���d�̍\���̃t�@�C��
load(text,'R0_TE_TE','A_tot','Abs','Abs_plots','wavelength');
load(text,'periodicity_x');

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
% %     xlabel('$\lambda \: \mathrm{(�m)}$','Interpreter','Latex')
%     ylabel('$A$','Interpreter','Latex')
%     box on
%     set(gca,'Fontsize',12)
%     set(gca,'XMinorTick','on','YMinorTick','on')
%     set(gcf,'color','w');

A{number}=1-R0_TE_TE;
A_GaAs{number}=Abs(3,:)+Abs(5,:)+Abs(4,:);%GaAs(QD�w)�̎w��
A_front_AlGaAs{number}=Abs(2,:);%���w�̎w��
A_back_AlGaAs{number}=Abs(6,:);%BSF�w�̎w��
A_Ag{number}=Abs(7,:);

% A_GaAs{number}=Abs(3,:)+Abs(5,:)+Abs(4,:);%GaAs(QD�w)�̎w��
% A_front_AlGaAs{number}=Abs(2,:);%���w�̎w��
% A_back_AlGaAs{number}=Abs(6,:)+Abs(7,:)+Abs(8,:)+Abs(9,:);%BSF�w�̎w��
% A_Ag{number}=Abs(10,:);
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

% m=size(wwavelength);
% m=m(1,1);
% multi=(m-1)/(npoints-1);
% irr_eV=irr_eV(1:multi:m,1);

P_kin(number)=nansum(A_GaAs_eV{number}(1:end-1).*irr_kin_eV(1:end-1).*(E(2:end)-E(1:end-1)));
J(number)=nansum(A_GaAs_eV{number}(1:end-1).*(irr_eV(1:end-1)./E(1:end-1)).*-(E(2:end)-E(1:end-1)));
J_EG=nansum((irr_eV(1:find(E<=E_G,1,'first')-2).*A_GaAs_eV{number}(1:find(E<=E_G,1,'first')-2))./E(1:find(E<=E_G,1,'first')-2).*-(E(2:find(E<=E_G,1,'first')-1)-E(1:find(E<=E_G,1,'first')-2)));

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
i=1;
% figure
% plot(wavelength,A{idx(i)},'linewidth',2)
% hold on
% plot(wavelength,A_total{idx(i)},'linewidth',2)
% plot(wavelength,A_front_AlGaAs{idx(i)},'linewidth',2)
% plot(wavelength,A_GaAs{idx(i)},'linewidth',2)
% plot(wavelength,A_back_AlGaAs{idx(i)},'linewidth',2)
% plot(wavelength,A_muc_Si{idx(i)},'linewidth',2)
% plot(wavelength,A_Ag{idx(i)},'linewidth',2)
% plot(lambda*1e6,A_GaAs_flat,'linewidth',2,'linestyle','-.','color','black')
% legend({'1-R','Total','Front AlGaAs','GaAs','Back AlGaAs','�c-Si','Ag','Flat mirror GaAs'})
% xlabel('Wavelength (�m)')
% ylabel('Absorptivity')
% box on
% set(gca,'Fontsize',12)
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gcf,'color','w');
% ylim([0 1])
% xlim([0.300 0.885])

colors=lines(7);

figure
plot(wavelength,A{idx(i)},'linewidth',2,'color',colors(2,:))
hold on
plot(wavelength,A_GaAs{idx(i)},'linewidth',2,'color',colors(1,:))
plot(wavelength,A_front_AlGaAs{idx(i)}+A_back_AlGaAs{idx(i)},'linewidth',2,'color',colors(5,:))
plot(wavelength,A_Ag{idx(i)},'linewidth',2,'color',[0.5 0.5 0.5])
% plot(lambda*1e6,A_GaAs_flat,'linewidth',2,'linestyle','-.','color',colors(1,:))
% legend({'1-R','Total','Front AlGaAs','GaAs','Back AlGaAs','�c-Si','Ag','Flat mirror GaAs'})
xlabel('Wavelength (��m)')
ylabel('Absorptivity')
box on
set(gca,'Fontsize',12)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w');
ylim([0 1])
xlim([0.400 1.2])

% load GaAs.dat
% lambda_models=GaAs(:,1)*1e-3;
% E_models=h*c./(q*lambda_models*1e-6);
% alpha_GaAs=4*pi*GaAs(:,3)./(lambda_models*1e-6);
% % A_single_pass=1-exp(-alpha_GaAs*thickness_GaAs);
% A_double_pass=1-exp(-2*alpha_GaAs*thickness_GaAs);
% A_double_pass_eV=interp1(E_models,A_double_pass,E);
% A_lamb=1-exp(-4*GaAs(:,2).^2.*alpha_GaAs*thickness_GaAs);
% A_lamb_eV=interp1(E_models,A_lamb,E);
% J_double_pass=nansum(A_double_pass_eV(1:end-1).*(irr_eV(1:end-1)./E(1:end-1)).*(E(2:end)-E(1:end-1)));
% J_lamb=nansum(A_lamb_eV(1:end-1).*(irr_eV(1:end-1)./E(1:end-1)).*(E(2:end)-E(1:end-1)));
% P_kin_double_pass=nansum(A_double_pass_eV(1:end-1).*irr_kin_eV(1:end-1).*(E(2:end)-E(1:end-1)));
% P_kin_lamb=nansum(A_lamb_eV(1:end-1).*irr_kin_eV(1:end-1).*(E(2:end)-E(1:end-1)));
% 
% 
% figure
% plot(wavelength,A_GaAs{idx(i)},'linewidth',2,'color',colors(1,:))
% hold on
% plot(lambda*1e6,A_GaAs_flat,'linewidth',2,'linestyle','-.','color',colors(1,:))
% % plot(lambda_models,A_single_pass,'linewidth',2,'color',colors(2,:))
% plot(lambda_models,A_double_pass,'--','linewidth',2,'color','black')
% plot(lambda_models,A_lamb,':','linewidth',3,'color','black')
% % legend({'1-R','Total','Front AlGaAs','GaAs','Back AlGaAs','�c-Si','Ag','Flat mirror GaAs'})
% xlabel('Wavelength (�m)')
% ylabel('Absorptivity')
% box on
% set(gca,'Fontsize',12)
% set(gca,'XMinorTick','on','YMinorTick','on')
% set(gcf,'color','w');
% ylim([0 1])
% xlim([0.300 0.885])
% text=['Cal at periodicity equal ',int2str(periodicity_x*1000),'um.mat'];
%     text='AbsRCWA_20nm_GaAs_check_Andrea.mat';
save(text);


