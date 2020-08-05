clear all
close all

%% Define parameters
h=6.626e-34; %Planck's constant
c=2.998e8; %Speed of light
q=1.602e-19; % Elementary charge

nb_lambda=400; %Number of wavelength
lambda=linspace(0.3,1,nb_lambda)*1e-6; %Wavelength vector
length_step=1e-9; %discretization step of the stack, for intensity calculation

cal_abs=1; % calculate the absorption and field intensity in each layer
cal_field=0; % calculate the field intensity in each layer
transform_E=0; % Transform the data into a function of energy

%% List of layers

n_0=1*ones(1,nb_lambda); %Incident medium
n_end=retindice_chen(lambda*1e6,2.7); %Last medium
% n_end=100i*ones(1,nb_lambda);

% list of indices of the stack

% n=[];
% n(1,:)=retindice_chen(lambda*1e6,23.33);
% n(2,:)=retindice_chen(lambda*1e6,23.11);
n(1,:)=retindice_chen(lambda*1e6,4.607);
% n(1,:)=3.5+0.01i*ones(1,nb_lambda);
n(2,:)=retindice_semicond(lambda*1e6,4.706);
n(3,:)=retindice_chen(lambda*1e6,4.607);
% n(5,:)=retindice_chen(lambda*1e6,4.702);
% n(6,:)=retindice_chen(lambda*1e6,4.621);
% n(4,:)=retindice_chen(lambda*1e6,23.11);
% n(4,:)=retindice_chen(lambda*1e6,2.7);

% list of thicknesses of the stack

% d(1)=70e-9;
% d(2)=40e-9;
d(1)=93e-9;
d(2)=200e-9;
d(3)=93e-9;
% d(4)=100e-9;
% d(5)=100e-9;
% d(6)=300e-9;

%% Initialization

nb_layers=size(n,1);
length_stack=sum(d);
nb_steps=round(d/length_step);
stack=linspace(0,length_stack,sum(nb_steps));
n_tot=[n_0; n; n_end]; %All indices, including media before and after the stack

Omega=cell(1,nb_lambda);
r=zeros(1,nb_lambda);
t=zeros(1,nb_lambda);
thetaz=zeros(nb_layers,sum(nb_steps),nb_lambda);
Delta=cell(nb_layers+1,nb_lambda);
Upsilon=cell(nb_layers,nb_lambda);
I_z=cell(1,nb_layers);
I_poynting=zeros(nb_layers+1,nb_lambda);
Omega_m_prime=cell(nb_layers+1,nb_lambda);
Omega_m=cell(nb_layers+1,nb_lambda);
% I_poynting_z=cell(1,nb_lambda);
% A_poynting_z=zeros(nb_layers,sum(nb_steps),nb_lambda);
% A_z=zeros(nb_layers,sum(nb_steps),nb_lambda);


%% Calculation of R and T

if isempty(n)
    r=(n_0-n_end)./(n_0+n_end);
    t=2*n_0./(n_0+n_end);
    R=abs(r).^2;
    T=(real(n_end)./real(n_0)).*abs(t).^2; % Transmittance
    A=1-R-T; % Absorption of the whole stack
    
else
    theta=2*pi*n.*d'./lambda; %Phase shift in a layer
    
    for j=1:nb_lambda
        for i=nb_layers:-1:1
            Delta{i+1,j}=(1/(2*n_tot(i+1,j)))*[n_tot(i+1,j)+n_tot(i+2,j) n_tot(i+1,j)-n_tot(i+2,j); n_tot(i+1,j)-n_tot(i+2,j) n_tot(i+1,j)+n_tot(i+2,j)];
            
            if i==nb_layers
                Omega_m{end,j}=(Delta{end,j}); % transfer matrix between layers 1 and i
            else
                Omega_m{i+1,j}=Delta{i+1,j}*Omega_m_prime{i+1,j};
            end
            Upsilon{i,j}=[exp(-1i*theta(i,j)) 0; 0 exp(1i*theta(i,j))];
            Omega_m_prime{i,j}=Upsilon{i,j}*Omega_m{i+1,j};
        end
        Delta{1,j}=(1/(2*n_tot(1,j)))*[n_tot(1,j)+n_tot(2,j) n_tot(1,j)-n_tot(2,j); n_tot(1,j)-n_tot(2,j) n_tot(1,j)+n_tot(2,j)];
        Omega_m{1,j}=Delta{1,j}*Omega_m_prime{1,j};
        Omega_m_prime{end,j}=[1 0;0 0];
        Omega{j}=Omega_m{1,j}; % Transfer matrix of the system
        r(j)=Omega{j}(2,1)/Omega{j}(1,1); % Reflection coefficient
        t(j)=1/Omega{j}(1,1); % Transmission coefficient
    end
    
    R=abs(r).^2; % Reflectance
    T=(real(n_end)./real(n_0)).*abs(t).^2; % Transmittance
    A=1-R-T; % Absorption of the whole stack
end
%% Calculation of absorptivity for each layer

if cal_abs==1
    for j=1:nb_lambda
        for i=1:nb_layers+1
            I_poynting(i,j)=abs(t(j))^2/real(n_0(j))*real(n_tot(i+1,j)*conj(Omega_m_prime{i,j}(1,1)+Omega_m_prime{i,j}(2,1))*(Omega_m_prime{i,j}(1,1)-Omega_m_prime{i,j}(2,1)));
        end
    end
    Abs_layer=I_poynting(1:nb_layers,:)-I_poynting(2:nb_layers+1,:);
end

%% Calculation of intensity of electric field with depth

if cal_field==1
    for j=1:nb_lambda
        for i=1:nb_layers
            for k=1:nb_steps(i)
                thetaz(i,k,j)=(k-1/2)*length_step*2*pi*n(i,j)/lambda(j);
                I_z{i}(k,j)=real(n(i,j))/real(n_0(j))*abs(t(j)*(Omega_m_prime{i,j}(1,1)*exp(1i*thetaz(i,k,j))+Omega_m_prime{i,j}(2,1)*exp(-1i*thetaz(i,k,j))))^2;
                A_z{i}(k,j)=(4*pi*imag(n(i,j))/lambda(j))*I_z{i}(k,j);
                LDOS_z{i}(k,j)=abs(1+Omega_m_prime{i,j}(2,1)*exp(-1i*thetaz(i,k,j))/(Omega_m_prime{i,j}(1,1)*exp(1i*thetaz(i,k,j))))^2;
            end
            for k=1:nb_steps(i)+1
                thetaz_poynting{i}(k,j)=(k-1)*length_step*2*pi*n(i,j)/lambda(j);
                I_poynting_z{i}(k,j)=abs(t(j))^2/real(n_0(j))*real(n_tot(i+1,j)*conj(Omega_m_prime{i,j}(1,1)*exp(1i*thetaz_poynting{i}(k,j))+Omega_m_prime{i,j}(2,1)*exp(-1i*thetaz_poynting{i}(k,j)))*(Omega_m_prime{i,j}(1,1)*exp(1i*thetaz_poynting{i}(k,j))-Omega_m_prime{i,j}(2,1)*exp(-1i*thetaz_poynting{i}(k,j))));
            end
            for k=1:nb_steps(i)
                A_poynting_z{i}(k,j)=(I_poynting_z{i}(k,j)-I_poynting_z{i}(k+1,j))/length_step;
            end
        end
    end
    
    I=I_z{1};
    A_local=A_z{1};
    LDOS=LDOS_z{1};
    A_poynting_local=A_poynting_z{1};
    I_poynting_local=I_poynting_z{1};
    for i=2:nb_layers
        I=cat(1,I,I_z{i});
        A_local=cat(1,A_local,A_z{i});
        LDOS=cat(1,LDOS,LDOS_z{i});
        A_poynting_local=cat(1,A_poynting_local,A_poynting_z{i});
        I_poynting_local=cat(1,I_poynting_local,I_poynting_z{i}(2:end,:));
    end
    
    I_poynting_local=(I_poynting_local(2:end,:)+I_poynting_local(1:end-1,:))/2;
    
end

%% Getting the energy vector

if transform_E
    E_A=fliplr(h*c./(q*lambda));
    A_GaAs_E=fliplr(Abs_layer(2,:));
    Abs_layer_E=fliplr(Abs_layer);
    A_total_E=fliplr(A+T);
    
    figure
    plot(E_A,A_total_E,'linewidth',3)
    xlim([1.3 3.5])
        save('Abs_TMM_AlGaAs_GaAs_HCSC_ELO_MBE2_Q3_2_2_wav_300_1000_habs_20_h_both_AlGaAs_82_Au_mirror.mat','E_A','A_GaAs_E','Abs_layer_E','A_total_E')
end


%% Plots

Icolors=lines(nb_layers);

figure
plot(lambda*1e6,1-R);
hold on
plot(lambda*1e6,T);
plot(lambda*1e6,A);
ylim([0 1])
xlabel('$\lambda \: \mathrm{(\mu m)}$','Interpreter','Latex')
set(gca,'Fontsize',12)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf,'color','w');
legend({'1-R','T','A'})
box on

if cal_abs==1
    figure
    plot(lambda*1e6,Abs_layer);
    xlabel('$\lambda \: \mathrm{(\mu m)}$','Interpreter','Latex')
    set(gca,'Fontsize',12)
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gcf,'color','w');
    legend({'layer 1','layer 2','layer 3'})
    box on
    ylim([0 1])
    
    
    %     figure
    %     hold on
    %     plot(lambda*1e6,I_poynting);
    %     xlabel('$\lambda \: \mathrm{(\mu m)}$','Interpreter','Latex')
    %     ylabel('$I_{Poynting}/I_0$','Interpreter','Latex')
    %     set(gca,'Fontsize',12)
    %     set(gca,'XMinorTick','on','YMinorTick','on')
    %     set(gcf,'color','w');
    %     legend({'Interface 1','Interface 2','Interface 3','Interface 4'})
    %     box on
    %     ylim([0 1])
end

if cal_field==1
    figure
    [C,h]=contourf(lambda*1e6,stack*1e6,I);
    hold all
    for i=1:nb_layers-1
        plot([lambda(1),lambda(end)]*1e6,[sum(d(1:i)) sum(d(1:i))]*1e6,'linewidth',1,'color','w')
    end
    xlabel('$\lambda \: \mathrm{(\mu m)}$','Interpreter','Latex')
    ylabel('depth $\: \mathrm{(\mu m)}$','Interpreter','Latex')
    set(gca,'Fontsize',12)
    set(gca,'YDir','reverse')
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gcf,'color','w');
    colormap('parula')
    h.LevelList=[0 1e-3 1e-2 1e-1 0.25 0.5 1 2 4 8];
    box on
    c=colorbar;
    c.Label.String = 'I_E/I_0';
    
    figure
    [C,h]=contourf(lambda*1e6,stack*1e6,I_poynting_local);
    hold all
    for i=1:nb_layers-1
        plot([lambda(1),lambda(end)]*1e6,[sum(d(1:i)) sum(d(1:i))]*1e6,'linewidth',1,'color','w')
    end
    xlabel('$\lambda \: \mathrm{(\mu m)}$','Interpreter','Latex')
    ylabel('depth $\: \mathrm{(\mu m)}$','Interpreter','Latex')
    set(gca,'Fontsize',12)
    set(gca,'YDir','reverse')
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gcf,'color','w');
    colormap('parula')
    h.LevelList=logspace(-15,0);
    set(gca,'colorscale','log')
    box on
    c=colorbar;
    c.Label.String = 'I_{Poynting}/I_0';
    
    
    %     figure
    %     plot(stack*1e6,A_local(:,100),'linewidth',3);
    %     hold on
    %     plot(stack*1e6,A_poynting_local(:,100),'o');
    %     xlabel('depth $\: \mathrm{(\mu m)}$','Interpreter','Latex')
    %     ylabel('Normalized absorbed power $\: \mathrm{(m^{-1})}$','Interpreter','Latex')
    %     set(gca,'Fontsize',12)
    %     set(gca,'XMinorTick','on','YMinorTick','on')
    %     set(gcf,'color','w');
    %     legend({'From electric field intensity','From Poynting vector'})
    %     box on
    %
    %     figure
    %     [C,h]=contourf(lambda*1e6,stack*1e6,LDOS);
    %     hold all
    %     for i=1:nb_layers-1
    %         plot([lambda(1),lambda(end)]*1e6,[sum(d(1:i)) sum(d(1:i))]*1e6,'linewidth',1,'color','w')
    %     end
    %     xlabel('$\lambda \: \mathrm{(\mu m)}$','Interpreter','Latex')
    %     ylabel('depth $\: \mathrm{(\mu m)}$','Interpreter','Latex')
    %     set(gca,'Fontsize',12)
    %     set(gca,'YDir','reverse')
    %     set(gca,'XMinorTick','on','YMinorTick','on')
    %     set(gcf,'color','w');
    %     colormap('parula')
    %     h.LevelList=logspace(-3,1);
    %     set(gca,'colorscale','log')
    %     box on
    %     c=colorbar;
    %     c.Label.String = 'LDOS';
    
end

