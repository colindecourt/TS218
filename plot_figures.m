clear all
close all
clc

%% Affichage de l'évlution du TEB SANS SYNCHRO TEMPORELLE
t01 = load('teb_0.1Ts_sans_syn.mat');
t033 = load('teb_0.33Ts_sans_syn.mat');
t035 = load('teb_0.35Ts_sans_syn.mat');
t04 = load('teb_0.4s_sans_syn.mat');


figure,
semilogy(t01.channel_params.EbN0dB, t01.ber, '-b^', ...
    t033.channel_params.EbN0dB, t033.ber, '-rd',...
    t035.channel_params.EbN0dB, t035.ber, '-gv', ...
    t04.channel_params.EbN0dB, t04.ber, 'kh-', ...
    t04.channel_params.EbN0dB(1:t04.c), t04.Pe(1:t04.c), '-pm',...
    'LineWidth', 3)

grid on;
title({'Evolution du TEB en fonction de rapport Eb/N0 pour \tau \in [0.1T_s, 0.33T_s, 0.35T_s, 0.4T_s]' ...
    '(sans synchronisation)'});
xlabel('Eb/N0');
ylabel('BER')
legend('\tau = 0.1T_s','\tau = 0.33T_s','\tau = 0.35T_s','\tau = 0.4T_s',...
    'Théorique','Location','SouthWest');

% Affichage de la perte de sensibilité
sensib = [0; 7.94-6.79; 28.76-6.79];
tau = [0; 0.1; 0.3];
figure,
plot(tau, sensib, '-b^', 'LineWidth', 3)
grid on
title({'Evolution de la perte de sensibilité en fonction de \tau' ...
    '(sans synchronisation)'})
ylabel('Perte de sensibilité')
xlabel('\tau*Ts')


%% Affichage de l'évlution du TEB AVEC SYNCHRO TEMPORELLE
t01p = load('teb_0.1ts_avec_synT_EbN0.mat');
t033p = load('teb_0.33ts_avec_synT_EbN0.mat');
t035p = load('teb_0.35ts_avec_synT_EbN0.mat');
t04p = load('teb_0.4ts_avec_synT_EbN0.mat');


figure,
semilogy(t01p.channel_params.EbN0dB, t01p.ber, '-b^', ...
    t033p.channel_params.EbN0dB, t033p.ber, '-rd',...
    t035p.channel_params.EbN0dB, t035p.ber, '-gv', ...
    t04p.channel_params.EbN0dB, t04p.ber, 'kh-', ...
    t04p.channel_params.EbN0dB(1:t04p.c), t04.Pe(1:t04p.c), '-pm',...
    'LineWidth', 3)

grid on;
title({'Evolution du TEB en fonction de rapport Eb/N0 pour \tau \in [0.1T_s, 0.33T_s, 0.35T_s, 0.4T_s]' ...
    '(avec synchronisation temporelle)'});
xlabel('Eb/N0');
ylabel('BER')
legend('\tau = 0.1T_s','\tau = 0.33T_s','\tau = 0.35T_s','\tau = 0.4T_s',...
    'Théorique','Location','SouthWest');

% Affichage de la perte de sensibilité
sensib = [0; 7.13-6.79; 6.97-6.79; 7.13-6.79; 7.13-6.79];
tau = [0; 0.1; 0.33; 0.35; 0.4];
figure,
plot(tau, sensib, '-b^', 'LineWidth', 3)
grid on
title({'Evolution de la perte de sensibilité en fonction de \tau' 'avec synchronisation'})
ylabel('Perte de sensibilité')
xlabel('\tau*Ts')


%% Affichage de l'évlution du TEB SANS SYNCHRO FREQUENTIELLE
phi10 = load('teb_phi10_sans_syn.mat');
phi30 = load('teb_phi30_sans_syn.mat');
phi50 = load('teb_phi50_sans_syn.mat');


figure,
semilogy(phi10.channel_params.EbN0dB, phi10.ber, '-b^', ...
    phi30.channel_params.EbN0dB, phi30.ber, '-rd',...
    phi50.channel_params.EbN0dB, phi50.ber, '-gv', ...
    phi50.channel_params.EbN0dB(1:phi50.c), phi50.Pe(1:phi50.c), '-pm',...
    'LineWidth', 3)

grid on;
title("Evolution du TEB en fonction de rapport Eb/N0 pour \phi \in [10°, 30°, 50°]");
xlabel('Eb/N0');
ylabel('BER')
legend('\phi = 10°','\phi = 30°','\phi = 50°',...
    'Théorique','Location','SouthWest');

% Affichage de la perte de sensibilité
% sensib = [0; 7.94-6.79; 28.76-6.79];
% tau = [0; 0.1; 0.3];
% figure,
% plot(tau, sensib, '-b^', 'LineWidth', 3)
% grid on
% title("Evolution de la perte de sensibilité en fonction de \tau")
% xlabel('Perte de sensibilité')
% ylabel('\tau*Ts')