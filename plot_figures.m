clear all
close all
clc

%% Affichage de l'évlution du TEB SANS SYNCHRO
t01 = load('teb_0.1Ts_sans_syn.mat');
t033 = load('teb_0.33Ts_sans_syn.mat');
t035 = load('teb_0.35Ts_sans_syn.mat');
t04 = load('teb_0.4s_sans_syn.mat');


c=0;
for i=1:length(t04.Pe)
    if t04.Pe(i) > 1e-6
        c = c+1;
    end
end

figure,
semilogy(t01.channel_params.EbN0dB, t01.ber, '-b^', ...
    t033.channel_params.EbN0dB, t033.ber, '-rd',... 
    t035.channel_params.EbN0dB, t035.ber, '-gv', ...
    t04.channel_params.EbN0dB, t04.ber, 'kh-', ...
    t04.channel_params.EbN0dB(1:c), t04.Pe(1:c), '-pm',...
    'LineWidth', 3) 

grid on;
title("Evolution du TEB en fonction de rapport Eb/N0 pour \tau \in [0.1T_s, 0.33T_s, 0.35T_s, 0.4T_s]");
xlabel('Eb/N0');
ylabel('BER')
legend('\tau = 0.1T_s','\tau = 0.33T_s','\tau = 0.35T_s','\tau = 0.4T_s',...
    'Théorique','Location','SouthWest');

%% Affichage de la perte de sensibilité
sensib = [0; 7.94-6.79; 28.76-6.79];
tau = [0; 0.1; 0.3];
figure,
plot(tau, sensib, '-b^', 'LineWidth', 3)
grid on 
title("Evolution de la perte de sensibilité en fonction de \tau")
xlabel('Perte de sensibilité')
ylabel('\tau*Ts')