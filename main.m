clear
close all
clc

%% Parametres des objets participant a l'emetteur
% Parametre du flux TS
% -------------------------------------------------------------------------
tx_vid_fname = 'tx_stream.ts';  % Fichier contenant le message � transmettre
rx_vid_prefix = 'rx_stream';  % Fichier contenant le message � transmettre
Fe = 1e6;
Ds = 250e3;
msg_oct_sz     = 188;
msg_bit_sz     = msg_oct_sz*8; % Taille de la payload des paquets en bits
pckt_per_frame = 8;
frame_oct_sz   = pckt_per_frame * msg_oct_sz;
frame_bit_sz   = 8*frame_oct_sz; % Une trame = 8 paquets
bool_store_rec_video = false;
mu =0;
Nfft = 1024;

% -------------------------------------------------------------------------

tau_init=0;
phi = 0; % en degr�
%% Cr�ation des structures de param�tres
waveform_params = configure_waveform(Fe, Ds); % Les parametres de la mise en forme
channel_params  = configure_channel(100,11021,60,1,0); % LEs param�tres du canal

%% Cr�ation des objets
[mod_psk, demod_psk]           = build_mdm(waveform_params); % Construction des modems
dvb_scramble                   = build_dvb_scramble(); %Construction du scrambler
[awgn_channel, doppler, channel_delay] = build_channel(channel_params, waveform_params); % Blocs du canal
stat_erreur = comm.ErrorRate('ReceiveDelay', 1504*8, 'ComputationDelay',1504*8*6); % Calcul du nombre d'erreur et du BER

raised_cos_filter_trans = comm.RaisedCosineTransmitFilter('Shape', 'Square root',...
    'RolloffFactor', 0.35, ...
    'FilterSpanInSymbols', 16,...
    'OutputSamplesPerSymbol',waveform_params.sim.Fse); % 16 est 2*Tg qui correspond au span

raised_cos_filter_rec = comm.RaisedCosineReceiveFilter('Shape', 'Square root', ...
    'RolloffFactor', 0.35,...
    'FilterSpanInSymbols', 16,...
    'InputSamplesPerSymbol',waveform_params.sim.Fse, ...
    'DecimationFactor',1);

synch = comm.CarrierSynchronizer(...
    'Modulation', 'QPSK', ...
    'SamplesPerSymbol', waveform_params.sim.Fse, ...
    'NormalizedLoopBandwidth', 0.005);

vid = dsp.VariableIntegerDelay('MaximumDelay', 1504*8);

% Conversions octet <-> bits
o2b = OctToBit();
b2o = BitToOct();

% Lecture octet par octet du fichier vid�o d'entree
message_source      = BinaryFileReader('Filename', tx_vid_fname, 'SamplesPerFrame', msg_oct_sz*pckt_per_frame, 'DataType', 'uint8');
% Ecriture octet par octet du fichier vid�o de sortie
message_destination = BinaryFileWriter('DataType','uint8');

%%
ber = zeros(1,length(channel_params.EbN0dB));
Pe = qfunc(sqrt(2*channel_params.EbN0));

for i_snr = 1:length(channel_params.EbN0dB)
    tau  = 0;
    
    if bool_store_rec_video
        message_destination.release;
        message_destination.Filename = [rx_vid_prefix, num2str(channel_params.EbN0dB(i_snr)),'dB.ts'];
    end
    
    
    awgn_channel.EbNo=channel_params.EbN0dB(i_snr);% Mise � jour du EbN0 pour le canal
    
    stat_erreur.reset; % reset du compteur d'erreur
    err_stat = [0 0 0];
    while (err_stat(2) < 100 && err_stat(3) < 1e6)
        message_source.reset;
        message_destination.reset;
        t0 = 0;
        while(~message_source.isDone)
            %% Emetteur
            tx_oct     = step(message_source); % Lire une trame
            tx_scr_oct = bitxor(tx_oct,dvb_scramble); % scrambler
            tx_scr_bit = step(o2b,tx_scr_oct); % Octets -> Bits
            tx_sym     = step(mod_psk,  tx_scr_bit); % Modulation QPSK
            Ns = length(tx_sym);
            gt = step(raised_cos_filter_trans,tx_sym);
            RL = abs(fftshift(fft(gt,1024))).^2;
            %% Canal
            tx_sps_dpl = step(doppler, gt); % Simulation d'un effet Doppler
            rx_sps_del = step(channel_delay, tx_sps_dpl, channel_params.Delai); % Ajout d'un retard de propagation
            rx_sps     = step(awgn_channel,channel_params.Gain * rx_sps_del); % Ajout d'un bruit gaussien
            
            %% Synchro fréquentielle
            [pxx,f] =pwelch(rx_sps.^4, hanning(Nfft), 0, Nfft, Fe, 'centered');
            [MAX, Ind] = max(pxx);
            t = t0+(1:length(rx_sps))/Fe;
            t0 = t(end) + 1/Fe;
            F_gross = f(Ind)/4;
            
            
            %% filtre adapté
            rl = step(raised_cos_filter_rec,rx_sps.*exp(-1i*2*pi*F_gross*t'));
            
            %% Synchro fine
            
            R = synch(rl);
            
            
            %% synchro temporelle
            rle = [zeros(waveform_params.sim.Fse,1); R; zeros(waveform_params.sim.Fse,1)];
            int_tau = floor(tau);
            te = waveform_params.sim.Fse + (1:waveform_params.sim.Fse:(1+(Ns-1)*waveform_params.sim.Fse)) + int_tau;
            frac_tau = tau - int_tau;
            r_n = rle(te)*(1-frac_tau) + rle(te + 1)*frac_tau;
            r_nd = 0.707*(sign(real(r_n)) + 1i * sign(imag(r_n)));
            err = r_n - r_nd;
            drl = 0.5*(1 - frac_tau)*(rle(te+1)-rle(te-1)) + 0.5*frac_tau*(rle(te+2)-rle(te));
            tau = tau - mu * real(err'*drl/Ns);
            
            
            %% Recepteur
            
            
            ejphi = 0;
            for ii = 1:8
                for i=1:4
                    ejphi = ejphi + r_n(16+i+(ii-1)*4*188)*(tx_sym(i+(ii-1)*4*188)');
                end
            end
            
            r_np = sqrt(32)*r_n/ejphi;
            
            rx_scr_llr = step(demod_psk,r_np);% Ce bloc nous renvoie des LLR (meilleur si on va interface avec du codage)
            rx_scr_bit = rx_scr_llr<0; % Bits
            rx_scr_bit_d = step(vid, rx_scr_bit, 1500*8);
            rx_scr_oct = step(b2o,rx_scr_bit_d); % Conversion en octet pour le scrambler
            % Attention � la synchro icilo
            rx_oct     = bitxor(rx_scr_oct,dvb_scramble); % descrambler
            
            %% Compare des erreurs binaires
            tx_bit     = step(o2b,tx_oct);
            rx_bit     = step(o2b,rx_oct);
            err_stat   = step(stat_erreur, tx_bit, rx_bit);
            
            
            %% Destination
            if bool_store_rec_video
                step(message_destination, rx_oct_dec); % Ecriture du fichier
            end
        end
    end
    
    % Minimisr fct periodique. Mais autant de minima que de periode lors de
    % la sgd. Rien ne dit qu'on convferge vers la m�me perdiode. On a une
    % rotation de pi/2 a chaque fois ou des fois.
    
    ber(i_snr) = err_stat(1);
    if ber(i_snr) > 1e-6
        figure(1),semilogy(channel_params.EbN0dB,ber);
        drawnow
    end
    
end
hold all
c=0;
for i=1:length(Pe)
    if Pe(i) > 1e-6
        c = c+1;
    end
end

semilogy(channel_params.EbN0dB(1:c),Pe(1:c));

figure;
plot(f,10*log(pxx/sum(pxx)))
hold on
plot(f,10*log(RL/sum(RL)))
xlabel('Frequency')
ylabel('Magnitude')
legend('Reponse en frequence du filtre de reception', 'Periodogramme du signal en sortie du canal')
title('Superposition periodogramme de Welch en sortie du canal et en sortie du filtre de reception')
% constellation r_n

%constellation a_n
scatterplot(tx_sym);
%% Trac� des constellations sans synchronisation pour EbN0dB = 100

%save('teb_phi50_sans_syn.mat','c', 'Pe','channel_params', 'ber');
