clear
close all
clc

%% Parametres des objets participant a l'emetteur
% Parametre du flux TS
% -------------------------------------------------------------------------
tx_vid_fname = 'tx_stream.ts';  % Fichier contenant le message ï¿½ transmettre
rx_vid_prefix = 'rx_stream';  % Fichier contenant le message ï¿½ transmettre

msg_oct_sz     = 188;
msg_bit_sz     = msg_oct_sz*8; % Taille de la payload des paquets en bits
pckt_per_frame = 8;
frame_oct_sz   = pckt_per_frame * msg_oct_sz;
frame_bit_sz   = 8*frame_oct_sz; % Une trame = 8 paquets
bool_store_rec_video = false;
% -------------------------------------------------------------------------

%% Création des structures de paramètres
waveform_params = configure_waveform(); % Les parametres de la mise en forme
channel_params  = configure_channel(0:10,0,0,1,0*(waveform_params.sim.Fse)); % LEs paramètres du canal

%% Création des objets
[mod_psk, demod_psk]           = build_mdm(waveform_params); % Construction des modems
dvb_scramble                   = build_dvb_scramble(); %Construction du scrambler
[awgn_channel, doppler, channel_delay] = build_channel(channel_params, waveform_params); % Blocs du canal
stat_erreur = comm.ErrorRate('ReceiveDelay', 1504*8, 'ComputationDelay',1504*8); % Calcul du nombre d'erreur et du BER

raised_cos_filter_trans = comm.RaisedCosineTransmitFilter('Shape', 'Square root',...
    'RolloffFactor', 0.35, ...
    'FilterSpanInSymbols', 16,...
    'OutputSamplesPerSymbol',waveform_params.sim.Fse); % 16 est 2*Tg qui correspond au span

raised_cos_filter_rec = comm.RaisedCosineReceiveFilter('Shape', 'Square root', ...
    'RolloffFactor', 0.35,...
    'FilterSpanInSymbols', 16,...
    'InputSamplesPerSymbol',waveform_params.sim.Fse, ...
    'DecimationFactor', waveform_params.sim.Fse);

vid = dsp.VariableIntegerDelay('MaximumDelay', 1504*8);

% Conversions octet <-> bits
o2b = OctToBit();
b2o = BitToOct();

% Lecture octet par octet du fichier vidéo d'entree
message_source      = BinaryFileReader('Filename', tx_vid_fname, 'SamplesPerFrame', msg_oct_sz*pckt_per_frame, 'DataType', 'uint8');
% Ecriture octet par octet du fichier vidéo de sortie
message_destination = BinaryFileWriter('DataType','uint8');

%%
ber = zeros(1,length(channel_params.EbN0dB));
Pe = qfunc(sqrt(2*channel_params.EbN0));

for i_snr = 1:length(channel_params.EbN0dB)
    if bool_store_rec_video
        message_destination.release;
        message_destination.Filename = [rx_vid_prefix, num2str(channel_params.EbN0dB(i_snr)),'dB.ts'];
    end
    
    awgn_channel.EbNo=channel_params.EbN0dB(i_snr);% Mise à jour du EbN0 pour le canal
    
    stat_erreur.reset; % reset du compteur d'erreur
    err_stat = [0 0 0];
    while (err_stat(2) < 100 && err_stat(3) < 1e6)
        message_source.reset;
        message_destination.reset;
        while(~message_source.isDone)
            %% Emetteur
            tx_oct     = step(message_source); % Lire une trame
            tx_scr_oct = bitxor(tx_oct,dvb_scramble); % scrambler
            tx_scr_bit = step(o2b,tx_scr_oct); % Octets -> Bits
            tx_sym     = step(mod_psk,  tx_scr_bit); % Modulation QPSK
            gt = step(raised_cos_filter_trans,tx_sym);
            %% Canal
            tx_sps_dpl = step(doppler, gt); % Simulation d'un effet Doppler
            rx_sps_del = step(channel_delay, tx_sps_dpl, channel_params.Delai); % Ajout d'un retard de propagation
            rx_sps     = step(awgn_channel,channel_params.Gain * rx_sps_del); % Ajout d'un bruit gaussien
            %% Recepteur
            ga = step(raised_cos_filter_rec, rx_sps);
            rx_scr_llr = step(demod_psk,ga);% Ce bloc nous renvoie des LLR (meilleur si on va interface avec du codage)
            rx_scr_bit = rx_scr_llr<0; % Bits
            rx_scr_bit_d = step(vid, rx_scr_bit, 1500*8);
            rx_scr_oct = step(b2o,rx_scr_bit_d); % Conversion en octet pour le scrambler
            % Attention à la synchro ici
            rx_oct     = bitxor(rx_scr_oct,dvb_scramble); % descrambler
            
            
            %% Compate des erreurs binaires
            tx_bit     = step(o2b,tx_oct);
            rx_bit     = step(o2b,rx_oct);
            err_stat   = step(stat_erreur, tx_bit, rx_bit);
            
            %% Destination
            if bool_store_rec_video
                step(message_destination, rx_oct_dec); % Ecriture du fichier
            end
        end
    end
    
    ber(i_snr) = err_stat(1);
    figure(1),semilogy(channel_params.EbN0dB,ber);
    drawnow
   
end
hold all
semilogy(channel_params.EbN0dB,Pe);

%% Tracé des constellations sans synchronisation pour EbN0dB = 100


