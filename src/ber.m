%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ATTENTION ATTENTION : Pour toute fonction de Matlab vous avez une aide en
%ligne accessible avec la commande : help nomdelafonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
tic

%%%%%%%%%%%%%%%%%%% Variables à utiliser %%%%%%%%%%%%%%%%%%%%%

nsb = 30		%Nombre de bits

fb = 1e5		%Fréquence des bits
Tb = 1/fb		%Période des bits
Eb = 1			%Energie par bit

T = 2*Tb		%Période des symboles
Nsymb = nsb/2   %Nombre de symboles

fc = 1e6		%Fréquence de la porteuse (carrier)

Ts = 1/(10*fc)	%Période d'échantillonnage : définit le pas de simulation
fs = 1/Ts		%Fréquence d'échantillonnage

npb = round(Tb/Ts)		%Nombre de points par bit

t = [Ts:Ts:nsb*Tb];     %Vecteur représentant l'axe du temps


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% EMETTEUR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Calcul d'une suite aleatoire de bits %%%%%%%%%%%%%%
seqbitsinit = 2*round(rand (1,nsb))-1;	%Suite de bits en NRZ (-1,1)

seqbits = sqrt(Eb) .* seqbitsinit;	%Mise en valeur reelle pour
									%l'energie de la suite de bits 

%-------------Demultiplexeur
path_I_tx = seqbits(1:2:length(seqbits)-1);	%Sous-suite de bits impairs
path_Q_tx = seqbits(2:2:length(seqbits));		%Sous-suite de bits pairs

qpsk_I_tx=path_I_tx(ones(2*npb,1),:);
qpsk_I_tx=qpsk_I_tx(:).';
qpsk_Q_tx=path_Q_tx(ones(2*npb,1),:);
qpsk_Q_tx=qpsk_Q_tx(:).';

%-------------Liberation de memoire (ces variables sont effaces)
%clear seqbits path_I path_Q;

%-------------Generation du signal QPSK
qpsk = 2*sqrt(Eb/T).* cos(2*pi*fc.*t - qpsk_Q_tx.*(pi/2 - qpsk_I_tx*pi/4));

%sim('QPSK_Arm', T);
%-------------Liberation de memoire (ces variables sont effaces)
%clear qpsk_I_tx qpsk_Q_tx;

fprintf(1,'\nTemps de génération de %.0f bits : %.6fs\n',nsb,toc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% CANAL DE TRANSMISSION %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RECEPTEUR %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
% Multiplieurs : la multiplication de vecteurs se fait avec l'opérateur .*
qpsk_I_rx = qpsk .* sqrt(2/T).*cos(2*pi*fc.*t);	%sortie du multiplieur du chemin I
qpsk_Q_rx = qpsk .* sqrt(2/T).*sin(2*pi*fc.*t);	%sortie du multiplieur du chemin Q

%-------------Liberation de memoire (ces variables sont effaces)
%clear qpsk;

% Integrateurs
rI = (reshape(qpsk_I_rx(1,:),length(t)/Nsymb,Nsymb))';
rintI = (Ts*npb)*cumsum(rI/npb,2);
intI = (reshape(rintI',1,length(t)));

rQ = (reshape(qpsk_Q_rx(1,:),length(t)/Nsymb,Nsymb))';
rintQ = (Ts*npb)*cumsum(rQ/npb,2);
intQ = (reshape(rintQ',1,length(t)));

fprintf(1,'\nTemps integration de %.0f bits : %.6fs\n',nsb,toc);
tic
%-------------Liberation de memoire (ces variables sont effaces)
%clear rI rintI rQ rintQ;

% Echantillonneurs
path_I_rx = intI(round(T/Ts-1):round(T/Ts):round(Nsymb*T/Ts-1));
path_Q_rx = intQ(round(T/Ts-1):round(T/Ts):round(Nsymb*T/Ts-1));

% Multiplexeur
nrz_I_rx = sign(path_I_rx);
nrz_Q_rx = sign(path_Q_rx);
cn = cat(1,nrz_I_rx,nrz_Q_rx);
resbits = reshape(cn,1,2*Nsymb);
	
if (resbits == seqbitsinit)
    fprintf(1,'ok')
else
    fprintf(1,'ko')
end
fprintf(1,'\nTemps de demodulation de %.0f bits : %.6fs\n',nsb,toc);

% Injection du bruit
 SNR = 22;
 NR_I = awgn(path_I_rx,SNR,0);
 NR_Q = awgn(path_Q_rx,SNR,0);
 
%Affichage

%scatter(path_I_rx,path_Q_rx,'d'),grid
%scatter(NR_I,NR_Q),grid;

plot_lims = [-2 2];
plot(NR_I,NR_Q, 'd');
hold on
plot(path_I_rx,path_Q_rx,'d');
xlim(plot_lims);
ylim(plot_lims);
title('QPSK constellation noise');
xlabel('real part');
ylabel('imaginary part');

%EVM
%https://fr.mathworks.com/help/comm/ref/evmmeasurement.html
% Ierr = path_I_rx - NR_I;
% Qerr = path_Q_rx - NR_Q;
% EVMref = (1/length(path_I_rx))*sum(power(path_I_rx,2)+power(path_Q_rx,2));
% EVMmes = (1/length(NR_I))*sum(power(Ierr,2)+power(Qerr,2));
% EVMpercentage = sqrt(EVMmes/EVMref)*100 %rms
% matlab RF toolbox
% evm = comm.EVM;
% rmsEVM1 = evm(nrz_I_rx,NR_I);
% release(evm);


%Affichage de l'affluence de SNR sur l'EVM

% Injection du bruit
y = 0;
for i = 1.0:1:200
    SNR = i;
    NR_I = awgn(path_I_rx,SNR,0);
    NR_Q = awgn(path_Q_rx,SNR,0)
    Ierr = path_I_rx - NR_I;
    Qerr = path_Q_rx - NR_Q;
    EVMref = (1/length(path_I_rx))*sum(power(path_I_rx,2)+power(path_Q_rx,2));
    EVMmes = (1/length(NR_I))*sum(power(Ierr,2)+power(Qerr,2));
    EVMpercentage = sqrt(EVMmes/EVMref)*100
    y(1,i) = EVMpercentage;
end

figure 
plot(y);


%%%%%%%%%%%%%%%%%%%%%%%%%%%% AFFICHAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AFFICHAGE : vous pouvez utiliser plot directement pour l'affichage d'une
% courbe dans une fenetre. Sinon vous pouvez l'utiliser avec subplot pour
% diviser la fenetre en plusieurs partie et afficher plusieurs courbes sans
% les superposer. HELP plot ou subplot dans Matlab pour plus de details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


