%ATTENTION ATTENTION : Pour toute fonction de Matlab vous avez une aide en
%ligne accessible avec la commande : help nomdelafonction
clear all;
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%% TRANSMISSION QPSK %%%%%%%%%%%%%%%%%%%%%%%%%

tic

%%%%%%%%%%%%%%%%%%% Variables � utiliser %%%%%%%%%%%%%%%%%%%%%

nsb = 100;		%Nombre de bits

fb = 1e5		%Fr�quence des bits
Tb = 1/fb		%P�riode des bits
Eb = 1			%Energie par bit

T = 2*Tb		%P�riode des symboles
Nsymb = nsb/2   %Nombre de symbolesseq

fc = 1e6		%Fr�quence de la porteuse (carrier)

Ts = 1/(100*fc)	%P�riode d'�chantillonnage : d�finit le pas de simulation
fs = 1/Ts		%Fr�quence d'�chantillonnage

npb = round(Tb/Ts)		%Nombre de points par bit

t = [Ts:Ts:nsb*Tb];     %Vecteur repr�sentant l'axe du temps

%%%%%%%%%%%%%%%%%%%%%%%% EMISSION %%%%%%%%%%%%%%%%%%%%%%%%%
%seqbitsinit = [-1 1 1 -1 1 -1 -1 -1 1 1];
%seqbitsinit = randi([-1 -1], 1,nsb)
seqbitsinit = ones(1,nsb) - floor(rand(1,nsb)*2)*2;
seqbits = sqrt(Eb) .* seqbitsinit;          

%-------------Demultiplexer
path_I = seqbits(1:2:length(seqbits)-1);	%Sous-suite de bits impairs
path_Q = seqbits(2:2:length(seqbits));		%Sous-suite de bits pairs

qpsk_I=path_I(ones(2*npb,1),:);
qpsk_I=qpsk_I(:).';
qpsk_Q=path_Q(ones(2*npb,1),:);
qpsk_Q=qpsk_Q(:).';

figure
plot(t,qpsk_I)


%tc=0:Ts:(10*Tb)-Ts;
 tc=0:Ts:(nsb*Tb)-Ts;%-------------Melangeurs
 qpsk_I_m = qpsk_I .*(sqrt(2/T)*cos(2*pi*fc*tc));                      % A completer
 qpsk_Q_m = qpsk_Q .*(sqrt(2/T)*sin(2*pi*fc*tc));
  
Tx=qpsk_I_m+qpsk_Q_m;%-------------Additionneurs


%-------------Affichage
figure
% subplot(6,1,1)
plot(tc,qpsk_I_m)
figure
% subplot(6,1,2)
plot(tc,qpsk_Q_m)
figure

% subplot(6,1,3)
plot(tc,Tx)
%subplot(6,1,4)
%plot(tc,qpsk_I)
%subplot(6,1,5)
%splot(tc,path_I)


% A completer

%%%%%%%%%%%%%%%%%%%%%%%% RECEPTION %%%%%%%%%%%%%%%%%%%%%%%%%

% Multiplieurs : la multiplication de vecteurs se fait avec l'op�rateur .*
% pi est une variable de Matlab, sqrt est la fonction racine carr�e
%qpskI = ;  %sortie du multiplieur du chemin I
%qpskQ = ;  %sortie du multiplieur du chemin Q
tc=0:Ts:(nsb*Tb)-Ts;
qpskI = Tx.*(sqrt(2/T)*cos(2*pi*fc*tc));
qpskQ = Tx.*(sqrt(2/T)*sin(2*pi*fc*tc));
num   = [Ts 0];     %Num�rateur de l'int�grateur
denom = [1 -1];     %D�nominateur de l'int�grateur
intI = [];
intQ = [];
for m = 1:1:Nsymb
    % Integrateurs
    YI = filter(num,denom,qpskI((m-1)*length(t)/Nsymb+1:m*length(t)/Nsymb));
    YQ = filter(num,denom,qpskQ((m-1)*length(t)/Nsymb+1:m*length(t)/Nsymb));
    intI = cat(2,intI,YI);  %Sortie de l'int�grateur du chemin I
    intQ = cat(2,intQ,YQ);  %Sortie de l'int�grateur du chemin Q
    
    % D�tecteur (comparaison avec 0) et Multiplexeur
    %ATTENTION les 0 ont la valeur -1 (les 1 sont des 1)
    %resbits(...)=
        
end
R_I=[];
R_Q=[];
periode = length(YI);
 td=periode:periode:length(intI);
R_I=intI(td);
R_Q=intQ(td);
SNR = 20;
NR_I = awgn(R_I,SNR,0);
NR_Q = awgn(R_Q,SNR,0);
%scatter(R_I,R_Q,'d'),grid;
scatter(NR_I,NR_Q),grid;
evm = comm.EVM;
rmsEVM1 = evm(R_I,NR_I);
release(evm);
% evm.ReferenceSignalSource = 'Estimated from reference constellation';
% evm.ReferenceConstellation = R_I;
% rmsEVM2 = evm(NR_I);
% [rmsEVM1 rmsEVM2]

% figure
% plot(R_I)
% path_I_rx = R_I(ones(npb,1),:);
% path_Q_rx = R_Q(ones(npb,1),:);
% path_I_rx=reshape(path_I_rx,1,prod(size(path_I_rx)));
% tb=0:5*Tb/(length(path_I_rx)-1):5*Tb;
% path_Q_rx=reshape(path_Q_rx,1,prod(size(path_Q_rx)));
% path_I_r=path_I_rx(1:npb:length(path_I_rx));
% path_Q_r=path_Q_rx(1:npb:length(path_Q_rx));
% nrz_I_rx = sign(path_I_r); 
% nrz_Q_rx = sign(path_Q_r);
% 
% 
% cn = cat(1,nrz_I_rx,nrz_Q_rx);
% resbits = reshape(cn,1,prod(size(cn)));
% figure
% plot(resbits)
% resbits = resbits(ones(npb,1),:);
% resbits=reshape(resbits,1,prod(size(resbits)));
% tbits=0:10*Tb/(length(resbits)-1):10*Tb;
% figure
% plot(tbits,resbits)


%if (resbits == seqbits)
%    fprintf(1,'ok')
%else
%    fprintf(1,'ko')
%end

% AFFICHAGE : vous pouvez utiliser plot directement pour l'affichage d'une
% courbe dans une fenetre. Sinon vous pouvez l'utiliser avec subplot pour
% diviser la fenetre en plusieurs partie et afficher plusieurs courbes sans
% les superposer. HELP plot ou subplot dans Matlab pour plus de details.

%figure
%subplot(3,1,1)
%plot(...)
% A completer

%%%%
%FTT%%
%%%%%
%Fs = 1000;            % Sampling frequency
%T = 1/Fs;             % Sampling period
L = 10000;             % Length of signal
t = (0:L-1)*Ts;        % Time vector


    t = 0:1/fs:1-1/fs;
    x =Tx;
    N = length(x);
    xdft = fft(x);
    xdft = xdft(1:N/2+1);
    psdx = (1/(fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:fs/length(x):fs/2;


    
