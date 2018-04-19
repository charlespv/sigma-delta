%-------------------------------------------------------------------
%			Trace.m
%------------------------------------------------------------------
echo off
%clear all
clear global
 close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 

%---------------------------
% Feedback coefficients
%---------------------------
% a1=0.0450; b1=a1;
% a2=0.2332; b2=a2;
% a3=0.5589; b3=a3;
% a4=1;

global yfen yfft Psig Pbruit L

N=2^14; L=N*1.2;
OSR_QPSK=64;
Beta=0.2;       %Raised Coine Roll-of Factor
fe=1;           %I am going to work on normalized sampling Frequency
Te=1/fe;        %Sample Time


C_aprox=[1.60695392426680e-05 2.95841123026030e-05 5.48262790927617e-05 9.19832588239474e-05 0.000143686658983094 0.000212359882657570 0.000299965229015387 0.000407761225657029 0.000536071876123679 0.000684099631525166 0.000849796105711060 0.00102983168952252 0.00121960996869630 0.00141347333081234 0.00160490224533182 0.00178687182515771 0.00195224251942712 0.00209417001193910 0.00220655685058343 0.00228443101865109 0.00232427271386166 0.00232427271386166 0.00228443101865109 0.00220655685058343 0.00209417001193910 0.00195224251942712 0.00178687182515771 0.00160490224533182 0.00141347333081234 0.00121960996869630 0.00102983168952252 0.000849796105711060 0.000684099631525166 0.000536071876123679 0.000407761225657029 0.000299965229015387 0.000212359882657570 0.000143686658983094 9.19832588239474e-05 5.48262790927617e-05 2.95841123026030e-05 1.60695392426680e-05];

%C_aprox = [1.01335456559147e-05 2.77320918336882e-05 6.36808243157193e-05 0.000127658824398318 0.000233296924959917 0.000397806209910628 0.000641970679163992 0.000989822267591615 0.00146795351978705 0.00210443075451286 0.00292731394861029 0.00396283609054794 0.00523331721401924 0.00675493852533352 0.00853553821547035 0.0105726010234585 0.0128516288437996 0.0153450744563062 0.0180119754626877 0.0207984334885326 0.0236389413860228 0.0264586063919067 0.0291761196236615 0.0317074007575902 0.0339696309439282 0.0358854695140235 0.0373871670562284 0.0384202598278323 0.0389466177016061 0.0389466177016061 0.0384202598278323 0.0373871670562284 0.0358854695140235 0.0339696309439282 0.0317074007575902 0.0291761196236615 0.0264586063919067 0.0236389413860228 0.0207984334885326 0.0180119754626877 0.0153450744563062 0.0128516288437996 0.0105726010234585 0.00853553821547035 0.00675493852533352 0.00523331721401924 0.00396283609054794 0.00292731394861029 0.00210443075451286 0.00146795351978705 0.000989822267591615 0.000641970679163992 0.000397806209910628 0.000233296924959917 0.000127658824398318 6.36808243157193e-05 2.77320918336882e-05 1.01335456559147e-05];


bin=N/4;
fc=(bin/N)*fe;
fc=fe/4;

%%%%For Standard with multi Channles, DS works on standard BW, up convvert the Standard%%%%%%%
%%%%%-------------------------------------------------------------%%%
% Ch_N=10;        %Number of Channles
% Tsymb=round((Ch_N*OSR_QPSK)*(1+Beta)*Te/2);            %Normally the BW of the rect bit has aprox value of 1/2Tsymobpe (sinc(fT), Raised coise makes the symbole takes the shape of sinc with nullings at the peak of next symbole to prevent ISI the new BW is (1/2Tsymb)(1+Beta), 
%                                                                                 %So its obvious that the minimum DS BW will be the BW mentioned above (channel BW), which means that DS will work on one channel and then up convert it
%                                                                             % All Standards are multi-channels, and as mentioned the channel BW is determined by the TSymbole, here we gona work for QPSK for RF BW of 100MHz, and consists of 10 Channels so we will calculate the channel BW and hence the Tsymbole
% BW_RF=fe/OSR_QPSK;                                                                           
% Samples_N=Tsymb/Te;     %number of samples per symbole of the Raise Coisne
% STD_BW=100E6
% Sampling_Freq=OSR_QPSK*STD_BW
% BW_RFF=Sampling_Freq*BW_RF
% Check=STD_BW==BW_RFF
% Channel_BW=BW_RFF/Ch_N
% BW_LPF=BW_RF;


%%For One Channel BW, DS works on channel BW, Up converting the Channel)%%%
%%%------------------------------------------------------------------%%%%
BW_CH=fe/(2*OSR_QPSK);
Tsymb=round((1/BW_CH)*(1+Beta)*Te/2);
Samples_N=Tsymb/Te;
Channel_BW_RF=5E6
Sampling_Freq=2*OSR_QPSK*Channel_BW_RF
BW_CHH=Sampling_Freq*BW_CH
Check=BW_CHH==Channel_BW_RF
BW_LPF=BW_CH*60;


N_b=12; %number of bits of for DS Quantizer

%amp=0.219;  %Normally the peak will be 1 because we are sending 1 and -1 but now as a since shape, but I but 0.18 coz I am not picking the right spanning factor which must be 2.
amp=0.55;

nom='QPSK_sim'; %input('Nom de fichier sans extension:', 's');

% Amp for the carrier (Modulation)
sin_amp = 1;              

GG=[1];

for compt=1:length(amp)
G=GG(1);    
amppp=amp(compt);
a_D=amppp;                      %amplitude for DS Quantizer
Delta=(2*a_D)/(2^N_b-1);     %Resolution
Bound1=-a_D-Delta/2:Delta:a_D+Delta/2;
code1=-a_D:Delta:a_D;

sim('QPSK_sim');
Sp=SPEC_QPSK(y.signals.values,N); % fonction qui effectue le fenetrage et calcule la densite spectrale de puissance
% Spp=SPEC_QPSK(X14.signals.values,N);
% Spp15=SPEC_QPSK(X15.signals.values,N);
Spy1=SPEC_QPSK(y1.signals.values,N);
% [val nc]=max(Spp);
%Sp=nowin(y,N);
F=Freq_QPSK(fe,N);
end

fs_bw=Sampling_Freq;           %this is the smapling frequency you gona use to obtain the desired BW for the certain OSR you Specified above, it is the smapling frequency of both DS and FIRDAC

%%%%%%%OOB for the WIFI 50MHz%%%%%%%%
% abo_fag=[0  49E6/fs_bw   51E6/fs_bw   100E6/fs_bw   150E6/fs_bw   200E6/fs_bw];
% abo_mag=[0      0           -20          -28         -40          -40    ];

%%%%%%%%OOB for UMTS 5MHz working, Including The In Band Mask %%%%%%%%%%
%abo_fag=[0 5e6/fs_bw     5e6/fs_bw  10e6/fs_bw    10e6/fs_bw    15e6/fs_bw  15e6/fs_bw           30e6/fs_bw   30e6/fs_bw        160e6/fs_bw      160e6/fs_bw      220e6/fs_bw  ];
%abo_mag=[0   0              -33          -33         -43              -43          0                     0           -85                  -85            -128            -128      ];
%abo_mag=10*log10(5e6)+abo_mag;

% %%%%%%%%OOB for UMTS 5MHz%%%%%%%%%%
% abo_fag=[0 5e6/fs_bw        30e6/fs_bw   30e6/fs_bw        160e6/fs_bw      160e6/fs_bw      220e6/fs_bw  ];
% abo_mag=[0    0                  0           -85                  -85            -128            -128      ];
% abo_mag=10*log10(5e6)+abo_mag;


% affichage du Spectre du signal de sortie

 Spdb=10*log10(Sp); %pricision
% Spdbb=10*log10(Spp);
 Spy1db=10*log10(Spy1);
% 
            figure(3)
%            plot(F,Spdbb,'b')
%             hold on
              plot(F,Spdb-25,'r')
%              grid
             hold on
            %plot(F,Spy1db,'m')
             
             plot(abo_fag,abo_mag,'b')
            



%%% EVM
%%%% After Modulation
I_out =  y.signals.values(1:Tsymb:length(y.signals.values));
Q_out = y1.signals.values(1:Tsymb:length(y1.signals.values));
I_out_init =  y.signals.values(1:Tsymb:length(y.signals.values));
Q_out_init = y1.signals.values(1:Tsymb:length(y1.signals.values));
for i=1:length(I_out)
    if (I_out(i) == 0)
        I_out(i) = 0; end 
    if (I_out(i) > 0)
        I_out(i) = 1; end 
    if (I_out(i) < 0)
        I_out(i) = -1; end 
    if (Q_out(i) == 0)
        Q_out(i) = 0; end 
    if (Q_out(i) > 0)
        Q_out(i) = 1; end 
    if (Q_out(i) < 0)
        Q_out(i) = -1; end
    end   
path_I_rx= X14.signals.values(1:Tsymb:length(X14.signals.values));
path_Q_rx= X15.signals.values(1:Tsymb:length(X15.signals.values));
% NR_I = y.signals.values;
% NR_Q = y1.signals.values;
small = min(I_out_init);
I_out_init = I_out_init ./ small;

small = min(Q_out_init);
Q_out_init = Q_out_init ./ small;

small = min(path_I_rx);
path_I_rx = path_I_rx ./ small;

small = min(path_Q_rx);
path_Q_rx = path_Q_rx ./ small;

NR_I = I_out_init;
NR_Q = Q_out_init;
Ierr = path_I_rx - NR_I;
Qerr = path_Q_rx - NR_Q;
EVMref = (1/length(path_I_rx))*sum(power(path_I_rx,2)+power(path_Q_rx,2));
EVMmes = (1/length(NR_I))*sum(power(Ierr,2)+power(Qerr,2));
EVMpercentage_AfterModulation = sqrt(EVMmes/EVMref)*100

plot(I_out_init,Q_out_init, 'h');
hold on
plot(I_out,Q_out, '*');
hold on
plot(path_I_rx,path_Q_rx,'d');
title('QPSK constellation with modulation');
xlabel('real part');
ylabel('imaginary part');
plot_lims = [-2 2];
xlim(plot_lims);
ylim(plot_lims);
figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% before Modulation %%%%%%%%%%%%%%
I_out =  y2.signals.values(1:Tsymb:length(y.signals.values));
Q_out = y3.signals.values(1:Tsymb:length(y.signals.values));
I_out_init =  y2.signals.values(1:Tsymb:length(y.signals.values));
Q_out_init = y3.signals.values(1:Tsymb:length(y.signals.values));
for i=1:length(I_out)
    if (I_out(i) == 0)
        I_out(i) = 0; end 
    if (I_out(i) > 0)
        I_out(i) = 1; end 
    if (I_out(i) < 0)
        I_out(i) = -1; end 
    if (Q_out(i) == 0)
        Q_out(i) = 0; end 
    if (Q_out(i) > 0)
        Q_out(i) = 1; end 
    if (Q_out(i) < 0)
        Q_out(i) = -1; end
    end   
path_I_rx= X14.signals.values(1:Tsymb:length(X14.signals.values));
path_Q_rx= X15.signals.values(1:Tsymb:length(X15.signals.values));
% NR_I = y.signals.values;
% NR_Q = y1.signals.values;

small = min(I_out_init);
I_out_init = I_out_init ./ small;

small = min(Q_out_init);
Q_out_init = Q_out_init ./ small;

small = min(path_I_rx);
path_I_rx = path_I_rx ./ small;

small = min(path_Q_rx);
path_Q_rx = path_Q_rx ./ small;

NR_I = I_out_init;
NR_Q = Q_out_init;
Ierr = path_I_rx - NR_I;
Qerr = path_Q_rx - NR_Q;
EVMref = (1/length(path_I_rx))*sum(power(path_I_rx,2)+power(path_Q_rx,2));
EVMmes = (1/length(NR_I))*sum(power(Ierr,2)+power(Qerr,2));
EVMpercentage_BeforeModulation = sqrt(EVMmes/EVMref)*100

% plot scope
% plot(path_I_rx)
% hold on
% plot(I_out)





plot(I_out_init,Q_out_init, 'h');
hold on
plot(I_out,Q_out, '*');
hold on
plot(path_I_rx,path_Q_rx,'d');
title('QPSK constellation without modulation');
xlabel('real part');
ylabel('imaginary part');
plot_lims = [-2 2];
xlim(plot_lims);
ylim(plot_lims);

%  figure(2) 
%   semilogx(F,Spdb,'r'); %pricision
%   hold on
%  semilogx(F,Spdbb,'b');
%  grid


% 
% Ifir_leakage=(max_UCs-No_unit_cells_aprox)*leakage_curr %both the leakgae current and max_UCs are in tamoora first few lines.
% Ifir_max=max(Bout1.signals.values);
% Ifir_min=min(Bout1.signals.values);
% Ifir_avg=(Ifir_max+Ifir_min)/2 
% Ifir_peak_single_side=Ifir_max-Ifir_avg
% Idiff_fir=2*Ifir_peak_single_side;
% Idc_req=2.2420e-04-Ifir_avg;
% modPFIRDAC_Idc_req=2.2420e-4-(Ifir_avg+Ifir_leakage)


