% -----------------------------------------
% SIGMA DELTA
% -----------------------------------------
%All simulation is done normalized to the sampling frequnecy (fs=1)
%This m-file able to synthesize NTf, Poles&Zeros Plot, Plotting NTF, plotting output spectrum, calculating SNR and calculating Coefficients for target NTF and modulator architecture.  

close all
Nbit=16;
order=2;            %Modulator Order
OSR=64;             %Over Samplitng Ration (modulator bandwidth BW is determined now equals to fs/2OSR & fs=1)
opt=0;              %Zeros optimization (0=no opt, 1=yes)
H_inf=0;          %Poles optimization for maximum flat out NTF and stability (1.5 Default)
f0=0;               %Resonance frequency normalized to fs (0=low pass)
nLev=2;             %Number of quantizer levels
Nfft=2^14;          %Number of FFT points
fin_bin=13;         %Input frequency bin (should be prime number and gives input frequency less than the bandwidth which equals fs/2OSR)
amp=0.55;            %Amplitude of input signal U
fin=fin_bin/Nfft;   %Input frequency (fin=fin_bin * fs/NFFT)
form='CIFB';        %Modulator Architecture used)
ntf=synthesizeNTF(order,OSR,opt,H_inf,f0);      %Relazing the noise tranfere function characteristics
% plotPZ(ntf)                                     %Poles and Zeros plot 
% figurI
f=linspace(0,0.5,Nfft/2+1);
z=exp(2i*pi*f);
% plot(f, dbv(evalTF(ntf,z)))                     %Plotting the NTF
% figure
t=[0:Nfft-1];                                   %Run time of input signal (should equals to (NFFT-1)*(1/fs)) 
u=amp*(nLev-1)*sin(2*pi*fin*t);                 %Input signal
v=simulateDSM(u,ntf,nLev);                      %Simulation of input signal with target Noise Transfere Function)
%window=blackman(Nfft);
%v_w=v.*window';
spec=fft(v)/(Nfft*(nLev-1)/2);                            %Performing FFt on the output
% spec=fft(v).*conj(fft(v))/Nfft; 
% spec=spec/max(spec);

snr=calculateSNR(spec(1:ceil(Nfft/(2*OSR))+1),fin_bin)     %Calculating SNR              
NBW=1/Nfft;                                                 %Frequency Bin in FFT (FFT resolution equals to fs/Nfft)
sqq=4*(evalTF(ntf,exp(2i*pi*f))/(nLev-1)).^2/3; 
% figure(3)
% hold on
% plot(f, dbv(spec(1:Nfft/2+1)),'b');                          %Plot of output spectrum
% hold on;

% plot(f, 20*log10(spec(1:Nfft/2+1)),'r')                          %Plot of output spectrum
% hold on;

% plot (f,dbp(sqq*NBW),'g','linewidth',1)                     %plotting NTF shape
% figure
% semilogx(f, dbv(spec(1:Nfft/2+1)))
grid
[a,g,b,c]=realizeNTF(ntf,form)                        %Calculation of Coefficient according to target NTf and momdulator architecture

amax = amp;
amin =0;
% Quantize Encoder
delta = (amax-amin)/(2^Nbit);
bound =[amin-(delta/2):delta:amax+(delta/2)];
codebook = [amin : delta : amax ];



%-------------------------------------------------------------------
%			Trace.m
%------------------------------------------------------------------
%echo off
%clear all
%clear global

%---------------------------
% Feedback coefficients
%---------------------------
a1=1;
a2=3/4;
%---------------------------
% Integrator Gain and pole
%---------------------------
fs = 1;
%Aint =0.5;
C = 1;
R = 2/C;
Ao=1000;
%GBW = 20*Aint*fs;
GBW=2*pi*5*fs;
p1 = GBW/Ao;

global yfen yfft Psig Pbruit L

N=16*1024; L=1.2*N;

fin=0;
ai=0;
% nom=input('Nom de fichier sans extension: ','s');
nom= 'SD_CIFB_2nd';
fe=1; fin=(13/N)*fe; Te=1/fe; 
%fin
BW=fe/(2*OSR)
bi=round(N*(BW/fe)+1)
nc=round(N*(fin/fe)+1) 
dN=5;
ai=0;
%amp=[0.1:0.05:1];
amp=[0.55];

if bi-dN<nc
bi-dN
nc
error
end 

%fircoeff=fircoeff_raw;
%%%%%%%OOB for the WIFI IEEE802.11a operating at 5GHz%%%%%%%%
% ch_bw is the channel Bandwidth it could be (10MHz, 20MHz, 50MHz, 90MHz)
% fs if the sampling frequency
% OSR = 32
% ch_bw = 50e+6;
% fs = 2*OSR*ch_bw; 
% abo_fag=[0 (ch_bw-1)/fs   (ch_bw+1)/fs   (2*ch_bw)/fs  (3*ch_bw)/fs   (4*ch_bw)/fs (4*ch_bw)/fs     0.5];
% abo_mag=[0      0           -20          -28         -40          -40               -128          -128  ];
% fircoeff = [1.01335456559147e-05 2.77320918336882e-05 6.36808243157193e-05 0.000127658824398318 0.000233296924959917 0.000397806209910628 0.000641970679163992 0.000989822267591615 0.00146795351978705 0.00210443075451286 0.00292731394861029 0.00396283609054794 0.00523331721401924 0.00675493852533352 0.00853553821547035 0.0105726010234585 0.0128516288437996 0.0153450744563062 0.0180119754626877 0.0207984334885326 0.0236389413860228 0.0264586063919067 0.0291761196236615 0.0317074007575902 0.0339696309439282 0.0358854695140235 0.0373871670562284 0.0384202598278323 0.0389466177016061 0.0389466177016061 0.0384202598278323 0.0373871670562284 0.0358854695140235 0.0339696309439282 0.0317074007575902 0.0291761196236615 0.0264586063919067 0.0236389413860228 0.0207984334885326 0.0180119754626877 0.0153450744563062 0.0128516288437996 0.0105726010234585 0.00853553821547035 0.00675493852533352 0.00523331721401924 0.00396283609054794 0.00292731394861029 0.00210443075451286 0.00146795351978705 0.000989822267591615 0.000641970679163992 0.000397806209910628 0.000233296924959917 0.000127658824398318 6.36808243157193e-05 2.77320918336882e-05 1.01335456559147e-05];
% Apass = 5 ; Astop = 105; Fpass = ch_bw/fs; Fstop = 1.5*abo_fag(3);
% Filter Order = 57

%%%%%%%%OOB for UMTS 5MHz channel BW operating at 1.98GHz%%%%%%%%%%
% ch_bw = 5E6;
% fs = 2*OSR*ch_bw;
abo_fag=[0  30E6/fs   30E6/fs        160E6/fs      160E6/fs      220E6/fs  0.5];
abo_mag=[0      0           -85          -85         -128          -128    -128    ];
% fircoeff=[1.60695392426680e-05 2.95841123026030e-05 5.48262790927617e-05 9.19832588239474e-05 0.000143686658983094 0.000212359882657570 0.000299965229015387 0.000407761225657029 0.000536071876123679 0.000684099631525166 0.000849796105711060 0.00102983168952252 0.00121960996869630 0.00141347333081234 0.00160490224533182 0.00178687182515771 0.00195224251942712 0.00209417001193910 0.00220655685058343 0.00228443101865109 0.00232427271386166 0.00232427271386166 0.00228443101865109 0.00220655685058343 0.00209417001193910 0.00195224251942712 0.00178687182515771 0.00160490224533182 0.00141347333081234 0.00121960996869630 0.00102983168952252 0.000849796105711060 0.000684099631525166 0.000536071876123679 0.000407761225657029 0.000299965229015387 0.000212359882657570 0.000143686658983094 9.19832588239474e-05 5.48262790927617e-05 2.95841123026030e-05 1.60695392426680e-05];
% Apass = 35 ; Astop = 98; Fpass = ch_bw/fs; Fstop = 1.45*abo_fag(3);
% Filter Order = 41

%%%%%%%%OOB for ISM Band and WiFi IEEE802.11g operating at 2.4GHz%%%%%%%%%%
ch_bw = 5E6;
fs = 2*OSR*ch_bw;
abo_fag=[0  50E6/fs   50E6/fs        150E6/fs      150E6/fs      250E6/fs  0.5];
abo_mag=[0      0           -85          -85         -128          -128    -128    ];
fircoeff=[1.60695392426680e-05 2.95841123026030e-05 5.48262790927617e-05 9.19832588239474e-05 0.000143686658983094 0.000212359882657570 0.000299965229015387 0.000407761225657029 0.000536071876123679 0.000684099631525166 0.000849796105711060 0.00102983168952252 0.00121960996869630 0.00141347333081234 0.00160490224533182 0.00178687182515771 0.00195224251942712 0.00209417001193910 0.00220655685058343 0.00228443101865109 0.00232427271386166 0.00232427271386166 0.00228443101865109 0.00220655685058343 0.00209417001193910 0.00195224251942712 0.00178687182515771 0.00160490224533182 0.00141347333081234 0.00121960996869630 0.00102983168952252 0.000849796105711060 0.000684099631525166 0.000536071876123679 0.000407761225657029 0.000299965229015387 0.000212359882657570 0.000143686658983094 9.19832588239474e-05 5.48262790927617e-05 2.95841123026030e-05 1.60695392426680e-05];
% Apass = 35 ; Astop = 98; Fpass = ch_bw/fs; Fstop = 1.45*abo_fag(3);
% Filter Order = 41

% mask_f=[0 5/640 30/640 30/640 130/640 130/640];
% mask_mag=[0   0  0 -85 -85 -124 ];


% u = 0.158;
% frac=(u-floor(u))*100
% arg2 = round(log2(frac - 2^(floor(log2(frac)))))
% arg3=ufi(u, length(de2bi(floor(u)))+arg2, arg2)


fircoeff_fp = 0;
for compt=1:length(fircoeff)
    u = fircoeff(compt);
    frac=(u-floor(u))*1E9;
    arg2 = round(log2(frac - 2^(floor(log2(frac)))));
    temp = ufi(u, length(de2bi(floor(u)))+arg2, arg2);
    fircoeff_fp(compt)=temp.data;
    sizebit(compt)=length(temp.bin);
end    
%---------------------------------


unitcell = sum(fircoeff);
small = min(fircoeff);
fircoeff_round = fircoeff ./ small; %normalizing by the smallest value
fircoeff_round = round(fircoeff_round);
maxmin_diff = peak2peak(fircoeff_round)



for compt=1:length(amp)
ai=amp(compt);
%---------------------------------------
% simulation simulink avec 
%Stop time = L
%Amplitude de la sinusoide = ai
%Frequence de la sinusoide = fin
%Taille du vecteur de sortie = N
%----------------------------------------
[TT,X,Y]=sim(nom,L);
Sp=SPEC(y,N); % fonction qui effectue le fenetrage et calcule la densite spectrale de puissance
%Sp=nowin(y,N);
F=Freq(fe,N);
snr(compt)=SNR(Sp,nc,dN,bi); % calcul du Rapport Signal Sur Bruit
end


%---------------------------------
for compt=1:length(amp)
ai=amp(compt);
%---------------------------------------
% simulation simulink avec 
%Stop time = L
%Amplitude de la sinusoide = ai
%Frequence de la sinusoide = fin
%Taille du vecteur de sortie = N
%----------------------------------------
[TT,X,Y]=sim(nom,L);
Sp1=SPEC(y1,N); % fonction qui effectue le fenetrage et calcule la densite spectrale de puissance
%Sp=nowin(y,N);
F=Freq(fe,N);
snr1(compt)=SNR(Sp,nc,dN,bi); % calcul du Rapport Signal Sur Bruit
end


for compt=1:length(amp)
ai=amp(compt);
%---------------------------------------
% simulation simulink avec 
%Stop time = L
%Amplitude de la sinusoide = ai
%Frequence de la sinusoide = fin
%Taille du vecteur de sortie = N
%----------------------------------------
[TT,X,Y]=sim(nom,L);
Sp2=SPEC(y2,N); % fonction qui effectue le fenetrage et calcule la densite spectrale de puissance
%Sp=nowin(y,N);
F=Freq(fe,N);
snr2(compt)=SNR(Sp,nc,dN,bi); % calcul du Rapport Signal Sur Bruit
end


%affichage du Spectre du signal de sortie
 figure(1)
 low = 15;
 Spdb=10*log10(Sp)- low;
 Spdb(14) = Spdb(14) + low;
 
 Spdb1=10*log10(Sp1)- low;
 Spdb1(14) = Spdb1(14) + low;
 
 Spdb2=10*log10(Sp2)- low;
 Spdb2(14) = Spdb2(14) + low;
%semilogx(F,Spdb);
plot(F,Spdb,'k');%y : raw coeff 
hold on
%plot(F,Spdb1,'c');%y1 : rounded coeff
hold on 
%plot(F,Spdb2,'g');%y2 : fixed point coeff
hold on 
plot(abo_fag,abo_mag,'r'); %with mask test
%legend('raw coeff','rounded coeff','fixed point coeff', ['MaxMin = ' num2str(maxmin_diff)])
grid
