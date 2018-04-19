%All simulation is done normalized to the sampling frequnecy (fs=1)
%This m-file able to synthesize NTf, Poles&Zeros Plot, Plotting NTF, plotting output spectrum, calculating SNR and calculating Coefficients for target NTF and modulator architecture.  

close all
order=4;            %Modulator Order
OSR=24;             %Over Samplitng Ration (modulator bandwidth BW is determined now equals to fs/2OSR & fs=1)
opt=0;              %Zeros optimization (0=no opt, 1=yes)
H_inf=1.5;          %Poles optimization for maximum flat out NTF and stability (1.5 Default)
f0=0;               %Resonance frequency normalized to fs (0=low pass)
nLev=2;             %Number of quantizer levels
Nfft=2^14;          %Number of FFT points
fin_bin=13;         %Input frequency bin (should be prime number and gives input frequency less than the bandwidth which equals fs/2OSR)
amp=0.7;            %Amplitude of input signal U
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
figure(3)
hold on
plot(f, dbv(spec(1:Nfft/2+1)),'b')                          %Plot of output spectrum
hold on;

% plot(f, 20*log10(spec(1:Nfft/2+1)),'r')                          %Plot of output spectrum
% hold on;

% plot (f,dbp(sqq*NBW),'g','linewidth',1)                     %plotting NTF shape
% figure
% semilogx(f, dbv(spec(1:Nfft/2+1)))
grid
[a,g,b,c]=realizeNTF(ntf,form)                        %Calculation of Coefficient according to target NTf and momdulator architecture

