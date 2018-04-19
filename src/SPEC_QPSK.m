%----------------------------
%	SPEC.m
%----------------------------
function SPEC_QPSK=SPEC_QPSK(y,N)

global yfen yfft SPEC_QPSK

fenetre=blackman(N);
%yfen=y.*fenetre;
%fenetre=hann(N);
yfen=y.*fenetre;
%yfen = y;
yfft=fft(yfen);
Py=yfft.*conj(yfft)/(length(yfft));
Py(N/2+1:N)=[ ];		
Py(1:N/2)=2*Py(1:N/2);		
Py=Py/(max(Py));		

SPEC_QPSK=Py;
