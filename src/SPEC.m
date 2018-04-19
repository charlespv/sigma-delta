%----------------------------
%	SPEC.m
%----------------------------
function SPEC=SPEC(y,N)

global yfen yfft SPEC

fenetre=blackman(N);
yfen=y.*fenetre;
%fenetre=hanning(N);
%yfen=y.*fenetre;
yfft=fft(yfen);
Py=yfft.*conj(yfft)/(length(yfft));
Py(N/2+1:N)=[ ];		
Py(1:N/2)=2*Py(1:N/2);		
Py=Py/(max(Py));		

SPEC=Py;
