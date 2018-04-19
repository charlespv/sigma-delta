%--------------------
%	SNR.m
%--------------------

function SNR=SNR(Sp,nc,dN,b)

global Psig Pbruit
Psig=Sp(nc);
Pbruit=sum(Sp(1:nc-dN))+sum(Sp(nc+dN:b));


SNR=10*log10(Psig/Pbruit);

