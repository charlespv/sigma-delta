%% SPECTRUM %%%%%

function Freq_QPSK=Freq_QPSK(fe,N)

Freq_QPSK=(fe/2)*(0:N/2-1)/(N/2);
%Freq=(fe/2)*(0:N/2-1);

