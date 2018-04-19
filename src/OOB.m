%%%%%%%OOB for the WIFI IEEE802.11a operating at 3.5GHz%%%%%%%%
% ch_bw is the channel Bandwidth it could be (10MHz, 20MHz, 50MHz, 90MHz)
%fs if the sampling frequency
abo_fag=[0 ch_bw-1/fs   ch_bw+1/fs   2*ch_bw/fs  3*ch_bw/fs_bw   4*ch_bw/fs_bw];
abo_mag=[0      0           -20          -28         -40          -40    ];



%%%%%%%%OOB for UMTS 5MHz channel BW operating at 1.98GHz%%%%%%%%%%
% abo_fag=[0  30E6/fs   30E6/fs        160E6/fs      160E6/fs      220E6/fs  0.5];
% abo_mag=[0      0           -85          -85         -128          -128    -128    ];


%%%%%%%%OOB for ISM Band and WiFi IEEE802.11g operating at 2.4GHz%%%%%%%%%%
% abo_fag=[0  50E6/fs   50E6/fs        150E6/fs      150E6/fs      250E6/fs  0.5];
% abo_mag=[0      0           -85          -85         -128          -128    -128    ];