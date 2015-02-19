% Program test_kcarta_norm.m
%
% Test Matlab version of kcarta
%
% $$$ The separation of gases in this program is not necessary, but is
% $$$ done to show how to get started on doing quick finite difference
% $$$ Jacobians
%
% $$$ This code is set up to compute a complete AIRS and IASI spectrum. If
% $$$ start doing smaller regions, the convolution code will break.

% Sample input file
% hfile = '../Data/desert_0725_2004.op.rtp';

%hfile = '/Users/strow/Desktop/test.rtp';
hfile = '/asl/s1/strow/rtprod_cris/2013/08/28/test.rtp';
% Profile index in hfile
ip = 60;

% Locations of compressed data file, solar data, etc.
datadir = '../Data';

% Use all Mfiles under Matlab/Kcarta (if using the Kcarta package)
addpath(genpath('../'));

tic

% Read the input RTP file
[head, hattr, prof, pattr] = rtpread(hfile);
btobs = real(rad2bt(head.vchan(:),prof.robs1));
btcal = rad2bt(head.vchan(:),prof.rcalc);
bias = btobs-btcal;
i = find(prof.plat > -30 & abs(bias(406,:)) < 1);

% Do all gases at once, no preprad
%gid = [1,2,3,4,5,6,7,8,9,10,11,12,18,22,23,26,27,51,52,53,54,56,57,60,61,62,63]; % nearly constant gases
gid = [1,2,3,4,5,6,9];

%for k=i
for k=i(1)
   ip = k
% Variables added for each chunk
allfreq        = [];
allrad         = [];
allrad_jac_co2 = [];


% Users responsibility to keep vv on kcarta chunk boundaries.  
% They are every 25 cm-1 starting at 605 cm-1
for vv = 605:25:2805
   % Chunk frequnecy
   vchunk = vv;

   kcprof       = op_rtp_to_lbl(ip,gid,head,prof);
   [absc, freq] = kcmixfast(kcprof, vchunk,datadir);

   rad = rtchunk(ip, prof, absc, freq, datadir);

   allrad = [allrad  rad'];
   allfreq = [allfreq freq];
end

allrad = allrad(:);

save_str = ['save kout' int2str(k) ' allrad'];
eval(save_str);

clear allrad

end


% 
% 
% toc
% 
% % To do AIRS SRF convolutions:
% tic
% % Be careful with sconv2, absolute path to SRF tables encoded
% [rairs,fairs] = sconv2(allrad',allfreq',1:2378);
% toc
% 
% % To do IASI SRF convolutions:
% tic
% [riasi, fiasi] = xfconvkc(allrad','iasi12992');
% % Could also use sconv2 with /asl/matlab/srftest/srftables_iasi_645_2760_v2.hdf'
% toc
% 
% % $$$ % To do CrIS
% % $$$ [rch, wch] = xfconvkc(allrad','crisB1');  hamming
% % Right now spectral SRF hdf lookup table does not exist for CrIS, but
% % is easily created if needed
