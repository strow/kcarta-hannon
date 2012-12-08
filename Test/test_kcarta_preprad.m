% Program test_kcarta_prerad.m
%
% Test Matlab version of kcarta
%
% $$$ This is a more complicated code than test_kcarta_norm since it uses
% $$$ preprad.m to pre-compute layer dependent variables that don't need to 
% $$$ be recomputed if you change the absorption coefficients, for doing 
% $$$ Jacobians, for example.  If you just want to do one radiances calculation
% $$$ follow the test_kcarta_norm.m example.    
%
% $$$ The separation of gases in this program is not necessary, but is
% $$$ done to show how to get started on doing quick finite difference
% $$$ Jacobians
%
% $$$ This code is set up to compute a complete AIRS and IASI spectrum. If
% $$$ start doing smaller regions, the convolution code will break.

% Sample input file
hfile = '../Data/desert_0725_2004.op.rtp';
% Profile index in hfile
ip = 1;

% Locations of compressed data file, solar data, etc.
datadir = '../Data';

% Use all Mfiles under Matlab/Kcarta (if using the Kcarta package)
addpath(genpath('../'));

% Variables added for each chunk
allfreq        = [];
allrad         = [];
allrad_jac_co2 = [];

tic

% Read the input RTP file
[head, hattr, prof, pattr] = rtpread(hfile);

% Group gas id's.  Mainly to separate ones that may change from constant ones
gid_other = [1]; % gases which vary by profile
gid_back  = [2,3,4,5,6,7,8,9,10,11,12,18,22,23,26,27,51,52,53,54,56,57,60,61,62,63]; % nearly constant gases

% Users responsibility to keep vv on kcarta chunk boundaries.  
% They are every 25 cm-1 starting at 605 cm-1
for vv = 605:25:2805

   % Chunk frequnecy
   vchunk = vv

   kcprof_back       = op_rtp_to_lbl(ip,gid_back,head,prof);
   [absc_back, freq] = kcmixfast(kcprof_back, vchunk,datadir);

   % Now do gid_other
   kcprof_other      = op_rtp_to_lbl(ip,gid_other,head,prof);

   [absc_other, freq] = kcmixfast(kcprof_other, vchunk,datadir);

   absc_bo = absc_back + absc_other;

   % Prepare for radiance calculations
   [lbot,blmult,secang,secsun,xsun,ermono,rplanck] = preprad(vchunk,1:length(freq),ip,prof,freq,datadir);

   % Adjust bottom layer
   indlay = 1:lbot;
   absc_bo(:,lbot)  = absc_bo(:,lbot)*blmult;

   % Compute radiance
   absc  = absc_bo; % + absc_co2;
   stemp = prof.stemp(ip);

   rad   = rtchunkfast(lbot,stemp,secang,secsun,xsun,ermono,rplanck,freq,absc);

%  Examples of perturbations without redoing absorption coeff calcs
%
% Do a surface T perturbation and radiance
% $$$ pertstemp=1;
% $$$ stemp = prof.stemp(ip) + pertstemp;
% $$$ rad_jac_stemp = rtchunkfast(lbot,stemp,secang,secsun,xsun,ermono,rplanck,freq,absc);
%
% Do a CO2 profile shift
% $$$ pertco2=1.01;
% $$$ absc = absc_bo + absc_co2*1.01;
% $$$ stemp = prof.stemp(ip);
% $$$ rad_jac_co2 = rtchunkfast(lbot,stemp,secang,secsun,xsun,ermono,rplanck,freq,absc);
% $$$ 
% $$$ allrad_jac_co2 = [allrad_jac_co2;  rad_jac_co2];

   allrad = [allrad  rad'];
   allfreq = [allfreq freq];

end

toc

% To do AIRS SRF convolutions:
tic
[rairs,fairs] = sconv2(allrad',allfreq',1:2378);
toc

% To do IASI SRF convolutions:
tic
[riasi, fiasi] = xfconvkc(allrad','iasi12992');
% Could also use sconv2 with /asl/matlab/srftest/srftables_iasi_645_2760_v2.hdf'
toc

% $$$ % To do CrIS
% $$$ [rch, wch] = xfconvkc(allrad','crisB1');  hamming
% Right now spectral SRF hdf lookup table does not exist for CrIS, but
% is easily created if needed
