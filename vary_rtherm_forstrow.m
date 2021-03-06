% Copyright 2012, Univ. Of Md, Balt. Co. Atmospheric Spectroscopy Laboratory
% kcarta is distributed under the terms of the GNU GPL v3

% this is the same as /home/sergio/KCMIX2/KCMIXCODE/vary_rtherm.m
% except we use TOA = layer 1, GND = layer 100

%% does more accurate background thermal : 
%%    acos(3/5) everywhere in the upper layers
%%    more selective at lower layers

% L. Strow, June 23, 2014; absc --> xabsc

% radiance calc along reflected thermal path
tspace = 2.7;
rthm = bt2rad(freq, tspace);

%%%% big assumption ... these are 100 AIRS layers
preslevels = load('../Data/airslevels.dat');
%%%% this was done for kCARTA, so 1 = GND, 100 = TOA

%nlevs = nlay + 1;
% LLS
[x1,x2]=size(xabsc);
nlay = x2;
nlevs = nlay + 1;

rCos = 3/5;

iThermalLayer = FindDiffusivityBdry(freq,preslevels,nlevs);
iThermalLayer = 101 - iThermalLayer + 1;  %% need to swap

iLL = nlay;

iN = length(freq);

sumk   = sum(xabsc(:,1:iLL)')';
sumkm1 = sum(xabsc(:,2:iLL)')';

% Declare planck array
rplanck = zeros(npts,lbot);

% loop downward over the layers
for ii = 1:lbot
  rplanck(:,ii) = bt2rad(freq,prof.ptemp(ii,ip));
end

for ii = 1:iThermalLayer
  angle(ii) = rCos;
  iiM1     = ii + 1;
  cos_ii   = rCos;
  cos_iiM1 = rCos;
  Temp     = prof.ptemp(ii);
  raT      = exp(-sumk/cos_ii);
  raTm1    = exp(-sumkm1/cos_iiM1);
  raPlanck = rplanck(:,ii);
  raEmission = (raTm1 - raT).*raPlanck;
  rthm       = raEmission + rthm;
  sumk   = sumkm1;
  sumkm1 = sumkm1 - xabsc(:,iiM1);
  end

for ii=iThermalLayer+1:iLL-1
  %fprintf(1,'b %3i   %8.6f \n',ii,Temp);
  iiM1     = ii + 1;
  cos_ii   = ExpInt3(sumk);
  cos_iiM1 = ExpInt3(sumkm1);
  Temp     = prof.ptemp(ii);
  raT      = exp(-sumk./cos_ii);
  raTm1    = exp(-sumkm1./cos_iiM1);
  raPlanck = rplanck(:,ii);
  raEmission = (raTm1 - raT).*raPlanck;
  rthm       = raEmission + rthm;
  sumk   = sumkm1;
  sumkm1 = sumkm1 - xabsc(:,iiM1);
  end

for ii=iLL:iLL
  %fprintf(1,'b %3i   %8.6f \n',ii,Temp);
  iiM1     = ii + 1;
  cos_ii   = ExpInt3(sumk);
  cos_iiM1 = ExpInt3(0);
  Temp     = prof.ptemp(ii);
  raT      = exp(-sumk./cos_ii);
  raTm1    = exp(-sumkm1/cos_iiM1);
  raPlanck = rplanck(:,ii);
  raEmission = (raTm1 - raT).*raPlanck;
  rthm       = raEmission + rthm;
  end
  
% Upwards reflected component at surface
rthm = rthm .* (1 - erfine(:,1));  % pi*rthm*(1-emis)/pi
