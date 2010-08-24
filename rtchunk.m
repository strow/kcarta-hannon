function [rad] = rtchunk(ip, prof, absc, freq, datadir);

% function [rad] = rtchunk(ip, prof, absc, freq);
%
% Calculate radiances for a "chunk" of kcarta optical depths.
%
% Input:
%    ip     - [1 x 1] index of desired profile in "prof"
%    prof   - RTP profile structure
%    absc   - [10000 x nlay] top-down optical depths
%    freq   - [10000 x 1] frequencies {cm^-1}
%
% Output:
%   rad     - [10000 x 1] radiance at top of atmosphere
% 

% Created: 14 Jun 2002, Howard Motteler
% Update: 16 Jun 2009, Scott Hannon -  rewrite for speed up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edit this section as needed

% Degrees to radians conversion
deg2rad = pi/180;

% directory for matlab solar spectral data
% solardir = '/home/motteler/kcmix2/solarV2';  
solardir = fullfile(datadir,'kcmix2_solarV2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[npts,nlay] = size(absc);

% make freq a column vector
freq = freq(:);

% ones column to match freq.
onescol = ones(length(freq), 1);

% interpolate emissivity & reflectivity to freq
nemis = prof.nemis(ip);
inde = 1:nemis;
X = prof.efreq(inde,ip);
Y = [prof.emis(inde,ip), prof.rho(inde,ip)];
erfine = interp1(X, Y, freq, 'linear', 'extrap');
jj = find(freq <= prof.efreq(1,ip));
erfine(jj,1) = prof.emis(1,ip);
erfine(jj,2) = prof.rho(1,ip);
jj = find(freq >= prof.efreq(nemis,ip)); 
erfine(jj,1) = prof.emis(nemis,ip);
erfine(jj,2) = prof.rho(nemis,ip);

% Determine bottom layer index and fractional multiplier
% Note: do not do fractional bottom layers with deltaP <= 5 mb
plevs = prof.plevs(1:prof.nlevs(ip),ip);
ii = min(find(plevs >= prof.spres(ip)-5));
lbot = min([ii - 1,nlay]);
blmult = (prof.spres(ip) - plevs(lbot))./(plevs(lbot+1)-plevs(lbot));
indlay = 1:lbot;
indlev = 1:(lbot+1);
zlay = 0.5*(prof.palts(indlay,ip) + prof.palts(indlay+1,ip));
xabsc = absc(:,indlay);
xabsc(:,lbot) = xabsc(:,lbot)*blmult;

% Local path zenith angle for each layer
if (prof.satzen(ip) >= 0)
  [zang] = sunang_conv(prof.satzen(ip),zlay);
else
   if (abs(prof.scanang(ip)) < 50)
      [zang] = vaconv(prof.scanang(ip),prof.zobs(ip),zlay);
   else
      error('No satzen or scanang')
   end
end
secang = sec(zang*deg2rad)'; %'

% --------------------------------
% Reflected thermal at arccos(3/5)
% --------------------------------
% Transmittance at diffusivity angle
sectherm = 5/3; % sec( acos(3/5) )
tran = exp(-xabsc .* sectherm);

% Radiance of empty space
rthm = bt2rad(freq, 2.7);

% Declare planck array
rplanck = zeros(npts,lbot);

% loop downward over the layers
for ii = 1:lbot
  rplanck(:,ii) = bt2rad(freq,prof.ptemp(ii,ip));
  rthm = rthm.*tran(:,ii) + rplanck(:,ii).*(1 - tran(:,ii));
end
% Note: above rthm is radiance at one representative angle; to convert
% this to the an approximate integral over all angles multiply by pi.

% Upwards reflected component at surface
rthm = rthm .* (1 - erfine(:,1));  % pi*rthm*(1-emis)/pi

% --------------------------------
% Atmospheric and surface emission
% --------------------------------
% Upward radiance at surface
rad = rthm + bt2rad(freq,prof.stemp(ip)).*erfine(:,1);

% Layer transmittances along view angle
tran = exp(-xabsc .* (ones(npts,1)*secang));

% Loop upward over the layers
for ii = lbot:-1:1
  rad = rad.*tran(:,ii) + rplanck(:,ii).*(1 - tran(:,ii));
end

% ---------------
% Reflected solar
% ---------------
if (prof.solzen(ip) >= 0 & prof.solzen(ip) < 90)

  cossun = cos(prof.solzen(ip)*deg2rad);
  dstsun = 1.496E+11;              % distance from earth to sun
  radsun = 6.956E+8;		   % radius of the sun
  omega = pi * (radsun/dstsun)^2;  % solid angle of sun from earth

  % load solar spectra for our chunk, defines the variables
  % sfrq and srad, the solar radiance data for this chunk
  eval(sprintf('load %s/srad%d', solardir, freq(1)));
  srad = srad(:)*1000;      % W to mW

  % Downward solar path angle
  [sunang] = sunang_conv(prof.solzen(ip),zlay);

  % Total round-trip solar secant angle
  secsun = secang + sec(sunang*deg2rad)'; %'

  % Surface-to-space optical depth
  solabs = sum(xabsc .* (ones(npts,1)*secsun), 2);

  % TOA reflected solar radiance
  rad = rad + cossun .* omega .* srad .* exp(-solabs) .* erfine(:,2);
end


%%% end of function %%%
