function [lbot,blmult,secang,secsun,xsun,ermono,rplanck] = ...
   preprad(vchunk,indv,ip,prof,fmono,datadir);

% function [lbot,blmult,secang,secsun,xsun,ermono,rplanck] = ...
%   preprad(vchunk,indv,ip,prof,fmono);
%
% Prepare for a radiance calculation by computing terms independent
% of the atmospheric transmittance.
%
% Input:
%    vchunk - [1 x 1] kCARTA chunk start frequency
%    indv   - [nmono x 1] desired freq indices in chunk
%    ip     - [1 x 1] index of desired profile in "prof"
%    prof   - RTP profile structure
%    fmono  - [nmono x 1] frequencies {cm^-1}
%
% Output:
%   lbot   - [1 x 1] bottom layer index
%   blmult - [1 x 1] bottom layer multiplier
%   secang - [1 x nlay] local path secant angle
%   secsun - [1 x nlay] round trip sun secant angle
%   xsun   - [nmono x 1] cossun*omega*hsun
%   ermono - [nmono x 2] emis (index=1) and rho (index=2) for fmono
%   rplanck- [nmono x nlay] Planck emission at layer mean temperature 
%

% Created: 17 Jun 2009, Scott Hannon -  based on rtchunk
% Update: 23 Jun 2009, S.Hannon - fix secang bug (was solzen not satzen)
% Update: 29 Jun 2009, S.Hannon - error traps for nemis=1 and no rho
%
% Copyright 2012, Univ. Of Md, Balt. Co. Atmospheric Spectroscopy Laboratory
% kcarta is distributed under the terms of the GNU GPL v3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edit this section as needed

% Degrees to radians conversion
deg2rad = pi/180;

% directory for matlab solar spectral data
solardir = fullfile(datadir,'kcmix2_solarV2'); % taro

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nlay = prof.nlevs(ip) - 1;
npts = length(indv);

% Interpolate emissivity & reflectivity to fmono
nemis = prof.nemis(ip);
if (nemis > 1)
   inde = 1:nemis;
   X = prof.efreq(inde,ip);
   if (isfield(prof,'rho'))
      Y = [prof.emis(inde,ip), prof.rho(inde,ip)];
   else
      Y = [prof.emis(inde,ip), (1-prof.emis(inde,ip))/pi];
   end
   fmono = fmono(:); % make sure fmono is [n x 1]
   ermono = interp1(X, Y, fmono, 'linear', 'extrap');
   jj = find(fmono <= prof.efreq(1,ip));
   ermono(jj,1) = prof.emis(1,ip);
   ermono(jj,2) = prof.rho(1,ip);
   jj = find(fmono >= prof.efreq(nemis,ip)); 
   ermono(jj,1) = prof.emis(nemis,ip);
   ermono(jj,2) = prof.rho(nemis,ip);
else
   ermono = ones(length(fmono),2);
   ermono(:,1) = prof.emis(1,ip);
   if (isfield(prof,'rho'))
      ermono(:,2) = prof.rho(1,ip);
   else
      ermono(:,2) = (1-prof.emis(1,ip))/pi;
   end
end

% Determine bottom layer index and fractional multiplier
% Note: do not do fractional bottom layers with deltaP <= 5 mb
plevs = prof.plevs(1:prof.nlevs(ip),ip);
ii = min(find(plevs >= prof.spres(ip)-5));
lbot = min([ii - 1,nlay]);
blmult = (prof.spres(ip) - plevs(lbot))./(plevs(lbot+1)-plevs(lbot));
indlay = 1:lbot;
plays = (plevs(indlay+1)-plevs(indlay)) ./ log(plevs(indlay+1)./plevs(indlay));

% Determine bottom fractional layer mean temperature
% Calc mean pressure for bottom fractional layer
rjunk1 = ( prof.spres(ip) - plevs(lbot) )/log( prof.spres(ip)/plevs(lbot) );
% Do interpolation for fractional bottom layer mean temperature
% assuming T is in linear in log(P)
rjunk2=( prof.ptemp(lbot,ip) - prof.ptemp(lbot-1,ip) )./ ...
   log( plays(lbot)/plays(lbot-1) ); % slope
tbot = rjunk2*log( rjunk1/plays(lbot-1) ) + prof.ptemp(lbot-1,ip);

% Local path zenith angle for each layer
% First calc layer mean altitudes
zlay = 0.5*(prof.palts(indlay,ip) + prof.palts(indlay+1,ip));
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

% Reflected solar
if (prof.solzen(ip) >= 0 & prof.solzen(ip) < 90)

  cossun = cos(prof.solzen(ip)*deg2rad);
  dstsun = 1.496E+11;              % distance from earth to sun
  radsun = 6.956E+8;		   % radius of the sun
  omega = pi * (radsun/dstsun)^2;  % solid angle of sun from earth

  % load solar spectra for our chunk, defines the variables
  % sfrq and srad, the solar radiance data for this chunk
  eval(sprintf('load %s/srad%d', solardir, vchunk));
  srad = srad(indv)*1000;      % W to mW

  % Downward solar path angle
  [sunang] = sunang_conv(prof.solzen(ip),zlay);

  % Total round-trip solar secant angle
  secsun = secang + sec(sunang*deg2rad)'; %'
  xsun = cossun.*omega.*srad(:);

else
  secsun = zeros(1,nlay);
  xsun = zeros(npts,1);
end

% Declare planck array
rplanck = zeros(npts,lbot);
% loop downward over the layers
for ii = 1:(lbot-1)
  rplanck(:,ii) = bt2rad(fmono,prof.ptemp(ii,ip));
end
rplanck(:,lbot) = bt2rad(fmono,tbot);

%%% end of function %%%
