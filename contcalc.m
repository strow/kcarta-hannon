function [absc] = contcalc(kcprof, freq, datadir)

% function [absc] = contcalc(kcprof, freq, copt);
%
% Calculate water absorption continuum over a requested 
% frequency interval by interpolating from tabulated data
%
% Input:
%   kcprof   - kcmix format profile data with fields:
%      glist - [1 x ngas] gas list {HITRAN}
%      mtemp - [nlay x 1] layer mean temperature {Kelvin}
%      mpres - [nlay x 1] layer mean pressure {atm}
%      gpart - [nlay x ngas] - layer partial pressure of H2O {atm}
%      gamnt - [nlay x ngas] - layer integrated H2O {moleculkes/cm^2}
%   freq   - [n x 1] desired output frequency grid {cm^-1}
%   copt   - optional function parameters:
%      cvers  - continuum version, default '4'
%      cswt   - self continuum adjustment weight, default 1
%      cfwt   - foreign continuum adjustment weight, default 1
%      cdir   - location of tabulated continuum data, default 
%            '/asl/data/kcarta/KCARTADATA/General/CKDieee_le'
%
% Output:
%   absc   - [n x nlay] optical depth
% 
% NOTE:
%   This code is significantly slower than a calculation from compact
%   tabulated coefficients
%

% Created: 14 Jun 2002 Howard Motteler
% Update: 12 Jun 2009 Scott Hannon - speed up code ~4x
% Update: 15 Jun 2009, S.Hannon - allow multiple gases in kcprof
% Update: 26 Jun 2009, S.Hannon - force tL outside ts range to use min/max 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% physcial constants
kAvog = 6.022045e26;
kPlanck2 = 1.4387863;

cdir = fullfile(datadir,'CKDieee_le'); % default dir

% Pass these in function call one day
% LLS: Scott said to sue cvers = 5, email of Aug, 21, 2010
cvers = '5';	% default CKD version
cswt = 1;	% default self-component weight
cfwt = 1;	% default foreign-component weight

indh2o = find(kcprof.glist == 1);
if (length(indh2o) ~= 1)
   error('Unable to find water in kcprof.glist')
end

% Read the continuum files
sfile = [cdir, '/CKDSelf', cvers, '.bin'];
ffile = [cdir, '/CKDFor', cvers, '.bin'];
[ks, junk, ts] = contread(sfile);
[kf, fcoarse, tf] = contread(ffile);
if (length(junk) ~= length(fcoarse))
   sfile
   ffile
   error('Mismatched length of self and foreign continuum files');
end
if (abs(fcoarse-junk) > 0.01)
   sfile
   ffile
   error('Mismatched frequency points in self and foreign continuum files');
end
nts = length(ts);

% Initialize output
nfreq = length(freq);
nlays = length(kcprof.mpres);
absc = zeros(nfreq, nlays);

%-------------------------------
% set up frequency interpolation 
%-------------------------------
df = diff(fcoarse);   % frequency step sizes
n = length(df);       % number of frequency steps
fmid = fcoarse(1:n) + df./2;                 % midpoints of step intervals
imid = interp1(fmid, 1:n, freq, 'nearest');  % indices of midpoints
ineed = unique([imid, imid(length(imid)) + 1]);
fw1 = (fcoarse(imid+1) - freq) ./ df(imid);  % weight for lower coarse point
fw2 = (freq - fcoarse(imid)) ./ df(imid);    % weight for upper coarse point
% Note: kfine = fw1 .* ks(imid) + fw2 * ks(imid+1);

% Declare work arrays
np1= n + 1;
cscoarse = zeros(1,np1);
cfcoarse = zeros(1,np1);
odcoarse = zeros(1,np1);

% loop on layers
for iL = 1 : nlays

  % Current layer mean temperature {Kelvin}
  tL = kcprof.mtemp(iL);

  % ------------------------------------------
  % temperature interpolate the self-continuum
  % ------------------------------------------
  % table index of greatest lower bound temp
  i1 = max([find(ts <= tL); 1]);
  % table index of least upper bound temp
  i2 = min([find(tL <= ts); nts]);
  % get temperature interpolation weights
  if ts(i2) ~= ts(i1)
    tw2 = (tL - ts(i1)) / (ts(i2) - ts(i1));
    tw1 = 1 - tw2;
  else
    tw2 = 1;
    tw1 = 0;
  end
  cscoarse(ineed) = ( (tw1.*ks(i1,ineed)) + (tw2.*ks(i2,ineed)) ).*cswt;

  % ---------------------------------------------
  % temperature interpolate the foreign continuum
  % ---------------------------------------------
  % table index of greatest lower bound temp
  % The foreign continuum has no temperature dependence
  %  i1 = max(find(tf <= tL));
  %  % table index of least upper bound temp
  %  i2 = min(find(tL <= tf));
  %  % get temperature interpolation weights
  %  if tf(i2) ~= tf(i1)
  %    tw2 = (tL - tf(i1)) / (tf(i2) - tf(i1));
  %    tw1 = 1 - tw2;
  %  else
  %    tw2 = 1;
  %    tw1 = 0;
  %  end
  %
  %  cfcoarse(ineed) = ( (tw1.*kf(i1,ineed)) + (tw2.*kf(i2,ineed)) ).*cfwt;
  %%%
  cfcoarse(ineed) = kf(1,ineed);
  %%%

  % -----------------------------------
  % combine self and foreign components  
  % -----------------------------------
  % scalar values for layer iL
  pL = kcprof.mpres(iL);		   % layer pressure, atms
  ppL = kcprof.gpart(iL,indh2o); 	   % water partial pressure, atms
  qL = kcprof.gamnt(iL,indh2o);            % water gas amount
  a1 = qL * kAvog * 296.0 / tL;
  a2 = kPlanck2 / (2 * tL);
  %
  odcoarse(ineed) = (cscoarse(ineed).*ppL + cfcoarse(ineed).*(pL - ppL)) ...
     .* fcoarse(ineed) .* tanh(a2.*fcoarse(ineed)) .* a1;

  % Interpolate from coarse to output freqs
  absc(:, iL) = (odcoarse(imid).*fw1 + odcoarse(imid+1).*fw2); %

end % loop on layers
