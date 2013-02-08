function [absc, fr] = kcmixfast(kcprof, vchunk, datadir);

% function [absc, fr] = kcmixfast(kcprof, vchunk)
% 
% Calculate a 25 cm^-1 "chunk" of mixed optical depths for a
% supplied profile from compressed lookup tables; includes
% water continuum if glist includes water.
%
% Input:
%    kcprof  - [structure] top-down AIRS layers profile with fields:
%       glist  -  [1 x ngas] HITRAN gas IDs
%       mpres  -  [nlay x 1] layers mean pressure {mb}
%       mtemp  -  [nlay x 1] layer mean temperature {K}
%       gamnt  -  [nlay x ngas] layer gas amounts {kilomole/cm^2}
%       gpart  -  [nlay x ngas] layer gas partial pressures {mb}
%    vchunk  - [1 x 1] start frequency of chunk {cm^-1}
%
% Output:
%    absc    - [10000 x nlay] optical depth
%    fr      - [10000 x 1] frequency {cm^-1}
%

% NOTES
%
% The reference profile is the profile that was used in generating
% the coefficient tabulation.  Layers in the reference profile must
% span layers in the input profile.
%
% We assume the reference profile and the coefficient database
% are reliable, so gasses that are not in the reference profile,
% or for which there is no data for a particular 25 cm^-1 interval, 
% are simply skipped.
%
% Modified so that it includes CO2 chi functions for the 2255,2280,
% 2355-2430 cm-1 regions
%
% Reset for package to (Data is passed in now as datadir)
%    'Data/v07.ieee-le/h2o.ieee-le' for H2O
%    'Data/kcarta/v24.ieee-le/co2.ieee-le' for CO2
%    'Data/v07.ieee-le/etc.ieee-le' for all others
%
% BUGS
%
% Various parameters, such as temperature and pressure offsets,
% and the compression exponent, are fixed in the code; it would
% be more general to read these in with the compressed data
%
% HISTORY:
% Created: Howard Motteler, circa 2002
% Update: Sergio Machado, 23 Nov 2005 - add CO2 chi adjustment?
% Update: Scott Hannon, 27 May 2009 - add "od_gind" so that chi will
%    only modify current gas CO2 rather than all mixed gases so far.
%    Make kpath manditory (was optional). Replace load of matlab
%    database with rdgaschunk_le for fortran binary database. Test
%    for any negative optical depth.
% Update: 15 Jun 2009, S.Hannon - "fast" variant created
% Update: 3 Jul 2009, S.Hannon - return fr even if no gases found
% Update: 6 Jul 2009, L. Strow - repackage, mostly data directories
%
% Copyright 2012, Univ. Of Md, Balt. Co. Atmospheric Spectroscopy Laboratory
% kcarta is distributed under the terms of the GNU GPL v3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pre-compression transform aka "k power"
kpow = 1/4;

% temperature tabulation offsets (Kelvin)
toffset = [-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50];

% pressure tabulation offsets for H2O (multiplier to ref partial pressure?)
poffset = [0.1, 1.0, 3.3, 6.7, 10.0];

% Reference profile
refp = 'refpro';

% kCARTA databases
kpathh2o = fullfile(datadir,'v07.ieee-le/h2o.ieee-le');
kpathco2 = fullfile(datadir,'v24.ieee-le/co2.ieee-le');
kpathetc = fullfile(datadir,'v07.ieee-le/etc.ieee-le');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check input
if (nargin ~= 3)
   error('Unexpected number of input arguments')
end
if (kcprof.mpres(2) < kcprof.mpres(1))
   error('kcprof must be a top-down AIRS layer profile')
end
ngas = length(kcprof.glist);  % number of gasses in input profile
nlay = length(kcprof.mpres);  % number of layers in input profile

% load reference profile (defines profile structure "refpro")
eval(sprintf('load %s', refp));

% initialize absorption output array
absc = zeros(1e4, nlay);

% initialize frequency array
fr = (vchunk + (0:9999)*.0025); %'

% Determine temperature interpolation weights for each layer
itwlo = zeros(nlay,1);
itwhi = zeros(nlay,1);
twlo  = zeros(nlay,1);
twhi  = zeros(nlay,1);
for iL = 1:nlay
   iLr = 101 - iL; % bottom-up layer index for kcarta and refpro
   tL = kcprof.mtemp(iL); % temperature of current layer
   tspan1 = toffset + refpro.mtemp(iLr);
   itlo(iL) = max([find(tspan1 <= tL), 1]);
   ithi(iL) = min([find(tL <= tspan1), length(tspan1)]);
   t11 = tspan1(itlo(iL));  % lower temperature bound
   t12 = tspan1(ithi(iL));  % upper temperature bound
   if t11 ~= t12
     % temperature interpolation weight
     twhi(iL) = (tL - t11) / (t12 - t11);
     twlo(iL) = 1 - twhi(iL);
   else 
     twhi(iL) = 1;
     twlo(iL) = 0;
   end
end

% loop on gasses in the supplied profile
for gind = 1:ngas
  gid = kcprof.glist(gind);

  % get file name of compressed data for this gas and chunk
  if (gid >= 3)
    cgfile = sprintf('%s/r%d_g%d.dat', kpathetc, vchunk, gid);
  else
    if (gid == 1)
      cgfile = sprintf('%s/r%d_g1.dat', kpathh2o, vchunk);
    else
      cgfile = sprintf('%s/r%d_g2.dat', kpathco2, vchunk);
    end
  end

  % index of current gas ID in the reference profile
  rgind = find(refpro.glist == gid);

  % check that we have reference and compressed data for this gas
  if ~isempty(rgind) & exist(cgfile) == 2
    
    % load compressed coefficient file, defines vars:
    %
    %   B       10000 x d         absorption basis
    %   fr          1 x 10000     associated frequencies
    %   gid         1 x 1         HITRAN gas ID
    %   kcomp       d x 100 x 11  compressed coefficients
    % or for H2O only
    %   kcomp       d x 100 x 11 x 5 compressed coefficients
    %
    [fr, fstep, toffset, kcomp, B, gid, ktype] = rdgaschunk_le(cgfile);
    [n, d] = size(B);

    % space for compact absorption coeffs for current gas
    kcmp1 = zeros(d, nlay);

    % loop on layer indices
    for iL = 1:nlay
      iLr = 101 - iL; % bottom-up layer index for kcarta and refpro
      if gid ~= 1
        % Gas is not water; do temperature interpolation
        kcmp1(:,iL) = kcomp(:,iLr,itlo(iL))*twlo(iL) + ...
                      kcomp(:,iLr,ithi(iL))*twhi(iL);
      else
        % Gas is water; do temperature and partial pressure interp

        % get partial pressure interval and interpolation weight
	% (partial pressures are tabulated in increasing order)
        qL = kcprof.gpart(iL,gind);       % partial pressure, this layer

	% get partial pressure tabulation bounding interval [q1, q2]
        qspan1 = poffset * refpro.gpart(iLr,rgind);
        iq11 = max([find(qspan1 <= qL), 1]);
        iq12 = min([find(qL <= qspan1), length(poffset)]);
        q11 = qspan1(iq11);
        q12 = qspan1(iq12);
        if q11 ~= q12
	  % partial pressure interpolation weight
          qw12 = (qL - q11) / (q12 - q11);
	  qw11 = 1 - qw12;
        else
          qw12 = 1;
	  qw11 = 0;
        end

	% interpolate temperature and partial pressure
	kcmp1(:,iL) = kcomp(:,iLr,itlo(iL),iq11)*twlo(iL)*qw11 + ...
                      kcomp(:,iLr,ithi(iL),iq11)*twhi(iL)*qw11 + ...
                      kcomp(:,iLr,itlo(iL),iq12)*twlo(iL)*qw12 + ...
                      kcomp(:,iLr,ithi(iL),iq12)*twhi(iL)*qw12;
      end % H2O test

      % scale interpolated compact absorptions by profile gas amount
      kcmp1(:,iL) = kcmp1(:,iL) .* ...
         (kcprof.gamnt(iL,gind)./refpro.gamnt(iLr,rgind)).^kpow;

    end % layer loop

    % accumulate expanded absorptions
    if kpow == 1/4
      od_gind = ((B * kcmp1).^ 2) .^ 2;  % faster when kpow = 1/4
    else
      od_gind = (B * kcmp1).^(1/kpow);  % the general case
    end

    if (gid == 2)
      iChi = -1;
      if vchunk == 2255
        chi = 'co2_4um_fudge_2255_a.txt';
        iChi = +1;
      elseif vchunk == 2280
        chi = 'co2_4um_fudge_2280_a.txt';
        iChi = +1;
      elseif vchunk == 2380
        chi = 'co2_4um_fudge_2380_b.txt';
        iChi = +1;
      elseif vchunk == 2405
        chi = 'co2_4um_fudge_2405_b.txt';
        iChi = +1;
        end
      if iChi > 0
        fprintf(1,'   CO2 chi function for %5i \n',floor(vchunk));
        chi = fullfile(datadir,['ChiFile/' chi]);
        chi = load(chi);
        chi = chi(:,2); chi = chi*ones(1,length(kcprof.mtemp));
        od_gind = od_gind.*chi;
      end
    end

    % Add current OD to running total
    ibad = find(od_gind < 0);
    od_gind(ibad) = 0;
    absc = absc + od_gind;

    % Add on water continuum
    if (gid == 1)
       [wcon] = contcalc(kcprof, fr,datadir);
       absc = absc + wcon;
    end

  else
      %disp(['Unable to read kcarta data for gas ' int2str(gid)])
  end % valid gas ID and compressed data existance check

end % gas loop

%%% end of function %%%
