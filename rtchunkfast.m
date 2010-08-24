function [rad] = rtchunkfast(lbot,stemp,secang,secsun,xsun,ermono,...
   rplanck,fmono,absc);

% function [rad] = rtchunkfast(lbot,stemp,secang,secsun,xsun,ermono,...
%   rplanck,fmono,absc);
%
% Calculate radiances for a "chunk" of kcarta optical depths using.
% Intended for use with "preprad.m" to prepare terms independent
% of the optical depth.
%
% Input:
%    lbot   - [1 x 1] bottom layer index
%    stemp -  [1 x 1] surface temperature
%    secang - [1 x nlay] local path secant angles
%    secsun - [1 x nlay] round trip sun secant angle
%    xsun   - [npts x 1] cossun*omega*hsun
%    ermono - [npts x 2] emis (1) and rho (2) points
%    rplanck- [npts x nlay] Planck radiance
%    fmono  - [npts x 1] frequency points
%    absc   - [npts x nlay] top-down optical depths
%
% Output:
%   rad     - [npts x 1] radiance at top of atmosphere
% 

% Created: 17 Jun 2009, Scott Hannon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare nadir optical depths
[npts,nlay] = size(absc);

% For debugging
% $$$ if (nlay < lbot)
% $$$   error('lbot > nlay for absc')
% $$$ end
% $$$ if (lbot < nlay)
% $$$   absc=absc(:,nlay);
% $$$ end

% --------------------------------
% Reflected thermal at arccos(3/5)
% --------------------------------
% Transmittance at diffusivity angle
sectherm = 5/3; % sec( acos(3/5) )
tran = exp(-absc .* sectherm);

% Radiance of empty space
fmono = fmono(:);
rthm = bt2rad(fmono, 2.7);

% loop downward over the layers
for ii = 1:lbot
  rthm = rthm.*tran(:,ii) + rplanck(:,ii).*(1 - tran(:,ii));
end
% Note: above rthm is radiance at one representative angle; to convert
% this to the an approximate integral over all angles multiply by pi.

% Upwards reflected component at surface
rthm = rthm .* (1 - ermono(:,1));  % pi*rthm*(1-emis)/pi

% --------------------------------
% Atmospheric and surface emission
% --------------------------------
% Upward radiance at surface
rad = rthm + bt2rad(fmono,stemp).*ermono(:,1);

% Layer transmittances along view angle
tran = exp(-absc .* (ones(npts,1)*secang));

% Loop upward over the layers
for ii = lbot:-1:1
  rad = rad.*tran(:,ii) + rplanck(:,ii).*(1 - tran(:,ii));
end

% ---------------
% Reflected solar
% ---------------
if (secsun(1) > 0)

  % round trip optical depth
  solabs = sum(absc .* (ones(npts,1)*secsun), 2);

  % TOA reflected solar radiance
  rad = rad + xsun .* exp(-solabs) .* ermono(:,2);
end

%%% end of function %%%
