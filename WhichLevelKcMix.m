function [iB] = WhichLevelKcMix(plev,nlevs,p)
%this is to find at which pressure level boundary a pressure "p" lies at

% Copyright 2012, Univ. Of Md, Balt. Co. Atmospheric Spectroscopy Laboratory
% kcarta is distributed under the terms of the GNU GPL v3


iB = nlevs;
if (plev(1) > plev(iB))
  plev = flipud(plev);
  end

while ((p <= plev(iB)) & (iB >= 1))
  iB = iB - 1;
  end
