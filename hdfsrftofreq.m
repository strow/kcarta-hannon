function [fchan, srf, inds, inde] = hdfsrftofreq(fmono, clist, sfile);

% function [fchan, srf, inds, inde] = hdfsrftofreq(fmono, clist, sfile);
%
% Read SRF file "sfile", pull out data for specified channels "clist", and
% return a area normalized SRF on the frequency grid "fmono".
%
% Input:
%    fmono : [npts x 1] monochromatic frequency points
%    clist : [1 x nchan] desired channel IDs
%    sfile : [string] name of HDF SRF file to read
%
% Output:
%    fchan : [nchan x 1] AIRS channel frequency
%    srf   : [npts x nchan] area normalized SRF on fmono grid
%    inds  : [1 x nchan] index in srf of start point
%    inde  : [1 x nchan] index in srf of end point
%

% Created: 16 Jun 2009, Scott Hannon - from Howard Mottelers "sconv2.m"
%
% Copyright 2012, Univ. Of Md, Balt. Co. Atmospheric Spectroscopy Laboratory
% kcarta is distributed under the terms of the GNU GPL v3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath h4toolsV201
addpath srftest

% guarantee input vectors are columns
fmono = fmono(:);
clist = clist(:);

% read the SRF data
[chanid, freq, width, fwgrid, srfval, fattr] = rdhdfsrf(sfile);

% use the subset of channels specified in clist
[id,indc,junk] = intersect(chanid,clist);
chanid = chanid(indc);
freq = freq(indc);
srfval = srfval(indc,:);
width = width(indc);
[nchan, nspts] = size(srfval);

% srfreq is an nchan x nspts array of frequency points for srfval
srfreq = (width * fwgrid) + freq * ones(1, nspts);

% set the output frequency list
fchan = freq;

% check that fmono is uniformly spaced and increasing
dfmono = diff(fmono);
dmax = max(dfmono);
dmin = min(dfmono);
if abs(dmax - dmin) > 10e-9
  error('fmono must be uniformly spaced')
end

% get the step size of fmono
dv = dfmono(1);
if dv <= 0
  error('fmono must be in increasing order')
end

% get the span of fmono
nfmono = length(fmono);
f1 = fmono(1);
f2 = fmono(nfmono);

% initialize output array
srf = zeros(nfmono,nchan);

% loop on SRFs in indc
for jj =1:nchan

  % get the frequency span of the current SRF
  v1 = srfreq(jj, 1);
  v2 = srfreq(jj, nspts);

  % if the SRF is outside fmono, skip this channel
  if v2 <= f1 | f2 <= v1
    fprintf(1, 'WARNING -- SRF id %d outside of input range\n', ...
	    chanid(jj));
    continue
  end

  % if the SRF overlaps fmono, lop it off to fit; 
  % give a warning message if we lop more than dv
  if v1 < f1
    if v1 < f1 - dv
      fprintf(1, 'WARNING -- truncating LHS of SRF id %d\n', ...
	      chanid(jj));
    end
    v1 = f1;
  end
  if f2 < v2
    if f2 + dv < v2
      fprintf(1, 'WARNING -- truncating RHS of SRF id %d\n', ...
	      chanid(jj));
    end
    v2 = f2;
  end

  % fmonod the indices of the current SRF in fmono 
  v1ind = ceil((v1-f1)/dv) + 1;
  v2ind = floor((v2-f1)/dv) + 1;

  % interpolate the SRF to a subinterval of the fmono grid
  s1 = interp1(srfreq(jj,:), srfval(jj,:), fmono(v1ind:v2ind), 'linear');

  % normalize the SRF
  s1 = s1 ./ sum(s1);

  % apply the SRF
  srf(v1ind:v2ind,jj) = s1;
  inds(jj) = v1ind;
  inde(jj) = v2ind;
end

%%% end of function %%%
