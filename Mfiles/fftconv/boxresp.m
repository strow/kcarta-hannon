function r = boxresp(v, L)

% function r = boxresp(v, L)
%
% boxcar response function
%
% inputs
%   v - wavenumbers
%   L - max path length
%
% output
%   r - boxcar response at v
%
% Copyright 2012, Univ. Of Md, Balt. Co. Atmospheric Spectroscopy Laboratory
% kcarta is distributed under the terms of the GNU GPL v3
%


if nargin == 1
  L = 1;
end

r = sinxx(2*pi*L*v);

