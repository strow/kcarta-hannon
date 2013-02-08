function y = boxapod(d, L)

% function y = boxapod(d, L)
%
% boxcar apodization
%
% inputs
%   d - distance; may be a vector
%   L - max path length
%
% output
%   y - apodization of d
%
% Copyright 2012, Univ. Of Md, Balt. Co. Atmospheric Spectroscopy Laboratory
% kcarta is distributed under the terms of the GNU GPL v3
%


if nargin == 1
  L = 1;
end

y = (abs(d) <= L);

