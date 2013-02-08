function c = hamapod(d, L)

% function c = hamapod(d, L)
%
% Hamming apodization
%
% inputs
%   d - distance; may be a vector
%   L - max path length
%
% output
%   c - apodization of d
%
% Copyright 2012, Univ. Of Md, Balt. Co. Atmospheric Spectroscopy Laboratory
% kcarta is distributed under the terms of the GNU GPL v3
%

if nargin == 1
  L = 1;
end

a = 0.23;

c = (abs(d) <= L) .* ((1 - 2*a) + 2*a*cos(pi*d/L));

