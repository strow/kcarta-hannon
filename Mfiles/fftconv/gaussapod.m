function y = gaussapod(d, L)

% function y = gaussapod(d, L)
%
% gassian apodization
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

c =  2* sqrt(-log(.4107)) / L;

a = d * c;

y = (0 <= d & d <= L) .* exp(-(a.*a));

