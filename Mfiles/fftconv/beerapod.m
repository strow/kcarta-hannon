function b = beerapod(d, L)

% function b = beerapod(d, L)
%
% Beer apodization
%
% inputs
%   d - distance; may be a vector
%   L - max path length
%
% output
%   b - apodization of d
% 
% NOTE: from the Wisconson group, via Larrabee
%
% Copyright 2012, Univ. Of Md, Balt. Co. Atmospheric Spectroscopy Laboratory
% kcarta is distributed under the terms of the GNU GPL v3
%

b = (abs(d) <= L) .* (1 - (d/L).^2).^2;

