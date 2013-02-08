function a = sinxx(x)
%
% Copyright 2012, Univ. Of Md, Balt. Co. Atmospheric Spectroscopy Laboratory
% kcarta is distributed under the terms of the GNU GPL v3
%

% function a = sinxx(x)
% 
% returns sin(x)/x, 1 at x = 0

d = x==0;
a = sin(x) ./ (x+d);
a = a .* (1-d) + d;

