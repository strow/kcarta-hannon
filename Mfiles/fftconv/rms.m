function r = rms(A)
%
% Copyright 2012, Univ. Of Md, Balt. Co. Atmospheric Spectroscopy Laboratory
% kcarta is distributed under the terms of the GNU GPL v3
%

% function r = rms(A)
%
% returns RMS average of vector or matrix A

[m,n] = size(A);

r = norm(A(:)) / sqrt(m*n);

