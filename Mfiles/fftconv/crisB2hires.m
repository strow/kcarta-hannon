
% CrIS specs for band 2
%
%
% Copyright 2012, Univ. Of Md, Balt. Co. Atmospheric Spectroscopy Laboratory
% kcarta is distributed under the terms of the GNU GPL v3
%

band = 2;

v1 = 1200;       % band low end
v2 = 1800;       % band high end

L1 = 0.8;        % longest path
dvc = 1/(2*L1);  % channel spacing

vlaser = 15780;  % gives integer dvc multiples
vsf = 2;	 % vlaser scaling factor

Lcut = L1;	 % path length saved

