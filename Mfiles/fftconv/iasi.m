
% IASI specs for 1-band spanning interferograms
%
%
% Copyright 2012, Univ. Of Md, Balt. Co. Atmospheric Spectroscopy Laboratory
% kcarta is distributed under the terms of the GNU GPL v3
%

band = 1;

v1 = 605;        % band low end
v2 = 2829.9975;  % band high end

L1 = 2.0;        % longest path

vlaser = 15780;  % gives integer dvc multiples
vsf = 2;	 % vlaser scaling factor
dvc = 1/(2*L1);  % channel spacing

Lcut = L1;	 % path length saved

