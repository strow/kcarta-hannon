
% IASI specs for 1-band spanning interferograms
%
%
% Copyright 2012, Univ. Of Md, Balt. Co. Atmospheric Spectroscopy Laboratory
% kcarta is distributed under the terms of the GNU GPL v3
%

band = 1;

v1 = 605;        % band low end {cm^-1}
v2 = 2829.9975;  % band high end {cm^-1}

L1 = 2.0;        % longest path {cm}

vlaser = 12992;  % laser frequency {cm^-1}
vsf = 2;	 % vlaser scaling factor
dvc = 1/(2*L1);  % channel spacing {cm^-1}

Lcut = L1;	 % path length saved {cm}


