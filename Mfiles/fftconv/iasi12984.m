
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

vlaser = 12984;  % laser frequency {cm^-1}
% vlaser is approximately 1.54 um but exact value is locked to
% some unspecified acetylene line.  The strongest line (at 296K)
% near 1.54 um is P23e at 1.53943 um = 6495.9119 cm^-1
% Code complains v2 > vmax with vlaser=6496, so double it to 12992
% Note that with vlaser=6496 cm^-1 and dvc=0.25 cm, the number
% of laser wavelengths in dvc is an exact integer
%    num_wavelengths = 1E-2 m/cm * dvc / ( 1E-6 m/um * 1/(vlaser/10000) )
%                    = 0.25 * 6496 = 1624 exactly
% This is true of any vlaser divisable by 4 (ie 1/0.25)
% 
vsf = 2;	 % vlaser scaling factor
dvc = 1/(2*L1);  % channel spacing {cm^-1}

Lcut = L1;	 % path length saved {cm}


