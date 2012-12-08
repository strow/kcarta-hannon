
% IASI in one band

band = 3;

  v1 = 605;        % band low end
  v2 = 2829.9975;  % band high end
  L1 = 2^15/15799  % L1 (scaled to vlaser)
  vsf = 2;	   % Vlaser scaling factor
  dvc = 1/(2*L1);

Lcut = L1;	   % path length saved

