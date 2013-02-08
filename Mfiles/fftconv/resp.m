function a = resp(v, L1, atype, aparg);

% function a = resp(v, L1, atype, aparg);
%
% select an apodization function
%
% Copyright 2012, Univ. Of Md, Balt. Co. Atmospheric Spectroscopy Laboratory
% kcarta is distributed under the terms of the GNU GPL v3
%

if nargin == 2
  fprintf(2, 'resp() warning: using default, kaiser-bessel #6\n');
  atype = 'kaiser-bessel';
  aparg = 6;
end

switch atype

  case {'boxcar', 'box'},		a = boxresp(v, L1);
  case {'triangle', 'tri'},		a = triresp(v, L1);
  case {'hamming', 'ham'},		a = hamresp(v, L1);
  case {'kaiser-bessel', 'kb'},		a = kbresp(v, L1, aparg);
  case {'norton-beer', 'nb'},		a = nbresp(v, L1, aparg);
  case {'cosine', 'cos'},		a = cosresp(v, L1);
  case {'beer'},			a = beerresp(v, L1);
  case {'gauss'},			a = gaussresp(v, L1);
  otherwise, error(['unknown apodization: ', atype]);

end

