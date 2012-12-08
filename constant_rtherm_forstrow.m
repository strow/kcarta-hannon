% --------------------------------
% Reflected thermal at arccos(3/5)
% --------------------------------
% Transmittance at diffusivity angle
sectherm = 5/3; % sec( acos(3/5) )
tran = exp(-xabsc .* sectherm);

% Radiance of empty space
rthm = bt2rad(freq, 2.7);

% Declare planck array
rplanck = zeros(npts,lbot);

% loop downward over the layers
for ii = 1:lbot
  rplanck(:,ii) = bt2rad(freq,prof.ptemp(ii,ip));
  rthm = rthm.*tran(:,ii) + rplanck(:,ii).*(1 - tran(:,ii));
end
% Note: above rthm is radiance at one representative angle; to convert
% this to the an approximate integral over all angles multiply by pi.

% Upwards reflected component at surface
rthm = rthm .* (1 - erfine(:,1));  % pi*rthm*(1-emis)/pi

