function [kcprof] = op_rtp_to_lbl(profnum, gasid, head, prof);

% function [kcprof] = op_rtp_to_lbl(profnum, gasid, head, prof);
%
% Pull out data for the specified profile and gas ID from a "layers"
% RTP profile structure, and return it in a new structure for use with
% the UMBC LBL and KCMIX codes.
%
% Input:
%    profnum  : [1 x 1] index of desired profile in "prof"
%    gasid    : [1 x ngas] gas ID
%    head     : RTP head structure
%    prof     : RTP prof structure
%
% Output:
%    kcprof : structure with the following fields:
%       glist : [1 x ngas] gas ID, same as "gasid"
%       mtemp : [nlays x 1] layer mean temperature {Kelvin}
%       mpres : [nlays x 1] layer pressure {atm}
%       gpart : [nlays x ngas] gas partial pressure {atm}
%       gamnt : [nlays x ngas] integrated gas amount {kilomoles/cm^2}
%

% Created 23 April 2002, Scott Hannon - based on "op_rtp_to_txt.m" for KCARTA
% Update: 06 Dec 2005, Sergio Machado - created non-txt version
% Update: 12 Jun 2009, S.Hannon - minor update and cleanup for rtpV201
% Update: 15 Jun 2009, S.Hannon - allow multiple gases in glist
%
% Copyright 2012, Univ. Of Md, Balt. Co. Atmospheric Spectroscopy Laboratory
% kcarta is distributed under the terms of the GNU GPL v3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check ptype
if (~isfield(head,'ptype'))
   error('head field ptype not found')
else
   if (head.ptype < 1)
      error('head field ptype must be a "layers" profile (ptype=1 or 2)')
   end
end

% Check glist
if (~isfield(head,'glist'))
   error('head field glist not found')
else
   [igas,indglist,indgasid] = intersect(head.glist,gasid);
   ngas = length(igas);
   if (length(gasid) ~= ngas)
      error('Unable to match all gases in gasid to head.glist')
   end
end

% Check gunit
if (~isfield(head,'gunit'))
   error('head field gunit not found')
else
   ii = unique( head.gunit(indglist) );
   if (length(ii) ~= 1 | ii(1) ~= 1)
      error('head.gunit must be 1 (molecules/cm^2) for all gases in gasid')
   end
end

% Check nlevs
if (~isfield(prof,'nlevs'))
   error('prof field nlevs not found')
end
nlevs=prof.nlevs(profnum);
nlays=nlevs - 1;

% Layer mean pressure {mb}
if (~isfield(prof,'plevs'))
   error('prof field plevs not found')
end
plevs = prof.plevs(1:nlevs,profnum);
ind1 = 1:nlays;
ind2 = 2:nlevs;
plays = (plevs(ind2)-plevs(ind1))./log(plevs(ind2)./plevs(ind1));

% Layer thickness {meters}
if (~isfield(prof,'palts'))
   error('prof field palts not found')
end
zhi = prof.palts(ind1,profnum);
zlo = prof.palts(ind2,profnum);
dz = abs(zhi - zlo);

% Layer mean temperature {Kelvin}
if (~isfield(prof,'ptemp'))
   error('prof field ptemps not found')
end
temp = prof.ptemp(ind1,profnum);

% Create output structure
kcprof.glist = gasid;
kcprof.mtemp = temp;
kcprof.mpres = plays./1013.25; % convert from mb to atm
kcprof.gpart = zeros(nlays,ngas);
kcprof.gamnt = zeros(nlays,ngas);

% Loop over gases
for ii=1:ngas
   gasstr=['gas_' int2str(gasid(ii))];
   if (isfield(prof,gasstr))

      % Get layer gas amount
      eval(['amount=prof.' gasstr '(ind1,profnum);']);

      % Calc partial pressure {atm}
      pp = amount .* temp ./ ( 2.6867775E+19 * 273.15 * 100*dz );

      kcprof.gpart(:,ii) = pp;
      kcprof.gamnt(:,ii) = amount./6.02214199E+26; % molecules to kilomoles

   else
      error(['prof does not contain field ' gasstr])
   end
end
%%% end of function %%%
