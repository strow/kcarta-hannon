% FASTKCMIX - last updated 17 Jun 2009, Scott Hannon
% WARNING! not yet been tested for accuracy/consistency
%
% Overview:
%    Compute optical depth and radiance for some "chunk"
%    of the kCARTA database. Most of these routines are derived from
%    Howard Motteler's kcmix code (circa 2002).  The code has been
%    simplified to increase speed, but at some loss in flexibility.
%    In particular, this code can only be used with RTP profiles on
%    the 100 AIRS layers. The klayers program should be run for ALL
%    required gases in the chunk of interest (see kCARTA database).
%
% Main routines:
%   kcmixfast.m - compute optical depth using kCARTA for one "chunk"
%   rtchunk.m - stand alone radiative transfer calculation
%   preprad.m - prepare for a radiance calculation by computing terms that
%      are independent of the atmospheric transmittance.  Intended for use
%      with "rtchunkfast.m".
%   rtchunkfast.m - raditive transfer calculation for use with "preprad.m"
%
% Support routines:
%   hdfsrftofreq.m - Read AIRS SRF file and pull out SRFs for selected channels
%   op_rtp_to_lbl.m - convert RTP profile for use with kCARTA
%   sunang_conv.m - convert surface zenith angle to local path angle
%   vaconv.m - convert satellite view angle to local path angle
%
% Low level routines:
%   contcalc.m - calculate water continuum optical depth
%   contread.m - read water continuum kCARTA database file
%   rdgaschunk_le.m - read little endian FORTRAN binary kCARTA database file
%
% Example programs:
%   Test/test_kcarta_norm.m
%   Test/test_kcarta_preprad.m
%
% Data files for example programs:
%      Data/desert_0725_2004.op.rtp
%
%%% end of file %%%
