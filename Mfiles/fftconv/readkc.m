function [data, wnums] = readkc(kfile, dfile)

% function readkc(kfile, dfile)
%
% readkc is a simple reader & unchunker for kcarta output files
%
% INPUTS
%
%   kfile  - kcarta output file
%   dfile  - readkc output file (optional)
%
% OUTPUTS
%
%   data   - a w by n array of data from kcarta
%   wnums  - a w by 1 vector of data wavenumbers
%
% If the input parameter dfile is specified, then the data array
% is written to file dfile, and the return values [data, wnums]
% are left unassigned.
%
% The w by n data array has one row for each output wavenumber and
% one column for each "output data block" row.  In practice, this
% means the data array columns are a concatenation of whatever was 
% asked for in the kcarta driver file.
% 
% LIMITATIONS
% 
% readkc() only works with the "short" kcarta header format.
%
% This version doesn't work with kcarta jacobian files, as these
% use a different header and slightly different output data blocks
% 
% Large data sets (e.g., mixed path sets of more than one or two
% "chunks") take up too much space to be returned in memory; for 
% such cases an output file should always be specified.

% H. Motteler, 8/24/98


[fin,msg] = fopen(kfile, 'r');
if fin == -1
  error(['error opening input file\n', msg]);
  end

%%%%%%%%%%%%%%%%%%%%%%
% READ HEADER RECORDS
%%%%%%%%%%%%%%%%%%%%%%

% MAIN HEADER

% version number
flen    = fread(fin, 1, 'integer*4');
version = fread(fin, 80, 'char');
version = setstr(version');
flen    = fread(fin, 1, 'integer*4');

% number of layers
flen   = fread(fin, 1, 'integer*4');
nlayer = fread(fin, 1, 'integer*4');
flen   = fread(fin, 1, 'integer*4');

% number of params 
flen    = fread(fin, 1, 'integer*4');
nparams = fread(fin, 1, 'integer*4');
flen    = fread(fin, 1, 'integer*4');
flen    = fread(fin, 1, 'integer*4');
rparams = fread(fin, nparams, 'real*4');
flen    = fread(fin, 1, 'integer*4');

% comment
flen    = fread(fin, 1, 'integer*4');
comment = fread(fin, 80, 'char');
comment = setstr(comment');
flen    = fread(fin, 1, 'integer*4');

% start, stop frequency
flen = fread(fin, 1, 'integer*4');
fmin = fread(fin, 1, 'real*4');
fmax = fread(fin, 1, 'real*4');
flen = fread(fin, 1, 'integer*4');

% low and high chunk index
flen      = fread(fin, 1, 'integer*4');
lowchunk  = fread(fin, 1, 'integer*4');
highchunk = fread(fin, 1, 'integer*4');
flen      = fread(fin, 1, 'integer*4');

% chunk size
% flen      = fread(fin, 1, 'integer*4');    
% chunksize = fread(fin, 1, 'integer*4');
% flen      = fread(fin, 1, 'integer*4');

% read iLongOrShort
flen  = fread(fin, 1, 'integer*4');
htype = fread(fin, 1, 'integer*4');
flen  = fread(fin, 1, 'integer*4');

if htype ~= -1
  error('readkc: can only do short headers');
  end

% if chunksize ~= 25
%   error('readkc: can only do 25 cm-1 chunks');
%   end

% GAS PATH HEADER

% number of paths to be output
flen    = fread(fin, 1, 'integer*4');
ngasout = fread(fin, 1, 'integer*4');
flen    = fread(fin, 1, 'integer*4');

if (ngasout > 0)
  flen     = fread(fin, 1, 'integer*4');
  gaspaths = fread(fin, ngasout, 'integer*4');
  flen     = fread(fin, 1, 'integer*4');
  end

% MIXED PATH HEADER

flen      = fread(fin, 1, 'integer*4');
nmixpaths = fread(fin, 1, 'integer*4');
flen      = fread(fin, 1, 'integer*4');

if (nmixpaths > 0)

  flen    = fread(fin, 1, 'integer*4');
  nmixout = fread(fin, 1, 'integer*4');
  flen    = fread(fin, 1, 'integer*4');

  if (nmixout > 0)
    flen     = fread(fin, 1, 'integer*4');
    mixpaths = fread(fin, nmixout, 'integer*4');
    flen     = fread(fin, 1, 'integer*4');
    end
  end

% ATMOSPHERE HEADER

% number of atmospheres in *RADFIL
flen   = fread(fin, 1, 'integer*4');
natmos = fread(fin, 1, 'integer*4');
flen   = fread(fin, 1, 'integer*4');

% max number of emissivity points per atm
flen  = fread(fin, 1, 'integer*4');
nemis = fread(fin, 1, 'integer*4');
flen  = fread(fin, 1, 'integer*4');

if (natmos > 0)

  for i = 1:natmos

    % read in atmosphere number, # of mixed paths in atmosphere
    flen = fread(fin, 1, 'integer*4');
    iAtm = fread(fin, 1, 'integer*4');
    iNumPathsAtmos = fread(fin, 1, 'integer*4');
    flen = fread(fin, 1, 'integer*4');

    if (iNumPathsAtmos > 0)

      %read in paths making up atmosphere
      flen = fread(fin, 1, 'integer*4');
      ia   = fread(fin, iNumPathsAtmos, 'integer*4');
      flen = fread(fin, 1, 'integer*4');

      flen   = fread(fin, 1, 'integer*4');
      rTbdy  = fread(fin, 1, 'real*4');
      rTinit = fread(fin, 1, 'real*4');
      rSat   = fread(fin, 1, 'real*4');
      rHgt   = fread(fin, 1, 'real*4');
      flen   = fread(fin, 1, 'integer*4');

      flen           = fread(fin, 1, 'integer*4');
      ikSolar        = fread(fin, 1, 'integer*4');
      rkSolarAngle   = fread(fin, 1, 'real*4');
      rkSolarRefl    = fread(fin, 1, 'real*4');
      ikThermal      = fread(fin, 1, 'integer*4');
      rkThermalAngle = fread(fin, 1, 'real*4');
      ikThermalJacob = fread(fin, 1, 'integer*4');
      flen           = fread(fin, 1, 'integer*4');

      flen = fread(fin, 1, 'integer*4');
      iEms = fread(fin, 1, 'integer*4');
      flen = fread(fin, 1, 'integer*4');

      %now read the freq, ems data points
      for jj=1:iEms
        flen = fread(fin, 1, 'integer*4');
        rff  = fread(fin, 1, 'real*4');
        ree  = fread(fin, 1, 'real*4');
        flen = fread(fin, 1, 'integer*4');
        end 

      %now read in layers to be output
      flen          = fread(fin, 1, 'integer*4');
      iNumLayersOut = fread(fin, 1, 'integer*4');
      flen          = fread(fin, 1, 'integer*4');

      if (iNumLayersOut > 0)
        % read in mixed paths to be output
        flen = fread(fin, 1, 'integer*4');
        ia   = fread(fin, iNumLayersOut, 'integer*4');
        flen = fread(fin, 1, 'integer*4');
        end

      if (iNumLayersOut > 0)
        % read in pressures at which radiances to be output
        flen = fread(fin, 1, 'integer*4');
        ra   = fread(fin, iNumLayersOut, 'real*4');
        flen = fread(fin, 1, 'integer*4');
        end
      end
    end
  end

%%%%%%%%%%%%%%%%%%%%%%
% DATA RECORD SUMMARY
%%%%%%%%%%%%%%%%%%%%%%

flen   = fread(fin, 1, 'integer*4');
nchunk = fread(fin, 1, 'integer*4'); % total number of chunks
nODBs  = fread(fin, 1, 'integer*4'); % number of (ODBs) "output data blocks"
flen   = fread(fin, 1, 'integer*4');

flen1    = fread(fin, 1, 'integer*4');
nODBrows = fread(fin, nODBs, 'integer*4'); % rows in each ODB
flen2    = fread(fin, 1, 'integer*4');

% sanity checks
if flen1 ~= flen2
  error('Fortran records out of phase!');
end

if nchunk ~= 1+highchunk-lowchunk
  error('readkc: chunk counts do not match!');
end

fprintf(2, 'readkc: %d chunks, %d ODBs, %d total rows\n', ...
	nchunk, nODBs, sum(nODBrows));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize output file or array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noutrow = nchunk * 10000;	% number of output rows
noutcol = sum(nODBrows);	% number of output columns
ODBrowdat = zeros(10000, 1);	% ODB row data buffer

if nargin == 1
 
  % initialize output arrays
  data = zeros(noutrow, noutcol);
  wnums = zeros(noutrow, 1);

else

  % pre-extend the output file
  totbytes = noutrow * noutcol * 4;
  fprintf(2, 'readkc: creating %d x %d element output array\n', ...
	  noutrow, noutcol);

  [fout,msg] = fopen(dfile, 'w');
  if fout == -1
    error(['error creating output file', msg]);
  end
  z = zeros(noutrow,1);
  for i = 1:noutcol
    if fwrite(fout, z, 'float') ~= noutrow
      error(['fwrite in extend failed, i=',num2str(i)]);
    end
  end
  fclose(fout);
  clear z

  % open the extended file in "update" mode
  [fout,msg] = fopen(dfile, 'r+');
  if fout == -1
    error(['error opening output file', msg]);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read "output data blocks"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for chunk = 1:nchunk

  cumODBrow = 1;

  for odb = 1:nODBs

    % read ODB header
    flen    = fread(fin, 1, 'integer*4');
    type    = fread(fin, 1, 'integer*4'); % ODB type
    subtype = fread(fin, 1, 'integer*4'); % ODB subtype
    nrow    = fread(fin, 1, 'integer*4'); % number of rows, this ODB
    flen    = fread(fin, 1, 'integer*4');

    flen    = fread(fin, 1, 'integer*4');
    ncols   = fread(fin, 1, 'integer*4'); % spectral pts in a chunk
    frlow   = fread(fin, 1, 'real*4');    % low freq of chunk
    frhigh  = fread(fin, 1, 'real*4');    % high freq of chunk
    frinc   = fread(fin, 1, 'real*4');    % frequency increment
    flen    = fread(fin, 1, 'integer*4');
  
    if chunk == 1
      fprintf(2, 'readkc: ODB type = %d, subtype = %d, rows = %d\n', ...
		type, subtype, nrow);
    end

    % sanity check
    if nrow ~= nODBrows(odb)
      fprintf(2, 'readkc: WARNING -- bad ODB row count,  ');
      fprintf(2, 'header = %-4d ODB = %-4d\n', nODBrows(odb), nrow);
    end

    % read a row of ODB data
    for ODBrow = 1:nODBrows(odb)
      flen1 = fread(fin, 1, 'integer*4');
      [ODBrowdat, count] = fread(fin, [10000,1], 'real*4');
      if count ~= 10000
	error(['fread failed, odb=',num2str(odb),' chunk=',num2str(chunk)]);
      end
      flen2 = fread(fin, 1, 'integer*4');

      if nargin == 1

        % write ODBrowdat to the data array
        outrows = (chunk-1) * 10000 + 1 : chunk * 10000;
        outcol = cumODBrow;
        data(outrows, outcol) = ODBrowdat;
        wnums(outrows) = frlow + (0:9999)*frinc;

      else

        % write ODBrowdat to the output file
        outpos = ((cumODBrow-1)*nchunk*10000 + (chunk-1)*10000)*4;
        if fseek(fout, outpos, -1) ~= 0
          error(['fseek failed, odb=',num2str(odb),' chunk=',num2str(chunk)]);
        end
        if fwrite(fout, ODBrowdat, 'real*4') ~= 10000
          error(['fwrite failed, odb=',num2str(odb),' chunk=',num2str(chunk)]);
        end     

      end

      cumODBrow = cumODBrow + 1;

    end % loop on current ODB rows

  end % loop on ODBs

end % loop on chunks

fclose(fin);

if nargin == 2 
  fclose (fout);
end

