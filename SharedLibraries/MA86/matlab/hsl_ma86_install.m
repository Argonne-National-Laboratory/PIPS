function hsl_ma86_install(varargin)
%
%   hsl_ma86_install installs hsl_ma86.
%
%   VARARGIN
%
%   hsl_ma86_install() installs hsl_ma86 and its Matlab Interface. It is
%      assumed that the BLAS and LAPACK routines provided by MATLAB
%      are used, and mex is configured to use your preferred compiler.
%      The test example is not run.
%
%   hsl_ma86_install(TEST) installs hsl_ma86 and its Matlab Interface and
%      optionally runs the test example. It is assumed that the BLAS and LAPACK
%      routines provided with the interface are used, and mex is configured to
%      use your preferred compiler.
%      If TEST <= 0, the test example is not run;
%      if TEST > 0, the test example is run and the user can compare the output
%      with that of the file OUT/install.output. Note that there may be very
%      small differences because of the arithmetic differing on different
%      computers.
%
%   hsl_ma86_install(TEST,LIBS) installs hsl_ma86 and its Matlab Interface
%      and optionally runs the test example. It is assumed that mex is
%      configured to use your preferred compiler.
%      If LIBS has the value 'matlab' this is equivalent to the
%        setting LIBS='-lmwlapack -lmwblas'. Use of this option on a 64-bit
%        platform will force use of 64-bit default integers and may impede
%        the performance of the HSL code.
%      Otherwise LIBS should be set to specify which BLAS to link against, and
%        may optionally specify the location of other libraries (eg libf95.a
%        libgcc.a if they are not on a default search path). Typically
%        this will take the form of LIBS='-L/path/to/blas -lblas'.
%
%   hsl_ma86_install(TEST,LIBS,MEXFLAGS) installs hsl_ma86 and its Matlab
%      Interface and optionally runs the test example. The contents of the
%      variable MEXFLAGS is passed to mex as follows.
%        mex $(MEXFLAGS) -c file.F90
%        mex $(MEXFLAGS) $(LIBS) -output foo.mex file.F90
%     If MEXFLAGS is not supplied it assumes the default value of
%       '-largeArrayDims' on a 64-bit platform and is empty on a 32-bit
%       platform.
%     If MEXFLAGS is supplied and BLAS='matlab' on a 64-bit machine, the user
%       must ensure that the relevant flag to force 64-bit default integers is
%       passed to the compiler by explictly setting FFLAGS='-i8 \$FFLAGS' (g95)
%       or '-fdefault-real-8 \$FFLAGS' (gfortran). (The \$FFLAGS is needed to
%       pick up extra options such as -fPIC from the mex configuration file)
%
%
%   See also  PATHTOOL, JAVAADDPATH, ADDPATH.
%   Also see DOC STARTUP.
%

% srcdir is the directory containing the Fortran source code for the package
srcdir = '../src';

nin  = nargin;
if (nin > 3),
   error('Too many arguments');
end;

v = getversion ;

%-------------------------------------------------------------------------------
% Check that Matlab version is recent enough
%-------------------------------------------------------------------------------

if (v(1) < 7 || v(2) <= 3)
	error('Matlab version too old - require version 7.4 (2007a) or newer');
end
if (v(2) < 6)
   warning('Old matlab version - we only test against 7.6 (2008a) or newer');
end


%-------------------------------------------------------------------------------
% Check that we are on a unix/linux machine
%-------------------------------------------------------------------------------

switch(lower(computer))
    case {'glnx86', 'glnxa64'}
        islinux = true;
        ismac = false;
    case 'maci64'
        islinux = false;
        ismac = true;
    otherwise
        error('Unsupported operating system');
end

% Determine machine architecture
is64 = (~isempty(strfind(computer,'64')));

%-------------------------------------------------------------------------------
% Initialise compile parameters
%-------------------------------------------------------------------------------


% Pick between g95 and gfortran options
if(v(2) <= 11 && islinux)
   % LINUX g95 for 2010b and earlier
   di8 = ' FFLAGS="-i8 \$FFLAGS"'; % option for 64-bit default integer
   g95_fix = true;
else
   % LINUX gfortran for 2011a and above
   % MAC gfortran always
   di8 = ' FFLAGS="-fdefault-integer-8 \$FFLAGS"'; % option for 64-bit d.i.
   g95_fix = false;
end

% Set defaults
LIBS = 'matlab';
if(is64)
   MEXFLAGS = '-largeArrayDims';
else
   MEXFLAGS = '';
end

% Process optional arguments
% LIBS
if (nin>=2)
   LIBS = varargin{2};
end

% Use matlab blas if requested. On 64-bit the requires 64-bit default integer.
% Note: MEXFLAGS can override di8 setting!
if(size(strfind(strtrim(LIBS),'matlab')==0))
   LIBS = '-lmwlapack -lmwblas';
   nometis = 1; % true
   if(is64)
      MEXFLAGS = strcat(MEXFLAGS, di8);
   end
else
   nometis = 0; % false
end

% MEXFLAGS
if (nin>=3)
   MEXFLAGS = varargin{3};
   g95_fix = false;
end

%-------------------------------------------------------------------------------
% Clean directory, then compile in good order
%-------------------------------------------------------------------------------
kk = do_unix_cmd('rm  *.o *.mod *.mexglx *.mexa64 2>/dev/null');
%mexcompile(MEXFLAGS, '-c', [srcdir,'/ddeps.f'], '')
%mexcompile(MEXFLAGS, '-c', [srcdir,'/zdeps.f'], '')
if(g95_fix)
   % g95 has a bug in it with iolength and -i8.
   % common90.g95.f90 contains an alternate version of hsl_zb01 that
   % works around this bug (in a manner that is not sadly not compatible
   % with other compilers!)
   do_unix_cmd('cp common90.g95.f90 common90.f90');
   mexcompile(MEXFLAGS, '-c', 'common90.f90', '')
   do_unix_cmd('rm common90.f90');
else
   mexcompile(MEXFLAGS, '-c', [srcdir,'/common90.f90'], '')
end
mexcompile(MEXFLAGS, '-c', [srcdir,'/common.f'], '')
mexcompile(MEXFLAGS, '-c', [srcdir,'/ddeps.f'], '')
mexcompile(MEXFLAGS, '-c', [srcdir,'/zdeps.f'], '')
mexcompile(MEXFLAGS, '-c', [srcdir,'/ddeps90.f90'], '')
mexcompile(MEXFLAGS, '-c', [srcdir,'/zdeps90.f90'], '')
mexcompile(MEXFLAGS, '-c', [srcdir,'/hsl_ma86d.f90'], '')
mexcompile(MEXFLAGS, '-c', [srcdir,'/hsl_ma86z.f90'], '')
mexcompile(MEXFLAGS, '-c', 'hsl_matlab.F90', '')
if(nometis)
   mexcompile(MEXFLAGS, '-c', [srcdir,'/fakemetis.f'], '')
   mexcompile(MEXFLAGS, '-output hsl_ma86_expert', 'hsl_ma86_expert.f90 hsl_matlab.o hsl_ma86d.o common90.o common.o ddeps.o ddeps90.o hsl_ma86z.o zdeps.o zdeps90.o fakemetis.o', LIBS)
else
   mexcompile(MEXFLAGS, '-output hsl_ma86_expert', 'hsl_ma86_expert.f90 hsl_matlab.o hsl_ma86d.o common90.o common.o ddeps.o ddeps90.o hsl_ma86z.o zdeps.o zdeps90.o', LIBS)
end

%-------------------------------------------------------------------------------
% Add to path and print nice message
%-------------------------------------------------------------------------------

currentpath = which (mfilename);
sep = find (currentpath == filesep);
currentpath = currentpath (1:sep(end));

fprintf ('Temporarily adding %s to your MATLAB path and JAVA path.\n', currentpath) ;
fprintf ('Do this permanently via pathtool.  Next, edit the file:\n') ;
fprintf ('%s file.\n', which ('classpath.txt')) ;
fprintf ('and add the line:\n') ;
fprintf ('%s\n', currentpath) ;
fprintf ('to the end of that file (which defines your JAVA class path).\n') ;

addpath (currentpath) ;
javaaddpath (currentpath) ;

fprintf ('\nAlternatively, add these two lines to your startup.m file:\n\n') ;
fprintf ('addpath (''%s'') ;\n', currentpath) ;
fprintf ('javaaddpath (''%s'') ;\n', currentpath) ;
fprintf ('\nSee also pathtool, javaaddpath, addpath, and "doc startup".\n\n') ;
fprintf ('Please cite HSL as\n\t"HSL, a collection of Fortran codes for large-scale scientific\n\t computation. See http://www.hsl.rl.ac.uk/"\n\n');
fprintf ('Please report any bugs to hsl@stfc.ac.uk\n\n');

if (nin >= 1 && varargin{1}>0),
   hsl_ma86_test;
end;
clear nin v currentpath sep mexpath MEXpath

%-------------------------------------------------------------------------------

function v = getversion
%GETVERSION determine the MATLAB version, and return it as a double.
v = sscanf (version, '%d.%d.%d') ;

%-------------------------------------------------------------------------------

function kk = do_unix_cmd(s)
% DO_UNIX_CMD executes the unix command that is supplied within the string s.
fprintf('%s\n',s);
kk = unix(s);

function mexcompile(MEXFLAGS, output, files, LIBS)
% Builds and executes string 'mex $(MEXFLAGS) $(LIBS) $(output) $(files)'
s = ['mex ', MEXFLAGS, ' ', LIBS, ' ', output, ' ', files];
fprintf('%s\n', s);
eval(s);
