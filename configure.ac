#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#    This file is part of ICTP RegCM.
#
#    ICTP RegCM is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ICTP RegCM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

AC_INIT([RegCM], [4.3.5.1], [ggiulian@ictp.it])

AC_SUBST([PACKAGE_URL], ["https://gforge.ictp.it/gf/project/regcm"])

#AC_PREREQ([2.64])

LT_INIT

AM_INIT_AUTOMAKE([foreign])

AC_PROG_RANLIB
AC_PROG_CC
AC_PROG_FC

AC_LANG([Fortran])
AC_OPENMP

AX_PROG_NC_CONFIG([NC_PREFIX=`nc-config --prefix` && NCAUTO="yes"],
                  [NC_PREFIX="${NETCDF}" && NCAUTO="no"])
AC_ARG_WITH([netcdf],
  AC_HELP_STRING([--with-netcdf],
            [Path to NetCDF installation (default: NETCDF environment)]),
  [NC_PREFIX="$withval" && NCAUTO="no"], [FAKE="no"]
)

AC_ARG_WITH([hdf5],
  AC_HELP_STRING([--with-hdf5],
            [Path to HDF5 installation (default: HDF5 environment)]),
  [HDF5_PREFIX="$withval"], [HDF5_PREFIX="${HDF5}"]
)

AC_ARG_WITH([szip],
  AC_HELP_STRING([--with-szip],
            [Path to SZIP installation (default: SZIP environment)]),
  [SZIP_PREFIX="$withval"], [SZIP_PREFIX="${SZIP}"]
)

PARALLEL="1"
AC_ARG_ENABLE(mpiserial,
  AC_HELP_STRING([--enable-mpiserial],
                 [Use the included MPI replacement library for single processor]),
  [DONOTCHECKMPI="DONOTCHECKMPI"; PARALLEL="0"]
)

USECPL="0"
#AM_CONDITIONAL([DO_COMPILE_LIBESMF], [test x"$USECPL" = x1])
AC_ARG_WITH([esmf],
  AC_HELP_STRING([--with-esmf=<yes|no|directory> (default=no)],
            [Path to ESMF installation (default: ESMF_LIB environment)]),
            [AS_CASE(["$with_esmf"],
                     [no],[AC_MSG_CHECKING([for ESMF library])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_CHECK_FILE([$ESMF_LIB/esmf.mk],
                           [AS_IF([test x"$PARALLEL" = x1],
                                  [ESMF_PREFIX="${ESMF_LIB}"
                                   USECPL="1"
                                   AM_CPPFLAGS="-DESMFCPL $AM_CPPFLAGS"],
                                  [AC_MSG_NOTICE([ESMF enabled code works only with MPI version])])])],
                     [*],[ESMF_PREFIX=$with_esmf
                          AS_IF([test -d "$ESMF_PREFIX"],
                                [AS_IF([test x"$PARALLEL" = x1],
                                       [USECPL="1" AM_CPPFLAGS="-DESMFCPL $AM_CPPFLAGS"],
                                       [AC_MSG_NOTICE([ESMF enabled code works only with MPI version])])],
                                [AC_MSG_NOTICE([$ESMF_PREFIX is not a directory! ESMF suppressed])])],
            [AC_MSG_CHECKING([for ESMF library])
             AC_MSG_RESULT([suppressed])]
  )]
)
AM_CONDITIONAL([DO_COMPILE_LIBESMF], [test x"$USECPL" = x1])
if eval "test x$USECPL = x1"; then
  AC_MSG_WARN([Compiling CPL enabled code using ESMF library])
  AC_SUBST(ESMF_PREFIX)
fi

AC_ARG_WITH([roms],
  AC_HELP_STRING([--with-roms=<no|directory> (default=no)],
            [Path to ROMS installation (default: none)]),
            [AS_CASE(["$with_roms"],
                     [no],[AC_MSG_CHECKING([for ROMS coupling])
                           AC_MSG_RESULT([suppressed])],
                     [yes],[AC_MSG_NOTICE([The ROMS installation directory must be given! ROMS coupling suppressed])],
                     [*],[ROMS_PREFIX=$with_roms
                          AS_IF([test -d "$ROMS_PREFIX"],
                                [AS_IF([test x"$PARALLEL" = x1 & test x"$USECPL" = x1],
                                       [USEROMS="1"],
                                       [AC_MSG_NOTICE([ROMS coupling works only with MPI and ESMF version])])],
                                [AC_MSG_NOTICE([$ROMS_PREFIX is not a directory! ROMS coupling suppressed])])],
            [AC_MSG_CHECKING([for ROMS coupling])
             AC_MSG_RESULT([suppressed])]
  )]
)
AM_CONDITIONAL([DO_COMPILE_ROMS], [test x"$USEROMS" = x1])
if eval "test x$USEROMS = x1"; then
  AC_CHECK_FILE([$ROMS_PREFIX/seaice.o],
                [AM_CPPFLAGS="-DROMSICE $AM_CPPFLAGS"; AC_MSG_WARN([Compiling with ROMS-ICE model])],
                AC_MSG_WARN([Compiling with ROMS model]))
  AC_SUBST(ROMS_PREFIX)
fi

AC_ARG_ENABLE(selectedrealkind,
  AC_HELP_STRING([--enable-selectedrealkind],
                 [define single precision and double precision numbers using the selected_real_kind function. Default uses the kind inquiry function.]),
                 [SRKDEF="SELECTEDREALKIND"]
)

AC_ARG_ENABLE(sequence,
  AC_HELP_STRING([--enable-sequence],
                 [Modify MCT types to make them contiguous in memory.]),
  [SRKDEF="SEQUENCE"],
)

USECLM="0"
USEMEGAN="0"
AC_ARG_ENABLE(clm,
  AC_HELP_STRING([--enable-clm],
                 [Supply this option if you plan on using CLM option.]),
  [USECLM="1"],
)
AM_CONDITIONAL(DO_COMPILE_LIBCLM, [test x"$USECLM" = x1])

AC_ARG_ENABLE(megan,
  AC_HELP_STRING([--enable-megan],
                 [Supply this option if you plan on using CLM option with MEGAN.]),
  [USEMEGAN="1"],
)

if eval "test x"$USECLM" = x1"; then
  AC_MSG_WARN([Compiling CLM enabled code])
  AM_CPPFLAGS="-DCLM $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
  if eval "test x"$USEMEGAN" = x1"; then
    AC_MSG_WARN([Compiling CLM with MEGAN enabled code])
    AM_CPPFLAGS="-DVOC $AM_CPPFLAGS"
    AC_SUBST(AM_CPPFLAGS)
  fi
fi

AC_ARG_ENABLE(nc4-gzip,
  AC_HELP_STRING([--enable-nc4-gzip],
                 [Supply this option if you plan on using NetCDF4/HDF5 output compression option.]),
  [NETCDF4_COMPRESS="1"],[NETCDF4_COMPRESS="0"]
)
if eval "test x$NETCDF4_COMPRESS = x1"; then
  AM_CPPFLAGS="-DNETCDF4_COMPRESS $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

AC_ARG_ENABLE(nc4-parallel,
  AC_HELP_STRING([--enable-nc4-parallel],
                 [Supply this option if you plan on using NetCDF4/HDF5 mpi output parallel output option.]),
  [NETCDF4_HDF5="1"],[NETCDF4_HDF5="0"]
)
if eval "test x$NETCDF4_HDF5 = x1"; then
  AM_CPPFLAGS="-DNETCDF4_HDF5 $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

DEBUG="0"
AC_ARG_ENABLE(debug,
  AC_HELP_STRING([--enable-debug],
                 [Enable debugging flags and per processor log file]),
  [DEBUG="1"]
)
if eval "test x$DEBUG = x1"; then
  AM_CPPFLAGS="-DDEBUG $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

if eval "test x$NCAUTO = xno"; then
  NC_INCLUDES="-I$NC_PREFIX/include"
else
  NCVER=`nc-config --version | cut -d " " -f 2`
  if eval "test x$NCVER = x4.0.1"; then
    NC_INCLUDES="`nc-config --fflags`"
  else
    # new netcdf library has this
    NC_INCLUDES="-I`nc-config --includedir`"
  fi
fi

if eval "test x$NCAUTO = xno"; then
  NC_LDFLAGS="-L$NC_PREFIX/lib"
  NC_LIBS="-lnetcdf"
else
  NC_LDFLAGS=""
  NC_LIBS="`nc-config --libs`"
fi
AC_SUBST(NC_PREFIX)

AC_LANG_PUSH([C])
AC_OPENMP
AC_C_BIGENDIAN
AC_LANG_POP([C])
__AC_FC_NAME_MANGLING

case $ac_cv_fc_mangling in
  "lower case, underscore, no extra underscore")
    name_mangling="FORTRAN_UNDERSCORE_";;
  "lower case, no underscore")
    name_mangling="FORTRAN_SAME";;
  "upper case, no underscore")
    name_mangling="FORTRAN_CAPS_";;
  "lower case, double underscore")
    name_mangling="FORTRAN_GNUF2C";;
  "lower case, no underscore, no extra underscore")
    name_mangling="FORTRAN_SAME";;
  *)
    name_mangling="FORTRAN_MANGLING_UNKNOWN"
    AC_MSG_WARN([UNKNOWN FORTRAN NAME-MANGLING SCHEME]);;
esac

RR_PATH_NETCDF
RR_PATH_NETCDF_F90

test "$osname"=NONE && osname=`uname -s`
OSNAME=$osname
AC_SUBST(OSNAME)

test "$machinename"=NONE && machinename=`uname -m`
fullhostname=`hostname -f`

# DEFINE VARIABLES ACCORDING TO OS AND COMPILER

if echo $osname | grep -i aix >/dev/null 2>&1; then
   SYSDEF="AIX"
fi
if echo $osname | grep -i darwin >/dev/null 2>&1; then
   SYSDEF="DARWIN"
fi
if echo $osname | grep -i unix_system_v >/dev/null 2>&1; then
   SYSDEF="UNIXSYSTEMV"
fi
if echo $osname | grep -i irix >/dev/null 2>&1; then
   SYSDEF="IRIX"
fi
if echo $osname | grep -i irix64 >/dev/null 2>&1; then
   SYSDEF="IRIX64"
fi
if echo $osname | grep -i linux >/dev/null 2>&1; then
   SYSDEF="LINUX"
fi
if echo $osname | grep -i osf1 >/dev/null 2>&1; then
   SYSDEF="OSF1"
fi
if echo $osname | grep -i super >/dev/null 2>&1; then
   SYSDEF="SUPERUX"
fi
if echo $osname | grep -i sun >/dev/null 2>&1; then
   SYSDEF="SUNOS"
fi
if echo $osname | grep -i t3e >/dev/null 2>&1; then
   SYSDEF="T3E"
fi
if echo $osname | grep -i unicos >/dev/null 2>&1; then
   SYSDEF="UNICOS"
fi
if test -z "$SYSDEF"; then
   AC_MSG_WARN([OPERATING SYSTEM UNKNOWN])
   SYSDEF="UNKNOWNOS"
fi

AM_CONDITIONAL(OSNAME_AIX, [test x"$osname" = xAIX])
AM_CONDITIONAL(OSNAME_LINUX, [test x"$osname" = xLinux])

AM_CONDITIONAL(COMPILER_IBM, [test x"$FC" = xxlf2003])
AM_CONDITIONAL(COMPILER_GNU, [test x"$FC" = xgfortran])
AM_CONDITIONAL(COMPILER_G95, [test x"$FC" = xg95])
AM_CONDITIONAL(COMPILER_SUN, [test x"$FC" = xsunf95])
AM_CONDITIONAL(COMPILER_PGI, [test x"$FC" = xpgf90 -o x"$FC" = xpgf95])
AM_CONDITIONAL(COMPILER_INTEL, [test x"$FC" = xifort])

AM_CONDITIONAL(REGCM_PARALLEL, [test x"$PARALLEL" = "x1"])

# Force Big endian I/O for known compilers
if eval "test x"$FC" = xgfortran"; then
  # Following check seems not to work on SuSE, need investigation
#  GFORTRAN_VERSION="`gfortran --version | head -1 | cut -d " " -f 4`"
#  GFORTRAN_MAJOR="`echo $GFORTRAN_VERSION | cut -d "." -f 1`"
#  GFORTRAN_MINOR="`echo $GFORTRAN_VERSION | cut -d "." -f 2`"
#  echo "gfortran is version $GFORTRAN_VERSION"
#  if eval "test "$GFORTRAN_MAJOR" -lt 4"; then
#    AC_MSG_WARN([We need at least gfortran 4.6.x.....])
#    exit 1
#  else
#    if eval "test "$GFORTRAN_MINOR" -lt 6"; then
#      AC_MSG_WARN([We need at least gfortran 4.6.x.....])
#      exit 1
#    fi
#  fi
  FCFLAGS="-fconvert=big-endian -fno-range-check"
  CPRDEF="GNU"
  AC_SUBST(FCFLAGS)
fi
if eval "test x"$FC" = xifort"; then
  FCFLAGS="-convert big_endian -assume byterecl -heap-arrays"
  CPRDEF="INTEL"
  AC_SUBST(FCFLAGS)
  LIBS=`echo $LIBS | sed -e 's/-Wl,--as-needed//g'`
  AC_SUBST(LIBS)
fi
if eval "test x"$FC" = xg95"; then
  FCFLAGS="-fendian=big -fmultiple-save"
  CPRDEF="GNU"
  AC_SUBST(FCFLAGS)
fi
if eval "test x"$FC" = xpgf90"; then
  FCFLAGS="-byteswapio"
  CPRDEF="PGI"
  AC_SUBST(FCFLAGS)
fi
if eval "test x"$FC" = xxlf2003"; then
  AM_CPPFLAGS="-WF,-DIBM $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
  CPRDEF="XLF"
fi
if eval "test x"$FC" = xsun95"; then
  CPRDEF="WORKSHOP"
fi

# Set debug/production flags
if eval "test x$DEBUG = x1"; then
  if eval "test x"$FC" = xg95"; then
    FCFLAGS="-O0 -g -Wall -Wextra -fbounds-check -ftrace=full $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$FC" = xgfortran"; then
    FCFLAGS="-O2 -g -ffpe-trap=zero,invalid,overflow,underflow -fbacktrace \
            -fdump-core -Wall -pedantic -Wconversion -Wsurprising -Waliasing \
            -Wampersand -Warray-bounds -Wcharacter-truncation \
            -Wline-truncation -Wintrinsics-std -Wno-tabs -Wunderflow \
            -Wunused-parameter -Wno-align-commons $FCFLAGS -fbounds-check \
            -fno-range-check"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$FC" = xifort"; then
    FCFLAGS="-O0 -g \
             -check arg_temp_created,bounds,format,output_conversion,uninit \
             -traceback -warn all -fp-model precise -fp-model source -ftz \
             -fimf-arch-consistency=true $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$FC" = xxlf2003"; then
    FCFLAGS="-O0 -qsuppress=1501-245 -bnoquiet -qcheck -qffltrap=nanq,overflow \
             -qmaxmem=-1 -qzerosize $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$FC" = xsunf95"; then
    FCFLAGS="-g -O0 -C -w3 -ansi -xfilebyteorder=big16:%all \
             -D__DATE__="'`date +%D`'" -D__TIME__="'`date +%T`'" $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$FC" = xpgf90 -o x"$FC" = xpgf95"; then
    FCFLAGS="-Mbounds -Mcache_align -Mchkfpstk -Mchkstk -Mdalign -Mdclchk \
             -Mdepchk -Miomutex -Mrecursive -Msave -Ktrap=fp -O0 -g $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
else
  if eval "test x"$FC" = xg95"; then
    FCFLAGS="-O3 $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$FC" = xgfortran"; then
    FCFLAGS="-O3 -mtune=native $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$FC" = xifort"; then
    FCFLAGS="-O3 -fp-model fast=2 -ip $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$FC" = xxlf2003"; then
    FCFLAGS="-qsuppress=1501-245 -bnoquiet -O3 -qstrict \
             -qmaxmem=-1 -qzerosize $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$FC" = xsunf95"; then
    FCFLAGS="-O3 -xfilebyteorder=big16:%all \
             -D__DATE__="'`date +%D`'" -D__TIME__="'`date +%T`'" $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$FC" = xpgf90 -o x"$FC" = xpgf95"; then
    FCFLAGS="-fastsse $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
fi

MCTPATH="${PWD}/external"
if eval "test x"$FC" = xsunf95"; then
  MCTFLAGS="-M$MCTPATH"
  AM_CPPFLAGS="$MCTFLAGS $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
else
  MCTFLAGS="-I$MCTPATH"
  AM_CPPFLAGS="$MCTFLAGS $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi
AM_LDFLAGS="-L$MCTPATH $AM_LDFLAGS"
AC_SUBST(AM_LDFLAGS)
if eval "test x"$PARALLEL" = x1"; then
  ACX_MPI
else
  FORT_SIZE="real4double8"
  MPIFC=$FC
  LIBS="$LIBS -lmpiser"
  AC_SUBST(MPIFC)
  AC_SUBST(FORT_SIZE)
  AM_CPPFLAGS="-DMPI_SERIAL $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

MPEUCFLAGS="-DSYS$SYSDEF -DCPR$CPRDEF -D$name_mangling"
MPEUFCFLAGS="-DSYS$SYSDEF -DCPR$CPRDEF"
if test -n "$SRKDEF"; then
  MPEUCFLAGS="$MPEUCFLAGS -D$SRKDEF"
  MPEUFCFLAGS="$MPEUFCFLAGS -D$SRKDEF"
fi
AC_SUBST(MPEUCFLAGS)
AC_SUBST(MPEUFCFLAGS)

F90FLAGS="$FCFLAGS"
AC_SUBST(F90FLAGS)

AC_SUBST(ESMFINCPTH)
AC_SUBST(ESMFLIBPTH)
if eval "test x$USECPL = x1"; then
  tmp_inc="`cat $ESMF_LIB/esmf.mk | grep ESMF_F90COMPILEPATHS | awk -F= '{print $2}'`"
  tmp_lib="`cat $ESMF_LIB/esmf.mk | grep ESMF_F90LINKPATHS | awk -F= '{print $2}'`"
  tmp_lib="$tmp_lib `cat $ESMF_LIB/esmf.mk | grep ESMF_F90ESMFLINKLIBS | awk -F= '{print $2}'`"
  ESMFINCPTH="$tmp_inc"
  ESMFLIBPTH="$tmp_lib"
fi

AC_SUBST(ROMS_ROOT)
if eval "test x$USEROMS = x1"; then
  ESMFINCPTH="$ESMFINCPTH -I$ROMS_PREFIX"
  ESMFLIBPTH="$ESMFLIBPTH -L$ROMS_PREFIX -lUTIL"
  ROMS_ROOT="$ROMS_PREFIX"
fi

REGCM_ROOT=${PWD}
AC_SUBST(REGCM_ROOT)
AC_PREFIX_DEFAULT(${PWD})

if eval "test x"$FC" = xsunf95"; then
  CPPFLAGS="$AM_CPPFLAGS -M$REGCM_ROOT/Share"
  AC_SUBST(CPPFLAGS)
else
  CPPFLAGS="$AM_CPPFLAGS -I$REGCM_ROOT/Share"
  AC_SUBST(CPPFLAGS)
fi
LIBS="-L$REGCM_ROOT/Share -lrcmlib $LIBS"
SYSLIBS="$LIBS"
AC_SUBST(SYSLIBS)


AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile external/Makefile Share/Makefile PreProc/Makefile \
                 PreProc/Terrain/Makefile PreProc/ICBC/Makefile \
                 PreProc/CH_EMISSIONS/Makefile \
                 PreProc/CH_EMISSIONS/interpolation/Makefile \
                 PreProc/CH_EMISSIONS/interpolation/tools/Makefile \
                 PreProc/CH_EMISSIONS/interpolation/tools/REGCM_grid/Makefile \
                 PreProc/CLM/Makefile PostProc/Makefile \
                 Main/mpplib/Makefile Main/netlib/Makefile \
                 Main/Makefile Main/batslib/Makefile Main/clmlib/Makefile \
                 Main/chemlib/Makefile Main/chemlib/GAS_CBMZ/Makefile \
                 Main/cumlib/Makefile Main/radlib/Makefile \
                 Main/radlib/RRTMG_SW/Makefile Main/pbllib/Makefile \
                 Tools/Programs/CheckSun/Makefile \
                 Tools/Programs/RegCM_read/Makefile \
                 Tools/Programs/SAV_read/Makefile])

AC_OUTPUT
