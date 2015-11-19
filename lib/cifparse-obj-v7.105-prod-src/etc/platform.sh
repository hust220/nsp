#!/bin/sh
#
#  File:  platform.sh
#  Original:  6-Aug-97  J. Westbrook
#
#  Copy the platform specific section of the Makefile in
#                  the current directory to:
#
#                      Makefile.platform. 
#
#
sysid="unknown"
release="0.0.0"

#
#  Use uname(1) to figure out what O/S is running.
#

case `uname -s` in
#   Check if it is a Windows 2000 platform
    CYGWIN_NT-5.0)
      sysid="cygwin"
      ;;

#   Check if it is a Windows XP platform
    CYGWIN_NT-5.1)
      sysid="cygwin"
      ;;

#   Check if it is a Mac OS platform
    Darwin)
#     First check if it is clang-xxxx
      clang_ver=`gcc --version | grep -e "clang-5"`
      if [[ ! -z $clang_ver ]]
      then 
        sysid="darwin4"
      else
#       Check if it is GCC version 4.x
        gcc_ver=`gcc --version | grep -e " 4\."`
        if [[ -z $gcc_ver ]]
        then
#         It is not GCC version 4.x. Check if it is GCC version 3.x
          gcc_ver=`gcc --version | grep -e " 3\."`
          if [[ -z $gcc_ver ]]
          then
#           It is not GCC version 3.x. Check if it is GCC version 2.x
            gcc_ver=`gcc --version | grep -e "2\."`
            if [[ -z $gcc_ver ]]
            then
#             It is not GCC version 2.x either. Production can not be compiled.
              sysid="unknown"
            else
#             It is GCC version 2.x
              sysid="darwin2"
            fi
          else
#           It is GCC version 3.x
            sysid="darwin3"
          fi
        else
#         It is GCC version 4.x
          sysid="darwin4"
        fi
      fi 
    ;;

#   Check if it is a Linux platform
    Linux)
#     Check if it is GCC version 4.x
      gcc_ver=`gcc --version | grep -e " 4\."`
      if [[ -z $gcc_ver ]]
      then
#       It is not GCC version 4.x. Check if it is GCC version 3.x
        gcc_ver=`gcc --version | grep -e " 3\."`
        if [[ -z $gcc_ver ]]
        then
#       It is not GCC version 3.x. Check if it is GCC version 2.x
          gcc_ver=`gcc --version | grep -e "2\."`
          if [[ -z $gcc_ver ]]
          then
#           It is not GCC version 2.x either. Production can not be compiled.
            sysid="unknown"
          else
#           It is GCC version 2.x
            sysid="gnu2"
          fi
        else
#         It is GCC version 3.x
          sysid="gnu3"
        fi
      else
#         It is GCC version 4.x
          sysid="gnu4"
      fi
    ;;

    SunOS)
	    sysid="sunos5"
      ;;
    IRIX)
	release=`/usr/bin/uname -r`
	case "$release" in
	   6.*) sysid="sgi6" ;;
	   5.*) sysid="sgi6" ;;
	   4.*) sysid="sgi6" ;;
	   *)   sysid="unknown" ;;
	esac
      ;;
    IRIX64)
	release=`/usr/bin/uname -r`
	case "$release" in
	   6.*) sysid="sgi6" ;;
	   5.*) sysid="sgi6" ;;
	   *)   sysid="unknown" ;;
	esac
      ;;

    HP-UX)
      case `uname -m` in
	  9000/7**)  sysid=unknown
	  ;;
	  *)         sysid=unknown
	  ;;
      esac
      ;;
    OSF1)
      sysid="osf"
      ;;
esac

#
#  If system type is unknown, write warning message and exit with 1.
#

if [ "$sysid" = "unknown" ] ; then
   echo "Warning: this seems to be an unsupported operating system."
   echo " "
   echo "Supported systems are:"
   echo "  SunOS ...... version 4.1.x and 5.2 or higher"
   echo "  Linux ...... any version "
   echo "  SGI IRIX ... version 5.3-6.4"
#  echo "  HP-UX ...... version 9.x, on HP9000/7xx computers"
   exit 1
fi
echo "Using platform specific configuration for system: $sysid"

#
#  Create source filename, test for presence.
#

source=make.platform.$sysid

if [ ! -r "$source" ] ; then
   echo "Can't locate the following file:"
   echo " "
   echo "   $source"
   echo " "
   exit 1
fi

#
#  Copy the file.
#

dest=Makefile.platform
#echo "Copying $source"
#echo "     to $dest"
rm -f $dest
cp $source $dest

# Check if make is GNU make
make_prog=`make --version 2>&1 | grep -e "GNU"`
if [ x"$make_prog" = x ]
then
   echo
   echo "Critical Build Error: GNU Make not accessible!"
   echo "    GNU Make must be accessible in order to build this software."
   echo "    If GNU make is already installed on this platform, make it"
   echo "    accessible to this user. If not installed, please install it."
   echo
   exit 1
fi

# Add source distribution preparation statements to the end of the
# platform makefile.
echo >> $dest
echo "# Added by platform.sh script" >> $dest
echo EXPORT=perl ../../etc/fileUpdate.pl >> $dest
echo EXPORT_LIST=../../etc/export_list >> $dest
echo EXPORT_DIR=export_dir >> $dest

echo "Platform configuration done."
exit 0
