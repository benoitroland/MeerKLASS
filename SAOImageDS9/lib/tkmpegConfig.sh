# tkmpegConfig.sh --
#
# This shell script (for sh) is generated automatically by tkmpeg's
# configure script.  It will create shell variables for most of
# the configuration options discovered by the configure script.
# This script is intended to be included by the configure scripts
# for tkmpeg extensions so that they don't have to figure this all
# out for themselves.  This file does not duplicate information
# already provided by tclConfig.sh, so you may need to use that
# file in addition to this one.
#
# The information in this file is specific to a single platform.

# tkmpeg's version number.
tkmpeg_VERSION='1.0'

# The name of the tkmpeg library (may be either a .a file or a shared library):
tkmpeg_LIB_FILE=libtkmpeg1.0.a

# String to pass to linker to pick up the tkmpeg library from its
# build directory.
tkmpeg_BUILD_LIB_SPEC='-L/usr/local/src/SAOImageDS9/tkmpeg -ltkmpeg1.0'

# String to pass to linker to pick up the tkmpeg library from its
# installed directory.
tkmpeg_LIB_SPEC='-L/usr/local/src/SAOImageDS9/lib/tkmpeg1.0 -ltkmpeg1.0'

# The name of the tkmpeg stub library (a .a file):
#tkmpeg_STUB_LIB_FILE=libtkmpegstub1.0.a

# String to pass to linker to pick up the tkmpeg stub library from its
# build directory.
#tkmpeg_BUILD_STUB_LIB_SPEC='-L/usr/local/src/SAOImageDS9/tkmpeg -ltkmpegstub1.0'

# String to pass to linker to pick up the tkmpeg stub library from its
# installed directory.
#tkmpeg_STUB_LIB_SPEC='-L/usr/local/src/SAOImageDS9/lib/tkmpeg1.0 -ltkmpegstub1.0'

# String to pass to linker to pick up the tkmpeg stub library from its
# build directory.
#tkmpeg_BUILD_STUB_LIB_PATH='/usr/local/src/SAOImageDS9/tkmpeg/libtkmpegstub1.0.a'

# String to pass to linker to pick up the tkmpeg stub library from its
# installed directory.
#tkmpeg_STUB_LIB_PATH='/usr/local/src/SAOImageDS9/lib/tkmpeg1.0/libtkmpegstub1.0.a'
