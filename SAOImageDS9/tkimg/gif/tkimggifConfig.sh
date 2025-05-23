# tkimggifConfig.sh --
#
# This shell script (for sh) is generated automatically by tkimggif's
# configure script.  It will create shell variables for most of
# the configuration options discovered by the configure script.
# This script is intended to be included by the configure scripts
# for tkimggif extensions so that they don't have to figure this all
# out for themselves.  This file does not duplicate information
# already provided by tclConfig.sh, so you may need to use that
# file in addition to this one.
#
# The information in this file is specific to a single platform.

# tkimggif's version number.
tkimggif_VERSION='1.4.13'

# The name of the tkimggif library (may be either a .a file or a shared library):
tkimggif_LIB_FILE=libtkimggif1.4.13.a

# String to pass to linker to pick up the tkimggif library from its
# build directory.
tkimggif_BUILD_LIB_SPEC='-L/usr/local/src/SAOImageDS9/tkimg/gif -ltkimggif1.4.13'

# String to pass to linker to pick up the tkimggif library from its
# installed directory.
tkimggif_LIB_SPEC='-L/usr/local/src/SAOImageDS9/lib/tkimggif1.4.13 -ltkimggif1.4.13'
