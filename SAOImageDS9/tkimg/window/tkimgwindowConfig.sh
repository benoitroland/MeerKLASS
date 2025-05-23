# tkimgwindowConfig.sh --
#
# This shell script (for sh) is generated automatically by tkimgwindow's
# configure script.  It will create shell variables for most of
# the configuration options discovered by the configure script.
# This script is intended to be included by the configure scripts
# for tkimgwindow extensions so that they don't have to figure this all
# out for themselves.  This file does not duplicate information
# already provided by tclConfig.sh, so you may need to use that
# file in addition to this one.
#
# The information in this file is specific to a single platform.

# tkimgwindow's version number.
tkimgwindow_VERSION='1.4.13'

# The name of the tkimgwindow library (may be either a .a file or a shared library):
tkimgwindow_LIB_FILE=libtkimgwindow1.4.13.a

# String to pass to linker to pick up the tkimgwindow library from its
# build directory.
tkimgwindow_BUILD_LIB_SPEC='-L/usr/local/src/SAOImageDS9/tkimg/window -ltkimgwindow1.4.13'

# String to pass to linker to pick up the tkimgwindow library from its
# installed directory.
tkimgwindow_LIB_SPEC='-L/usr/local/src/SAOImageDS9/lib/tkimgwindow1.4.13 -ltkimgwindow1.4.13'
