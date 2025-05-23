# tksvgConfig.sh --
#
# This shell script (for sh) is generated automatically by tksvg's
# configure script.  It will create shell variables for most of
# the configuration options discovered by the configure script.
# This script is intended to be included by the configure scripts
# for tksvg extensions so that they don't have to figure this all
# out for themselves.  This file does not duplicate information
# already provided by tclConfig.sh, so you may need to use that
# file in addition to this one.
#
# The information in this file is specific to a single platform.

# tksvg's version number.
tksvg_VERSION='0.7'

# The name of the tksvg library (may be either a .a file or a shared library):
tksvg_LIB_FILE=libtksvg0.7.a

# String to pass to linker to pick up the tksvg library from its
# build directory.
tksvg_BUILD_LIB_SPEC='-L/usr/local/src/SAOImageDS9/tksvg -ltksvg0.7'

# String to pass to linker to pick up the tksvg library from its
# installed directory.
tksvg_LIB_SPEC='-L/usr/local/src/SAOImageDS9/lib/tksvg0.7 -ltksvg0.7'

# The name of the tksvg stub library (a .a file):
#tksvg_STUB_LIB_FILE=libtksvgstub0.7.a

# String to pass to linker to pick up the tksvg stub library from its
# build directory.
#tksvg_BUILD_STUB_LIB_SPEC='-L/usr/local/src/SAOImageDS9/tksvg -ltksvgstub0.7'

# String to pass to linker to pick up the tksvg stub library from its
# installed directory.
#tksvg_STUB_LIB_SPEC='-L/usr/local/src/SAOImageDS9/lib/tksvg0.7 -ltksvgstub0.7'

# String to pass to linker to pick up the tksvg stub library from its
# build directory.
#tksvg_BUILD_STUB_LIB_PATH='/usr/local/src/SAOImageDS9/tksvg/libtksvgstub0.7.a'

# String to pass to linker to pick up the tksvg stub library from its
# installed directory.
#tksvg_STUB_LIB_PATH='/usr/local/src/SAOImageDS9/lib/tksvg0.7/libtksvgstub0.7.a'
