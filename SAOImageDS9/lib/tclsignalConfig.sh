# tclsignalConfig.sh --
#
# This shell script (for sh) is generated automatically by tclsignal's
# configure script.  It will create shell variables for most of
# the configuration options discovered by the configure script.
# This script is intended to be included by the configure scripts
# for tclsignal extensions so that they don't have to figure this all
# out for themselves.  This file does not duplicate information
# already provided by tclConfig.sh, so you may need to use that
# file in addition to this one.
#
# The information in this file is specific to a single platform.

# tclsignal's version number.
tclsignal_VERSION='1.4.4.1'

# The name of the tclsignal library (may be either a .a file or a shared library):
tclsignal_LIB_FILE=libtclsignal1.4.4.1.a

# String to pass to linker to pick up the tclsignal library from its
# build directory.
tclsignal_BUILD_LIB_SPEC='-L/usr/local/src/SAOImageDS9/tclsignal -ltclsignal1.4.4.1'

# String to pass to linker to pick up the tclsignal library from its
# installed directory.
tclsignal_LIB_SPEC='-L/usr/local/src/SAOImageDS9/lib/tclsignal1.4.4.1 -ltclsignal1.4.4.1'

# The name of the tclsignal stub library (a .a file):
#tclsignal_STUB_LIB_FILE=libtclsignalstub1.4.4.1.a

# String to pass to linker to pick up the tclsignal stub library from its
# build directory.
#tclsignal_BUILD_STUB_LIB_SPEC='-L/usr/local/src/SAOImageDS9/tclsignal -ltclsignalstub1.4.4.1'

# String to pass to linker to pick up the tclsignal stub library from its
# installed directory.
#tclsignal_STUB_LIB_SPEC='-L/usr/local/src/SAOImageDS9/lib/tclsignal1.4.4.1 -ltclsignalstub1.4.4.1'

# String to pass to linker to pick up the tclsignal stub library from its
# build directory.
#tclsignal_BUILD_STUB_LIB_PATH='/usr/local/src/SAOImageDS9/tclsignal/libtclsignalstub1.4.4.1.a'

# String to pass to linker to pick up the tclsignal stub library from its
# installed directory.
#tclsignal_STUB_LIB_PATH='/usr/local/src/SAOImageDS9/lib/tclsignal1.4.4.1/libtclsignalstub1.4.4.1.a'
