# TclXML combo package index file - handcrafted
#
# $Id: pkgIndex.tcl.in,v 1.13 2003/12/03 20:06:34 balls Exp $

namespace eval ::xml {
    variable _init 0
}
namespace eval ::xml::libxml2 {
    variable _init 0
}
namespace eval ::dom {
    variable _init 0
}
namespace eval ::dom::libxml2 {
    variable _init 0
}
namespace eval ::xslt {
    variable _init 0
}

# From http://wiki.tcl.tk/9427
proc ::xml::_platform {} {
    global tcl_platform

    set plat [lindex $tcl_platform(os) 0]
    set mach $tcl_platform(machine)
    switch -glob -- $mach {
	sun4* {
	    set mach sparc
	}
	intel -
	i*86* {
	    set mach x86
	}
	{Power Macintosh} {
	    set mach ppc
	}
    }
    return "$plat-$mach"
}

proc ::xml::_loadlib {dir package version} {
    global tcl_platform

    set lib $package[info sharedlibextension]
    return [file join $dir [_platform] $lib]
}

namespace eval ::xml {
    variable pkginit

    if {![info exists pkginit]} {
	set pkginit 0
    }

    # Try to locate the binary library:
    # 1) Using TEA conventions
    # 2) Using StarKit conventions
    # 3) Using platform-specific conventions

    proc pkgload {dir {binary 0} {cmd {}}} {
	variable pkginit

	if {$pkginit} {return {}}

	namespace eval :: [format {
	    package require xmldefs 3.2
	    package require xml::tcl 3.2
	    # TEA style
	    if {[catch {load [file join [list %s] libtclxml3.2.a] Tclxml}]} {
		# StarKit style
		if {[catch {load [::xml::_loadlib [list %s] Tclxml 3.2]}]} {
		    # Mac OS X frameworks are different
	            if {[catch {load [file join [list %s] .. .. Tclxml] Tclxml}]} {
		        # Unable to load binary implmentation,
		        # just use pure-Tcl implmentation instead
		        if {$binary} {
			    return -code error "unable to load shared library"
		        }
	            } else {
		        set ::xml::libxml2::_init 1
		        set ::dom::libxml2::_init 1
		        set ::xslt::_init 1
		        source [file join [list %s] tcldom-libxml2.tcl]
		        source [file join [list %s] tclxslt-libxslt.tcl]
	            }
	        } else {
		    set ::xml::libxml2::_init 1
		    set ::dom::libxml2::_init 1
		    set ::xslt::_init 1
		    source [file join [list %s] tcldom-libxml2.tcl]
		    source [file join [list %s] tclxslt-libxslt.tcl]
	        }
	    } else {
		set ::xml::libxml2::_init 1
		set ::dom::libxml2::_init 1
		set ::xslt::_init 1
		source [file join [list %s] tcldom-libxml2.tcl]
		source [file join [list %s] tclxslt-libxslt.tcl]
	    }
	    package require xml::tclparser 3.2
	    package provide tclparser 3.2
	    package provide xml::libxml2 3.2
	    package provide xml 3.2
	    package provide dom 3.2
	    package provide dom::libxml2 3.2
	    package provide xslt 3.2
	    package provide xslt::libxslt 3.2

	    set pkginit 1
	} $dir $dir $dir $dir $dir $dir $dir $dir $dir $dir]

	eval $cmd
    }
}

package ifneeded xml::tcl     3.2 [list source [file join $dir xml__tcl.tcl]]
package ifneeded sgmlparser   1.1       [list source [file join $dir sgmlparser.tcl]]
package ifneeded xpath        1.0       [list source [file join $dir xpath.tcl]]
package ifneeded xmldep       1.0       [list source [file join $dir xmldep.tcl]]

# Requesting a specific package means we want it to be the default parser class.

package ifneeded xml::libxml2 3.2 [list ::xml::pkgload $dir 1 {::xml::parser default libxml2}]

# tclparser works with either xml::c or xml::tcl
package ifneeded tclparser 3.2 [list ::xml::pkgload $dir 0 {
    ::xml::parser default tclparser
    package provide tclparser 3.2
}]

# use tcl only (mainly for testing)
package ifneeded puretclparser 3.2 "
    package require xml::tcl       3.2
    package require xmldefs
    package require xml::tclparser 3.2
    package provide puretclparser  3.2
"

# Requesting the generic package leaves the choice of default parser automatic

package ifneeded xml 3.2 [list ::xml::pkgload $dir 0]
package ifneeded dom 3.2 [list ::xml::pkgload $dir 0]
package ifneeded dom::libxml2 3.2 [list ::xml::pkgload $dir 1]
package ifneeded xslt 3.2 [list ::xml::pkgload $dir 1]
package ifneeded xslt::libxslt 3.2 [list ::xml::pkgload $dir 1]

package ifneeded xmlswitch 3.2 [list source [file join $dir xmlswitch.tcl]]

package ifneeded xslt::cache 3.2 [list source [file join $dir xsltcache.tcl]]
package ifneeded xslt::utilities 1.2 [list source [file join $dir utilities.tcl]]
package ifneeded xslt::process 1.1 [list source [file join $dir process.tcl]]
package ifneeded xslt::resources 1.3 [list source [file join $dir resources.tcl]]

if {[info tclversion] <= 8.0} {
    package ifneeded sgml           1.9       [list source [file join $dir sgml-8.0.tcl]]
    package ifneeded xmldefs        3.2 [list source [file join $dir xml-8.0.tcl]]
    package ifneeded xml::tclparser 3.2 [list source [file join $dir tclparser-8.0.tcl]]
} else {
    package ifneeded sgml           1.9       [list source [file join $dir sgml-8.1.tcl]]
    package ifneeded xmldefs        3.2 [list source [file join $dir xml-8.1.tcl]]
    package ifneeded xml::tclparser 3.2 [list source [file join $dir tclparser-8.1.tcl]]
}


