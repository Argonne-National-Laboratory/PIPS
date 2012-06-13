
# MAKE-PACKAGE
# Builds pkgIndex.tcl
# This enables the package require statement in Tcl

set name     $env(LEAF_PKG)
set version  0.0
set leaf_so  $env(LEAF_SO)
set leaf_tcl $env(LEAF_TCL)

puts [ ::pkg::create -name $name -version $version \
           -load $leaf_so ]
#  -source $leaf_tcl
