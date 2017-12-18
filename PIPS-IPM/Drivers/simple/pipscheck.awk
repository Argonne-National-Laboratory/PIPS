#!/usr/bin/awk -f

BEGIN { 
   runs = 0
   eps = 0.5
   fails = 0
   file = ARGV[1]
} 
{
    if( $1 == "EOF" )
    {
       print "#runs: "i
       exit 1
    }
    
    if( $1 != "#" )
    {
       runs++
       pipscall = "gamspips.sh -DIR=pipstmp"$1" "$2" "$3" "$4
       #system(pipscall)
       
       # check .out file       
       cmd = "checkoutfile.awk pipstmp"$1"/pips.out"
       cmd | getline newresult
       
       # check .sol file
       cmd = "checksolfile.awk del.sol"
       cmd | getline oldresult
       
       # compare results
       if( newresult == "fail" )
       {
          print $1" failed"
          fails++
       }
       else if( oldresult == "notfound" )
       {
          print $1" no sol"
       }
       else if( newresult - oldresult < eps && oldresult - newresult < eps )
       {
          print $1" ok"
       }
       else
       {
          print $1" failed (wrong obj)"
          fails++
       }
    }
}
END {
    if( fails > 0 )
    {
       print "++++++++++++++++++\ncheck FAILED: nfails: "fails" (of "runs" runs)\n++++++++++++++++++"
    }
    else
    {
       print "++++++++++++++++++\ncheck passed ("runs" runs)\n++++++++++++++++++"
    }
}