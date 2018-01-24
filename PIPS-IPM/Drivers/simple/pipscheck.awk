#!/usr/bin/awk -f

BEGIN { 
   runs = 0
   eps = 0.5
   fails = 0
   split(ARGV[1], array, ".")
   testfile = array[1]
} 
{
    if( $1 == "EOF" )
    {
       exit 0
    }
    
    if( $1 != "#" )
    {
       instancename = $1
       instancedir = "pipstmp"instancename 
       runs++
       pipscall = "sh gamspips.sh -DIR="instancedir" "$2" "$3" "$4" "$5
       system(pipscall)
       
       # check .out file       
       cmd = "awk -f checkoutfile.awk "instancedir"/pips.out"
       cmd | getline newresult
       
       # check .sol file
       cmd = "awk -f checksolfile.awk -v instance="instancename" "testfile".sol"
       cmd | getline oldresult
       
       # compare results
       if( newresult == "fail" )
       {
          print $1" failed"
          fails++
          failsarr[fails] = instancename" not solved to optimality"
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
          failsarr[fails] = instancename" has wrong objective"
       }
    }
}
END {
    if( fails > 0 )
    {
       print "\n++++++++++++++++++\ncheck FAILED: nfails: "fails" (of "runs" runs)\n++++++++++++++++++"
       for( i = 1; i <= fails; i++ )
           print failsarr[i] "\n"
    }
    else
    {
       print "\n++++++++++++++++++\ncheck passed ("runs" runs)\n++++++++++++++++++"
    }
}
