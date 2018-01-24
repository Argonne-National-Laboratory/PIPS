#!/usr/bin/awk -f

BEGIN { 
   runs = 0
   eps = 0.5
   fails = 0
   split(ARGV[1], array, ".")
   testfile = array[1]
   
   resfile = testfile".res"
   
   print "TIME       ITERATIONS      STATE             OBJ   \n" > resfile
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
       pipscall = "sh gamspipstime.sh -DIR="instancedir" "$2" "$3" "$4" "$5
       system(pipscall)
       
       # check .out file       
       cmd = "awk -f checkoutfile2.awk "instancedir"/pips.out"
       cmd | getline result
       
       # check result
       print result"\n" > testfile".res"
    }
}
END {
    system("less "resfile)
}