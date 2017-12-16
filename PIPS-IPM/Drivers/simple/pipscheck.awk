#!/usr/bin/awk -f

BEGIN { 
   i = 0
} 
{
    if( $1 == "EOF" )
    {
       print "#runs: "i
       exit 1
    }
    
    if( $1 != "#" )
    {
       i++
       pipscall = "gamspips.sh -DIR=pipstmp"$1" "$2" "$3" "$4
       print pipscall
       system(pipscall)
    }
}