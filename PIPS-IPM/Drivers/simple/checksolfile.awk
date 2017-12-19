#!/usr/bin/awk -f

BEGIN { 
   status = "fail"
} 
{
    { 
        if( $1 == instance )
        {
           status = "ok"
           val = $2
        }
    }
}
END {
    if( status == "fail" ) 
       print "fail"
    else
       print val
}