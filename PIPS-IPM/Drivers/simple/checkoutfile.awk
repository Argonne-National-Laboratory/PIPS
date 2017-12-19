#!/usr/bin/awk -f

BEGIN { 
   val = 0
   status = "fail"
} 
{
    { 
        if( $0 ~ /*** SUCCESSFUL TERMINATION ***/ )
           status = "ok"
            
        if( $1 == "Objective:" )    
           val = $2
    }
}
END {
    if( status == "fail" ) 
       print "fail"
    else
       print val
}