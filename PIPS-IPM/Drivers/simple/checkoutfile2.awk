#!/usr/bin/awk -f

BEGIN { 
   val = 0
   status = "fail"
   time = -1.0
   iterations = -1;
} 
{
    { 
        if( $0 ~ /*** SUCCESSFUL TERMINATION ***/ )
           status = "ok"
            
        if( $1 == "Objective:" )    
           val = $2
           
        if( $0 ~ "elapsed" )
           time = substr($3,1,7) 
           
        if( $2 == "Iteration" )
           iterations = $3 
    }
}
END {
       print time"            "iterations"         "status"         "val
}