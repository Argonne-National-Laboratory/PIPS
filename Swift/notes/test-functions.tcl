
package require functions 0.0

# readSolution "problem.data" "solution.data" 3

set d [ Data_make_test ]
set length [ Data_length_get $d ]
puts $length
Data_free $d

set d [ readSolution X X 3 ]


set t [ Data_double_get $d 0 ]
puts $t

set L [ functions::double_list $d ]
puts $L
