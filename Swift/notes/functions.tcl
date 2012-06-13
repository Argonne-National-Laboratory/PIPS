
# Turbine wrappers

namespace eval functions {

    proc readSolution { stack b dataPath solutionPath scen } {

        rule "readSolution-$b" [ list $dataPath $solutionPath $scen ] \
            $turbine::WORK \
            "readSolution_body $b $dataPath $solutionPath $scen"
    }

    proc readSolution_body { b dataPath solutionPath scen } {

        set dp_value [ turbine::retrieve_string $dataPath ]
        # set sp_value [ turbine::retrieve_string $solutionPath ]
        set s_value [ turbine::retrieve_integer $scen ]

        functions::readSolution
    }

    proc double_list { data } {
        # length in bytes
        set length [ Data_length_get $data ]
        # divide to get count of doubles
        set n [ expr $length / 8 ]
        set result [ list ]
        for { set i 0 } { $i < $n } { incr i } {
            lappend result [ Data_double_get $data $i ]
        }
        return $result
    }
}

