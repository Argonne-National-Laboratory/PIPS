
# ROUNDING.TCL
#
# Swift/T dataflow wrappers

package require data 0.0

namespace eval rounding {

    proc evaluateRecourseLP_turbine { stack output inputs } {

        # Set up inputs
        set dataPath          [ lindex $inputs 0 ]
        set nScen             [ lindex $inputs 1 ]
        set scen              [ lindex $inputs 2 ]
        set candidateSolution [ lindex $inputs 3 ]

        # Issue data dependent function call
        turbine::rule "evaluateRecourseLP-$output" $inputs \
            $turbine::WORK \
            "rounding::evaluateRecourseLP_body $output $dataPath $nScen $scen $candidateSolution"
    }

    proc evaluateRecourseLP_body { result dataPath nScen scen candidateSolution } {

        # Retrieve inputs
        set dp_value  [ turbine::retrieve_string $dataPath ]
        set ns_value  [ turbine::retrieve_integer $nScen ]
        set s_value   [ turbine::retrieve_integer $scen ]
        set cs        [ adlb::retrieve_blob $candidateSolution ]
        set pointerCS [ lindex $cs 0 ]
        set pointerCS [ Data_cast_to_pointer $pointerCS ]
        set lengthCS  [ lindex $cs 1 ]

        # Call C++ function
        set r_value [ evaluateRecourseLP $dp_value $ns_value $s_value $pointerCS $lengthCS ]

        # Store output
        turbine::store_float $result $r_value
    }

    proc readConvSolution_turbine { stack output inputs } {

        # Set up inputs
        set dataPath [ lindex $inputs 0 ]
        set convPath [ lindex $inputs 1 ]

        # Issue data dependent function call
        turbine::rule "readConvSolution-$output" $inputs \
            $turbine::WORK \
            "rounding::readConvSolution_body $output $dataPath $convPath"
    }

    proc readConvSolution_body { b dataPath convPath } {

        # Retrieve inputs
        set dp_value [ turbine::retrieve_string $dataPath ]
        set cp_value [ turbine::retrieve_string $convPath ]

        # Call C++ function
        set d [ readConvSolution $dp_value $cp_value ]

        # Pack and store outputs
        set pointer [ data::pointer $d ]
        set length  [ data::length  $d ]
        turbine::store_blob $b $pointer $length
        Data_free $d
    }

    proc roundSolution_turbine { stack output inputs } {

        # Set up inputs
        set convSolution [ lindex $inputs 0 ]
        set cutoff       [ lindex $inputs 1 ]

        # Issue data dependent function call
        turbine::rule  "round-$output" $inputs \
            $turbine::WORK \
            "rounding::roundSolution_body $output $convSolution $cutoff"
    }

    proc roundSolution_body { b convSolution cutoff } {

        # Retrieve inputs
        set cs        [ adlb::retrieve_blob $convSolution ]
        set pointerCS [ lindex $cs 0 ]
        set pointerCS [ Data_cast_to_pointer $pointerCS ]
        set lengthCS  [ lindex $cs 1 ]
        set c_value [ turbine::retrieve_float $cutoff ]

        # Call C++ function
        set d [ roundSolution $pointerCS $lengthCS $c_value ]

        # Pack and store outputs
        set pointer [ data::pointer $d ]
        set length  [ data::length  $d ]
        turbine::store_blob $b $pointer $length
        Data_free $d
    }
}
