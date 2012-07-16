
# Turbine wrappers

package require data 0.0

namespace eval rounding_functions {

    proc evaluateRecourseLP_turbine { stack output inputs } {

        set dataPath          [ lindex $inputs 0 ]
        set nScen             [ lindex $inputs 1 ]
        set scen              [ lindex $inputs 2 ]
        set candidateSolution [ lindex $inputs 3 ]

        turbine::rule "evaluateRecourseLP-$output" $inputs \
            $turbine::WORK \
            "rounding_functions::evaluateRecourseLP_body $output $dataPath $nScen $scen $candidateSolution"
    }

    proc evaluateRecourseLP_body { b dataPath nScen scen candidateSolution } {

        puts evaluateRecourseLP_body

        set dp_value [ turbine::retrieve_string $dataPath ]
        set ns_value [ turbine::retrieve_integer $nScen ]
        set s_value  [ turbine::retrieve_integer $scen ]

        set d [ evaluateRecourseLP $dp_value $ns_value $s_value 0 ]
        set pointer [ data::pointer $d ]
        set length  [ data::length  $d ]
        adlb::store_blob $b $pointer $length
        Data_free $d
    }
}
