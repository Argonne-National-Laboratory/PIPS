
# Turbine wrappers

namespace eval functions {

    proc readSolution_turbine { stack output inputs } {

        set dataPath     [ lindex $inputs 0 ]
        set solutionPath [ lindex $inputs 1 ]
        set scen         [ lindex $inputs 2 ]

        turbine::rule "readSolution-$output" $inputs \
            $turbine::WORK \
            "functions::readSolution_body $output $dataPath $solutionPath $scen"
    }

    proc readSolution_body { b dataPath solutionPath scen } {

        puts readSolution_body

        set dp_value [ turbine::retrieve_string $dataPath ]
        # set sp_value [ turbine::retrieve_string $solutionPath ]
        set s_value [ turbine::retrieve_integer $scen ]

        set d [ readSolution X X 3 ]
        puts OKAY
    }
}
