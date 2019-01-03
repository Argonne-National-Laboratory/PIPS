/*
 * Scaler.C
 *
 *  Created on: 19.12.2017
 *      Authors: Daniel Rehfeldt, Svenja Uslu
 */

#include "Scaler.h"

Scaler::Scaler(Data* prob, bool bitshifting, bool usesides)
: problem(prob), do_bitshifting(bitshifting), with_sides(usesides)
{
}

Scaler::~Scaler()
{
}
