#ifndef BLOCKSQP_DEFS_H
#define BLOCKSQP_DEFS_H

// SYSTEM HEADER FILES
#include "string.h"
#include "math.h"
#include "stdio.h"
#include <set>
#include <limits>

namespace blockSQP
{

static double const myInf = std::numeric_limits<double>::infinity();    ///< Used to mark sparse zeros in Jacobian

typedef char NAMESTR[256];
typedef char PATHSTR[4096];

} // namespace blockSQP

#endif
