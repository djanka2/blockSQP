#ifndef DEFS_H
#define DEFS_H

// SYSTEM HEADER FILES
#include "string.h"
#include "math.h"
#include "stdio.h"
//#include "stdlib.h"
#include <set>
#include <limits>
static double const myInf = std::numeric_limits<double>::infinity();    ///< Used to mark sparse zeros in Jacobian


#define MAXNAME 60
typedef char NAMESTR[MAXNAME];

#define MAXPATH 4096
typedef char PATHSTR[MAXPATH];


#endif
