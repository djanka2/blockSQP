#ifndef BLOCKSQP_DEFS_H
#define BLOCKSQP_DEFS_H

// SYSTEM HEADER FILES
#include "math.h"
#include "stdio.h"
#include <set>

//#define MYDEBUGLEVEL 1  ///< Print debug information (impairs performance)
//#define PARALLELQP ///< Always solve two QPs (in parallel)

// Choice of QP solver
//#define QPSOLVER_QPOASES_DENSE
//#define QPSOLVER_QPOASES_SPARSE
#define QPSOLVER_QPOASES_SCHUR

#define NEWQPLOOP

/*--------------------------------------------------------------------*/

#include "qpOASES.hpp"

#if defined QPSOLVER_QPOASES_SPARSE || defined QPSOLVER_QPOASES_SCHUR
#define QPSOLVER_SPARSE
#else
#define QPSOLVER_DENSE
#endif


namespace blockSQP
{

typedef char PATHSTR[4096];

} // namespace blockSQP

#endif
