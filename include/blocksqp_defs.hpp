#ifndef BLOCKSQP_DEFS_H
#define BLOCKSQP_DEFS_H

// SYSTEM HEADER FILES
#include "math.h"
#include "stdio.h"
#include <set>
#include <limits>
#include <string>

//#define MYDEBUG   ///< Print debug information (impairs performance)

// Choice of QP solver
//#define QPSOLVER_QPOPT
//#define QPSOLVER_QPOASES_DENSE
//#define QPSOLVER_QPOASES_SPARSE
#define QPSOLVER_QPOASES_SCHUR

/*--------------------------------------------------------------------*/

#ifdef QPSOLVER_QPOPT
#include "qpopt.h"
#elif defined QPSOLVER_QPOASES_SPARSE || defined QPSOLVER_QPOASES_DENSE || defined QPSOLVER_QPOASES_SCHUR
#define QPSOLVER_QPOASES
#include "qpOASES.hpp"
#endif

#if defined QPSOLVER_QPOASES_SPARSE || defined QPSOLVER_QPOASES_SCHUR
#define QPSOLVER_SPARSE
#else
#define QPSOLVER_DENSE
#endif


namespace blockSQP
{

static double const myInf = std::numeric_limits<double>::infinity();    ///< Used to mark sparse zeros in Jacobian

typedef char NAMESTR[256];
typedef char PATHSTR[4096];

} // namespace blockSQP

#endif
