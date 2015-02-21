#include "blocksqp.hpp"

namespace blockSQP
{

/**
 * Standard Constructor:
 * Default settings
 */
SQPoptions::SQPoptions()
{
    // 0: no output, 1: normal output, 2: verbose output
    printLevel = 2;
    // 1: (some) colorful output
    printColor = 1;

    //eps = 2.2204e-16;
    eps = 1.0e-15;
    inf = 1.0e20;
    opttol = 1.0e-5;
    nlinfeastol = 1.0e-5;

    // General algorithmic options

    // 0: no globalization, 1: filter line search
    globalization = 1;

    // 0: no feasibility restoration phase 1: if line search fails, start feasibility restoration phase
    restoreFeas = 1;

    // 0: globalization is always active, 1: take a full step at first SQP iteration, no matter what
    skipFirstGlobalization = false;

    // 0: one update for large Hessian, 1: apply updates blockwise, 2: 2 blocks: 1 block updates, 1 block Hessian of obj.
    blockHess = 1;

    // after too many consecutive skipped updates, Hessian block is reset to (scaled) identity
    maxConsecSkippedUpdates = 100;

    // for which blocks should second derivatives be provided by the user:
    // 0: none, 1: for the last block, 2: for all blocks
    whichSecondDerv = 1;

    // 0: initial Hessian is diagonal matrix, 1: scale initial Hessian according to Nocedal p.143,
    // 2: scale initial Hessian with Oren-Luenberger factor 3: geometric mean of 1 and 2
    // 4: centered Oren-Luenberger sizing according to Tapia paper
    // 5: Byrd-Lu scaling (for L-SR1 only)
    hessScaling = 2;
    fallbackScaling = 4;
    iniHessDiag = 1.0;

    // Activate damping strategy for BFGS (if deactivated, BFGS might yield indefinite updates!)
    hessDamp = 1;

    // Damping factor for Powell modification of BFGS updates ( between 0.0 and 1.0 )
    hessDampFac = 0.2;

    // 0: constant, 1: SR1, 2: BFGS (damped), 3: [not used] , 4: finiteDiff, 5: Gauss-Newton
    hessUpdate = 1;
    fallbackUpdate = 2;

    // 0: full memory updates 1: limited memory
    hessLimMem = 1;

    // memory size for L-BFGS/L-SR1 updates
    hessMemsize = 20;

    // maximum number of line search iterations
    maxLineSearch = 30;

    // if step has to be reduced in too many consecutive iterations, feasibility restoration phase is invoked
    maxConsecReducedSteps = 100;

    // maximum number of second-order correction steps
    maxSOCiter = 3;

    // Oren-Luenberger scaling parameters
    colEps = 0.1;
    colTau1 = 0.5;
    colTau2 = 1.0e4;

    // Filter line search parameters
    gammaTheta = 1.0e-5;
    gammaF = 1.0e-5;
    kappaSOC = 0.99;
    kappaF = 0.999;
    thetaMax = 1.0e7;       // reject steps if constr viol. is larger than thetaMax
    thetaMin = 1.0e-5;      // if constr viol. is smaller than thetaMin require Armijo cond. for obj.
    delta = 1.0;
    sTheta = 1.1;
    sF = 2.3;
    eta = 1.0e-4;

    // Inertia correction for filter line search and indefinite Hessians
    kappaMinus = 0.333;
    kappaPlus = 8.0;
    kappaPlusMax = 100.0;
    deltaH0 = 1.0e-4;
}


/**
 * Some options cannot be set together, resolve here
 */
void SQPoptions::optionsConsistency()
{
    // If we compute second constraints derivatives switch to finite differences Hessian (convenience)
    if( whichSecondDerv == 2 )
    {
        hessUpdate = 4;
        blockHess = 1;
    }

    // If we don't use limited memory BFGS we need to store only one vector.
    if( !hessLimMem )
        hessMemsize = 1;

#ifndef QPSOLVER_QPOASES_SCHUR
    if( hessUpdate == 1 )
    {
        printf( "SR1 update only works with qpOASES Schur complement version. Using BFGS updates instead.\n" );
        hessUpdate = 2;
        hessScaling = fallbackScaling;
    }
#endif

    // Don't do analytical Hessian for standard OC problems
    //if( vars->objLo < -1.0 )
        //objSecondDerv = 0;
}

} // namespace blockSQP
