#include "blocksqp.hpp"

namespace blockSQP
{

/**
 * Standard Constructor:
 * Default settings
 */
SQPoptions::SQPoptions()
{
    eps = 2.2204e-16;
    inf = 1.0e20;
    opttol = 1.0e-5;
    nlinfeastol = 1.0e-5;

    // General algorithmic options

    // 0: no globalization, 1: filter line search
    globalization = 1;

    // 0: no feasibility restoration phase 1: if line search fails, start feasibility restoration phase
    restoreFeas = 0;

    // 0: globalization is always active, 1: take a full step at first SQP iteration, no matter what
    skipFirstGlobalization = false;

    // 0: one update for large Hessian, 1: apply updates blockwise, 2: 2 blocks: 1 block updates, 1 block Hessian of obj.
    blockHess = 1;

    // after too many consecutive skipped updates, Hessian block is reset to (scaled) identity
    maxConsecSkippedUpdates = 5;

    // (only for Multiple Shooting OED) second derivative of objective: via standard update (0) or exact (1)
    objSecondDerv = 0;

    // only for simplified Multiple Shooting OED: second derivative of constraints to form exact Hessian
    conSecondDerv = 0;

    // 0: initial Hessian is diagonal matrix, 1: scale initial Hessian according to Nocedal p.143,
    // 2: scale initial Hessian with Oren-Luenberger factor 3: geometric mean of 1 and 2
    // 4: centered Oren-Luenberger sizing according to Tapia paper
    // 5: Byrd-Lu scaling (for L-SR1 only)
    hessScaling = 2;
    fallbackScaling = 4;
    iniHessDiag = 1.0;

    // Damping factor for Powell modification of BFGS updates ( between 0.0 and 1.0 )
    hessDamp = 0.2;

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

    // maximum number of SOC line search iterations
    maxSOCiter = 2;

    // Oren-Luenberger scaling parameters
    colEps = 0.1;
    colTau1 = 0.5;
    colTau2 = 1.0e4;

    /*
     * Magic parameters (mostly from IPOPT paper [Waechter, Biegler 2006])
     */
    // Filter line search parameters
    gammaTheta = 1.0e-5;
    gammaF = 1.0e-5;
    kappaSOC = 0.99;
    kappaF = 0.999;
    thetaMax = 1.0e7; // reject steps if constr viol. is larger than thetaMax
    thetaMin = 1.0e-5; // if constr viol. is smaller than thetaMin require Armijo cond. for obj.
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
    if( conSecondDerv )
    {
        hessUpdate = 4;
        blockHess = 1;
    }

    // If we don't use limited memory BFGS we need to store only one vector.
    if( !hessLimMem )
        hessMemsize = 1;


    // Don't do analytical Hessian for standard OC problems
    //if( vars->objLo < -1.0 )
        //objSecondDerv = 0;
}

} // namespace blockSQP
