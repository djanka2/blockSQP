#include "blocksqp.hpp"
#include "blocksqp_general_purpose.hpp"

namespace blockSQP
{

#ifdef NEWQPLOOP
/**
 * Inner loop of SQP algorithm:
 * Solve a sequence of QPs until pos. def. assumption is satisfied.
 */
int SQPmethod::solveQP( Matrix &deltaXi, Matrix &lambdaQP, int flag )
{
    int l, maxQP = 1;
    if( param->globalization == 1 && flag == 0 && stats->itCount > 1 )
        if( param->hessUpdate == 1 || param->hessUpdate == 4 )
            maxQP = 2;

    /*
     * Prepare for qpOASES
     */

    // Setup QProblem data
    qpOASES::Matrix *A;
    qpOASES::SymSparseMat *H;
    if( flag == 0 )
        A = new qpOASES::SparseMatrix( prob->nCon, prob->nVar,
                             vars->jacIndRow, vars->jacIndCol, vars->jacNz );
    double *g = vars->gradObj.ARRAY();
    double *lb = vars->deltaBl.ARRAY();
    double *lu = vars->deltaBu.ARRAY();
    double *lbA = vars->deltaBl.ARRAY() + prob->nVar;
    double *luA = vars->deltaBu.ARRAY() + prob->nVar;

    // Save active set from the last successful iteration
    qpOASES::Bounds B;
    qp->getBounds( B );
    qpOASES::Constraints C;
    qp->getConstraints( C );

    // qpOASES options
    qpOASES::Options opts;
    if( flag == 0 && maxQP == 2 )
        opts.enableInertiaCorrection = qpOASES::BT_FALSE;
    else
        opts.enableInertiaCorrection = qpOASES::BT_TRUE;
    opts.enableEqualities = qpOASES::BT_TRUE;
    opts.initialStatusBounds = qpOASES::ST_INACTIVE;
    opts.printLevel = qpOASES::PL_NONE;
    opts.numRefinementSteps = 2;
    opts.epsLITests =  2.2204e-08;
    qp->setOptions( opts );

    /// \todo only need active set
    // Store last successful QP in temporary storage
    qpSave = *qp;

    // Other variables for qpOASES
    double cpuTime = 10000.0;
    int maxIt = (flag==0) ? 2500 : 100;
    qpOASES::SolutionAnalysis solAna;
    qpOASES::returnValue ret;

    /*
     * QP solving loop for convex combinations (sequential)
     */
    for( l=0; l<maxQP; l++ )
    {
        /*
         * Set Hessian
         */
        if( l > 0 )
        {// If the solution of the first QP was rejected, consider second Hessian
            stats->qpResolve++;
            *qp = qpSave;

            // Compute fallback update only once
            if( l == 1 )
            {
                vars->hess = vars->hess2;

                // Limited memory: compute fallback update only when needed
                if( param->hessLimMem )
                {
                    // If last block contains exact Hessian, we need to copy it
                    if( param->whichSecondDerv == 1 )
                        for( int i=0; i<vars->hess[prob->nBlocks-1].M(); i++ )
                            for( int j=i; j<vars->hess[prob->nBlocks-1].N(); j++ )
                                vars->hess2[prob->nBlocks-1]( i,j ) = vars->hess1[prob->nBlocks-1]( i,j );

                    stats->itCount--;
                    int hessDampSave = param->hessDamp;
                    param->hessDamp = 1;
                    calcHessianUpdateLimitedMemory( param->fallbackUpdate, param->fallbackScaling );
                    param->hessDamp = hessDampSave;
                    stats->itCount++;


                }
                // Full memory: both updates must be computed in every iteration
            }

            // 'Nontrivial' convex combinations
            if( maxQP > 2 )
            {
                /* Convexification parameter: mu_l = l / (maxQP-1).
                 * Compute it only in the first iteration, afterwards update
                 * by recursion: mu_l/mu_(l-1) */
                double mu = (l==1) ? 1.0 / (maxQP-1) : ((double) l)/((double) (l-1));
                double mu1 = 1.0 - mu;
                for( int iBlock=0; iBlock<vars->nBlocks; iBlock++ )
                    for( int i=0; i<vars->hess[iBlock].M(); i++ )
                        for( int j=i; j<vars->hess[iBlock].N(); j++ )
                        {
                            vars->hess2[iBlock]( i,j ) *= mu;
                            vars->hess2[iBlock]( i,j ) += mu1 * vars->hess1[iBlock]( i,j );
                        }
            }
        }

        /*
         * Call qpOASES
         */
        if( flag == 0 )
        {
            // Set sparse Hessian
            vars->convertHessian( prob, param->eps, vars->hess, vars->hessNz,
                              vars->hessIndRow, vars->hessIndCol, vars->hessIndLo );
            H = new qpOASES::SymSparseMat( prob->nVar, prob->nVar,
                                       vars->hessIndRow, vars->hessIndCol, vars->hessNz );
            H->createDiagInfo();

            // Call qpOASES
            if( qp->getStatus() == qpOASES::QPS_HOMOTOPYQPSOLVED ||
                qp->getStatus() == qpOASES::QPS_SOLVED )
            {
                ret = qp->hotstart( H, g, A, lb, lu, lbA, luA, maxIt, &cpuTime );
            }
            else
                ret = qp->init( H, g, A, lb, lu, lbA, luA, maxIt, &cpuTime );
        }
        else if( flag == 1 ) // Second order correction: H and A do not change
            ret = qp->hotstart( g, lb, lu, lbA, luA, maxIt, &cpuTime );

        /*
         * Check assumption (G3*) if nonconvex QP was solved
         */
        if( l < maxQP-1 && flag == 0 )
        {
            // QP was solved successfully and positive curvature after removing bounds
            if( ret == qpOASES::SUCCESSFUL_RETURN &&
                solAna.checkCurvatureOnStronglyActiveConstraints( qp ) == qpOASES::SUCCESSFUL_RETURN )
            {
                stats->qpIterations = maxIt + 1;
                break;
            }
            else
            {
                // Statistics: save iterations from unsuccessful QP
                if( ret == qpOASES::RET_SETUP_AUXILIARYQP_FAILED )
                    stats->qpIterations2++;
                else
                    stats->qpIterations2 += maxIt + 1;
                // Count total number of rejected SR1 updates
                stats->rejectedSR1++;
            }
        }
        else // Convex QP was solved
            stats->qpIterations += maxIt + 1;

    } // End of QP solving loop

    /*
     * Post-processing
     */

    // Read solution
    qp->getPrimalSolution( deltaXi.ARRAY() );
    qp->getDualSolution( lambdaQP.ARRAY() );
    vars->qpObj = qp->getObjVal();

    // Compute constrJac*deltaXi, need this for second order correction step
    #ifdef QPSOLVER_SPARSE
    sparseAtimesb( vars->jacNz, vars->jacIndRow, vars->jacIndCol, deltaXi, vars->AdeltaXi );
    #else
    Atimesb( vars->constrJac, deltaXi, vars->AdeltaXi );
    #endif

    // Print qpOASES error code, if any
    if( ret != qpOASES::SUCCESSFUL_RETURN && flag == 0 )
        printf( "qpOASES error message: \"%s\"\n", qpOASES::getGlobalMessageHandler()->getErrorCodeMessage( ret ) );

    // Point Hessian again to the desired Hessian
    vars->hess = vars->hess1;

    // For full Hessian: Restore fallback Hessian if convex combinations
    // were used during the loop
    if( !param->hessLimMem && maxQP > 2 && flag == 0 )
    {
        double mu = 1.0 / ((double) (l));
        double mu1 = 1.0 - mu;
        for( int iBlock=0; iBlock<vars->nBlocks; iBlock++ )
            for( int i=0; i<vars->hess[iBlock].M(); i++ )
                for( int j=i; j<vars->hess[iBlock].N(); j++ )
                {
                    vars->hess2[iBlock]( i,j ) *= mu;
                    vars->hess2[iBlock]( i,j ) += mu1 * vars->hess1[iBlock]( i,j );
                }
    }

    // Return code depending on qpOASES returnvalue
    /* 0: Success
     * 1: Maximum number of iterations reached
     * 2: Unbounded
     * 3: Infeasible
     * 4: Other error */
    if( ret == qpOASES::SUCCESSFUL_RETURN )
        return 0;
    else if( ret == qpOASES::RET_MAX_NWSR_REACHED )
        return 1;
    else if( ret == qpOASES::RET_HESSIAN_NOT_SPD ||
             ret == qpOASES::RET_HESSIAN_INDEFINITE ||
             ret == qpOASES::RET_INIT_FAILED_UNBOUNDEDNESS ||
             ret == qpOASES::RET_QP_UNBOUNDED ||
             ret == qpOASES::RET_HOTSTART_STOPPED_UNBOUNDEDNESS )
        return 2;
    else if( ret == qpOASES::RET_INIT_FAILED_INFEASIBILITY ||
             ret == qpOASES::RET_QP_INFEASIBLE ||
             ret == qpOASES::RET_HOTSTART_STOPPED_INFEASIBILITY )
        return 3;
    else
        return 4;
}
#endif

/**
 * Always solve 2 QPs: one with the original, one with the fallback method
 */
int SQPmethod::solveQP2( Matrix &deltaXi, Matrix &lambdaQP, int flag )
{
    #ifndef QPSOLVER_SPARSE
    printf( "solveQP2 not implemented for dense qpOASES. Abort.\n" );
    return 4;
    #endif

    int i, maxQP = 1;
    if( param->globalization == 1 && flag == 0 && stats->itCount > 1 )
        if( param->hessUpdate == 1 || param->hessUpdate == 4 )
            maxQP = 2;

    qpOASES::Matrix *A;
    qpOASES::SolutionAnalysis solAna;
    qpOASES::SymmetricMatrix *H1, *H2;
    qpOASES::Options opts, opts2;
    qpOASES::returnValue ret, ret2;
    double *lb, *lu, *lbA, *luA;
    double cpuTime, cpuTime2;
    bool enableQPLoop = false;
    int maxIt, maxIt2;

    // set options for qpOASES
    if( flag == 0 && maxQP == 2 )
        opts.enableInertiaCorrection = qpOASES::BT_FALSE;
    else
        opts.enableInertiaCorrection = qpOASES::BT_TRUE;
    opts.enableEqualities = qpOASES::BT_TRUE;
    opts.initialStatusBounds = qpOASES::ST_INACTIVE;
    opts.printLevel = qpOASES::PL_NONE;
    opts.numRefinementSteps = 2;
    opts.epsLITests =  2.2204e-08;
    opts2 = opts;
    opts2.enableInertiaCorrection = qpOASES::BT_TRUE;
    cpuTime = cpuTime2 = 10000.0;

    // Create Jacobian, is the same for both QPs
    if( flag == 0 )
    {
        maxIt = maxIt2 = 2500;
        A = new qpOASES::SparseMatrix( prob->nCon, prob->nVar, vars->jacIndRow, vars->jacIndCol, vars->jacNz );
    }
    else
        maxIt = 100;

    // Set step bounds
    lb = vars->deltaBl.ARRAY();
    lu = vars->deltaBu.ARRAY();

    lbA = vars->deltaBl.ARRAY() + prob->nVar;
    luA = vars->deltaBu.ARRAY() + prob->nVar;

    //printf("max threads = %i\n", omp_get_max_threads());
    int saveNumThreads = omp_get_max_threads();
    bool success = false;
    omp_set_num_threads( 2 );
    #pragma omp parallel for default(shared) private(i)
    for( i=0; i<maxQP; i++ )
    {
        //printf("QP %i, process number: %i\n",i+1, omp_get_thread_num());
        if( i == 0 )
        {
            if( flag == 0 )
            {
                // Convert 1st Hessian to sparse format
                vars->convertHessian( prob, param->eps, vars->hess1, vars->hessNz,
                                      vars->hessIndRow, vars->hessIndCol, vars->hessIndLo );
                H1 = new qpOASES::SymSparseMat( prob->nVar, prob->nVar,
                                               vars->hessIndRow, vars->hessIndCol, vars->hessNz );
            }

            // Call qpOASES for indefinite Hessian
            qp->setOptions( opts );
            if( qp->getStatus() == qpOASES::QPS_HOMOTOPYQPSOLVED || qp->getStatus() == qpOASES::QPS_SOLVED )
            {
                if( flag == 0 )
                    ret = qp->hotstart( H1, vars->gradObj.ARRAY(), A, lb, lu, lbA, luA, maxIt, &cpuTime );
                else if( flag == 1 ) // Second order correction: H and A do not change
                    ret = qp->hotstart( vars->gradObj.ARRAY(), lb, lu, lbA, luA, maxIt, &cpuTime );
            }
            else
            {
                ret = qp->init( H1, vars->gradObj.ARRAY(), A, lb, lu, lbA, luA, maxIt, &cpuTime );
            }

            // Statistics
            if( ret == qpOASES::RET_SETUP_AUXILIARYQP_FAILED )
                stats->qpIterations += 1;
            else
                stats->qpIterations += maxIt + 1;

            // If the first QP was solved successfully check if assumption (G3*) is satisfied
            if( ret == qpOASES::SUCCESSFUL_RETURN && maxQP == 2 )
            {
                // Remove some of the bounds. Is it still positive definite?
                ret = solAna.checkCurvatureOnStronglyActiveConstraints( qp );
            }

            if( ret != qpOASES::SUCCESSFUL_RETURN && flag == 0 )
            {
                stats->rejectedSR1++;
                stats->qpResolve++;
            }

            if( ret == qpOASES::SUCCESSFUL_RETURN )
            {
                success = true;
                #pragma omp flush (success)
            }
        }
        else if( i == 1 )
        {
            #pragma omp flush (success)
            if( !success )
            {
                // Convert 2nd Hessian to sparse format
                vars->convertHessian( prob, param->eps, vars->hess2, vars->hessNz2,
                                      vars->hessIndRow2, vars->hessIndCol2, vars->hessIndLo2 );

                H2 = new qpOASES::SymSparseMat( prob->nVar, prob->nVar,
                                               vars->hessIndRow2, vars->hessIndCol2,
                                               vars->hessNz2 );

                // Call qpOASES for positive definite Hessian
                qp2->setOptions( opts2 );
                if( qp2->getStatus() == qpOASES::QPS_HOMOTOPYQPSOLVED || qp2->getStatus() == qpOASES::QPS_SOLVED )
                {
                    ret2 = qp2->hotstart( H2, vars->gradObj.ARRAY(), A, lb, lu, lbA, luA, maxIt2, &cpuTime2 );
                }
                else
                {
                    ret2 = qp2->init( H2, vars->gradObj.ARRAY(), A, lb, lu, lbA, luA, maxIt2, &cpuTime2 );
                }

                // Statistics
                if( ret2 == qpOASES::RET_SETUP_AUXILIARYQP_FAILED )
                    stats->qpIterations2 += 1;
                else
                    stats->qpIterations2 += maxIt2 + 1;
            }
        }
    }
    omp_set_num_threads( saveNumThreads );

    // If first QP was successful, take it, else take the result from fallback
    if( ret == qpOASES::SUCCESSFUL_RETURN || maxQP == 1 )
    {
        *qp2 = *qp;
        ret2 = ret;
    }
    else
    {
        *qp = *qp2;
        ret = ret2;
    }

    // Read solution
    qp->getPrimalSolution( deltaXi.ARRAY() );
    qp->getDualSolution( lambdaQP.ARRAY() );

    // Compute constrJac*deltaXi, need this for second order correction step
    sparseAtimesb( vars->jacNz, vars->jacIndRow, vars->jacIndCol, deltaXi, vars->AdeltaXi );

    /* Print qpOASES error code (if postprocessQP is called, don't print error from first QP or SOC) */
    if( ret != qpOASES::SUCCESSFUL_RETURN && flag == 0 )
        printf( "solveQP: \"%s\"\n", qpOASES::getGlobalMessageHandler()->getErrorCodeMessage( ret ) );

    /* 0: Success
     * 1: Maximum number of iterations reached
     * 2: Unbounded
     * 3: Infeasible
     * 4: Other error */
    if( ret == qpOASES::SUCCESSFUL_RETURN )
        return 0;
    else if( ret == qpOASES::RET_MAX_NWSR_REACHED )
        return 1;
    else if( ret == qpOASES::RET_HESSIAN_NOT_SPD || ret == qpOASES::RET_HESSIAN_INDEFINITE ||
             ret == qpOASES::RET_INIT_FAILED_UNBOUNDEDNESS || ret == qpOASES::RET_QP_UNBOUNDED || ret == qpOASES::RET_HOTSTART_STOPPED_UNBOUNDEDNESS )
        return 2;
    else if( ret == qpOASES::RET_INIT_FAILED_INFEASIBILITY || ret == qpOASES::RET_QP_INFEASIBLE || ret == qpOASES::RET_HOTSTART_STOPPED_INFEASIBILITY )
        return 3;
    else
        return 4;
}


#ifndef NEWQPLOOP
/**
 * Call external solver qpOASES. The classes and methods are declared in qpOASES.hpp
 * flag = 0: First QP in a major iteration (default)
 * flag = 1: QP for second order correction step
 */
int SQPmethod::solveQP( Matrix &deltaXi, Matrix &lambdaQP, int flag )
{
    Matrix jacT;
    qpOASES::Matrix *A;
    qpOASES::SymmetricMatrix *H;
    qpOASES::Options opts;
    qpOASES::returnValue ret;
    double *lb, *lu, *lbA, *luA;
    double cpuTime;
    bool enableQPLoop = false;
    int maxIt;

    // QP loop is expensive, do it only when required:
    // if filter line search, if QP is not a SOC, first iteration is pos def by construction
    if( param->globalization == 1 && flag == 0 && stats->itCount > 1 )
        // SR1 update or BFGS without damping
        if( param->hessUpdate == 1 || param->hessUpdate == 4 || (param->hessUpdate == 2 && !param->hessDamp) )
            enableQPLoop = true;

    // set options for qpOASES
    opts.enableInertiaCorrection = qpOASES::BT_TRUE;
    opts.enableEqualities = qpOASES::BT_TRUE;
    opts.initialStatusBounds = qpOASES::ST_INACTIVE;
    opts.printLevel = qpOASES::PL_NONE;
    //opts.printLevel = qpOASES::PL_HIGH;
    opts.numRefinementSteps = 2;
    opts.epsLITests =  2.2204e-08;
    cpuTime = 10000.0;

    if( flag == 0 )
        maxIt = 2500;
    else // don't want to spend too much time for second order correction
        maxIt = 100;

    if( flag == 0 )
    {
#ifdef QPSOLVER_SPARSE
        // Convert Hessian to sparse format
        //vars->convertHessian( prob, param->eps );
        vars->convertHessian( prob, param->eps, vars->hess, vars->hessNz,
                              vars->hessIndRow, vars->hessIndCol, vars->hessIndLo );
        H = new qpOASES::SymSparseMat( prob->nVar, prob->nVar,
                                       vars->hessIndRow, vars->hessIndCol,
                                       vars->hessNz );

        A = new qpOASES::SparseMatrix( prob->nCon, prob->nVar, vars->jacIndRow, vars->jacIndCol, vars->jacNz );
#else
        // Convert diagonal block Hessian to double array
        int count = 0;
        int blockCnt = 0;
        for( int i=0; i<prob->nVar; i++ )
            for( int j=0; j<prob->nVar; j++ )
            {
                if( i == vars->blockIdx[blockCnt+1] )
                    blockCnt++;
                if( j >= vars->blockIdx[blockCnt] && j < vars->blockIdx[blockCnt+1] )
                    vars->hessNz[count++] = vars->hess[blockCnt]( i - vars->blockIdx[blockCnt], j - vars->blockIdx[blockCnt] );
                else
                    vars->hessNz[count++] = 0.0;
            }
        H = new qpOASES::SymDenseMat( prob->nVar, prob->nVar, prob->nVar, vars->hessNz );

        // transpose Jacobian (qpOASES needs row major arrays)
        Transpose( vars->constrJac, jacT );
        A = new qpOASES::DenseMatrix( prob->nCon, prob->nVar, prob->nVar, jacT.ARRAY() );
#endif
    }

    // Set step bounds
    lb = vars->deltaBl.ARRAY();
    lu = vars->deltaBu.ARRAY();

    lbA = vars->deltaBl.ARRAY() + prob->nVar;
    luA = vars->deltaBu.ARRAY() + prob->nVar;

#if (MYDEBUGLEVEL >= 3)
    stats->dumpQPCpp( prob, vars, qp );
#endif

    /// \todo For now: just duplicate the whole qp object (we probably only need the active set)
    if( enableQPLoop )
    {// Store last successful QP in temporary storage
        qpSave = *qp;
        opts.enableInertiaCorrection = qpOASES::BT_FALSE;
    }

    // Call qpOASES
    qp->setOptions( opts );
    if( qp->getStatus() == qpOASES::QPS_HOMOTOPYQPSOLVED || qp->getStatus() == qpOASES::QPS_SOLVED )
    {
        if( flag == 0 )
            ret = qp->hotstart( H, vars->gradObj.ARRAY(), A, lb, lu, lbA, luA, maxIt, &cpuTime );
        else if( flag == 1 ) // Second order correction: H and A do not change
            ret = qp->hotstart( vars->gradObj.ARRAY(), lb, lu, lbA, luA, maxIt, &cpuTime );
    }
    else
        ret = qp->init( H, vars->gradObj.ARRAY(), A, lb, lu, lbA, luA, maxIt, &cpuTime );

    // Statistics
    if( ret == qpOASES::RET_SETUP_AUXILIARYQP_FAILED )
        stats->qpIterations += 1;
    else
        stats->qpIterations += maxIt + 1;

    /* If filter line search with indefinite Hessian approximations is used
     * we may have to convexify and re-solve QP. (expensive!) */
    if( enableQPLoop )
        ret = QPLoop( opts, ret, deltaXi, lambdaQP, vars->gradObj.ARRAY(), A, lb, lu, lbA, luA );

    // Read solution
    qp->getPrimalSolution( deltaXi.ARRAY() );
    qp->getDualSolution( lambdaQP.ARRAY() );
    vars->qpObj = qp->getObjVal();

    // Compute constrJac*deltaXi, need this for second order correction step
    #ifdef QPSOLVER_SPARSE
    sparseAtimesb( vars->jacNz, vars->jacIndRow, vars->jacIndCol, deltaXi, vars->AdeltaXi );
    #else
    Atimesb( vars->constrJac, deltaXi, vars->AdeltaXi );
    #endif

    /* Print qpOASES error code (if postprocessQP is called, don't print error from first QP or SOC) */
    if( ret != qpOASES::SUCCESSFUL_RETURN && flag == 0 )
        printf( "solveQP: \"%s\"\n", qpOASES::getGlobalMessageHandler()->getErrorCodeMessage( ret ) );

    /* 0: Success
     * 1: Maximum number of iterations reached
     * 2: Unbounded
     * 3: Infeasible
     * 4: Other error */
    if( ret == qpOASES::SUCCESSFUL_RETURN )
        return 0;
    else if( ret == qpOASES::RET_MAX_NWSR_REACHED )
        return 1;
    else if( ret == qpOASES::RET_HESSIAN_NOT_SPD || ret == qpOASES::RET_HESSIAN_INDEFINITE ||
             ret == qpOASES::RET_INIT_FAILED_UNBOUNDEDNESS || ret == qpOASES::RET_QP_UNBOUNDED || ret == qpOASES::RET_HOTSTART_STOPPED_UNBOUNDEDNESS )
        return 2;
    else if( ret == qpOASES::RET_INIT_FAILED_INFEASIBILITY || ret == qpOASES::RET_QP_INFEASIBLE || ret == qpOASES::RET_HOTSTART_STOPPED_INFEASIBILITY )
        return 3;
    else
        return 4;
}
#endif


#ifdef QPSOLVER_QPOASES_SCHUR
qpOASES::returnValue SQPmethod::QPLoop( qpOASES::Options opts, qpOASES::returnValue ret, Matrix &deltaXi, Matrix &lambdaQP,
                                        double *g, qpOASES::Matrix *A, double *lb, double *lu, double *lbA, double *luA )
{
    int iQP, maxIt, hessDampSave;
    qpOASES::SymmetricMatrix *H;
    qpOASES::SolutionAnalysis solAna;
    qpOASES::returnValue retval;
    double mu1, mu2, cpuTime;

    // If the first QP was solved successfully check if assumption (G3*) is satisfied
    if( ret == qpOASES::SUCCESSFUL_RETURN )
    {
        // Remove some of the bounds. Is it still positive definite?
        if( solAna.checkCurvatureOnStronglyActiveConstraints( qp ) == qpOASES::SUCCESSFUL_RETURN )
            return qpOASES::SUCCESSFUL_RETURN;
    }

    // Statistics: save iterations from unsuccessful QP
    stats->qpIterations2 = stats->qpIterations;
    stats->qpIterations = 0;
    stats->rejectedSR1++;

    /* If QP was not successfully solved or assumption (G3*) is violated, we have to convexify and re-solve */
    for( iQP=0; iQP<param->maxConvQP; iQP++ )
    {

        mu1 = (iQP+1.0) / param->maxConvQP;
        mu2 = 1.0 - mu1;

        if( param->fallbackUpdate == 0 )
        {/* 1.) New Hessian: (1-mu)*H_SR1 + mu*(scale*I) */
            if( param->hessLimMem )
            {
                stats->itCount--;
                calcHessianUpdateLimitedMemory( param->hessUpdate, param->hessScaling, mu1 );
                stats->itCount++;
            }
            else
                printf( "Convex combination of SR1 and scaled identity not implemented for full-space Hessian!\n" );
        }
        else
        {/* 2.) New Hessian: (1-mu)*H_SR1 + mu*H_BFGS */

            // Limited memory: call update routine again, completely rebuild Hessian
            if( param->hessLimMem )
            {
                stats->itCount--;
                hessDampSave = param->hessDamp;
                param->hessDamp = 1;
                calcHessianUpdateLimitedMemory( param->fallbackUpdate, param->fallbackScaling );
                param->hessDamp = hessDampSave;
                stats->itCount++;

                if( iQP != param->maxConvQP-1 )
                {
                    // Store BFGS Hessian
                    SymMatrix *hessSave;
                    hessSave = new SymMatrix[vars->nBlocks];
                    int i, j;
                    for( int iBlock=0; iBlock<vars->nBlocks; iBlock++ )
                    {
                        hessSave[iBlock].Dimension(vars->hess[iBlock].M());
                        for( i=0; i<vars->hess[iBlock].M(); i++ )
                            for( j=i; j<vars->hess[iBlock].N(); j++ )
                                hessSave[iBlock]( i, j ) = vars->hess[iBlock]( i, j );
                    }

                    // SR1 update
                    stats->itCount--;
                    calcHessianUpdateLimitedMemory( param->hessUpdate, param->hessScaling );
                    stats->itCount++;

                    // Combine both

                    for( int iBlock=0; iBlock<vars->nBlocks; iBlock++ )
                        for( i=0; i<vars->hess[iBlock].M(); i++ )
                            for( j=i; j<vars->hess[iBlock].N(); j++ )
                                vars->hess[iBlock]( i, j ) = mu2*vars->hess[iBlock]( i, j ) +
                                                             mu1*hessSave[iBlock]( i, j );

                    //delete hessSave;
                }
            }
            else
            { // Full memory: set hess pointer to hess2, update is automatically maintained
                /// \todo convex combination is not yet implemented for the full memory case!
                vars->hess = vars->hess2;
            }
        }

        // Convert Hessian to sparse format
        //vars->convertHessian( prob, param->eps );
        vars->convertHessian( prob, param->eps, vars->hess, vars->hessNz,
                              vars->hessIndRow, vars->hessIndCol, vars->hessIndLo );
        H = new qpOASES::SymSparseMat( prob->nVar, prob->nVar, vars->hessIndRow,
                                       vars->hessIndCol, vars->hessNz );

        #if (MYDEBUGLEVEL >= 3)
        stats->dumpQPCpp( prob, vars, qp );
        #endif

        // Solve QP: Warmstart with the last successfully solved QP
        maxIt = 2500;
        cpuTime = 10000.0;
        *qp = qpSave;
        if( iQP == param->maxConvQP-1 )
        {   // Just to be safe in case there are neg eigenvalues due to ill-conditioning, switch on inertia correction
            opts.enableInertiaCorrection = qpOASES::BT_TRUE;
            qp->setOptions( opts );
        }
        retval = qp->hotstart( H, g, A, lb, lu, lbA, luA, maxIt, &cpuTime );

        // Statistics
        if( retval == qpOASES::RET_SETUP_AUXILIARYQP_FAILED )
            stats->qpIterations = 1;
        else
            stats->qpIterations = maxIt + 1;
        stats->qpResolve++;

        // Reset pointer to SR1 matrix
        if( !param->hessLimMem )
            vars->hess = vars->hess1;

        // Last Hessian is positive definite by construction, don't need to test
        if( iQP != param->maxConvQP-1 )
        {
            // Check if assumption (G3*) is satisfied
            if( retval == qpOASES::SUCCESSFUL_RETURN )
            {
                // Remove some of the bounds. Is it still positive definite?
                if( solAna.checkCurvatureOnStronglyActiveConstraints( qp ) == qpOASES::SUCCESSFUL_RETURN )
                    return qpOASES::SUCCESSFUL_RETURN;
            }
            // Statistics: save iterations from unsuccessful QP
            stats->qpIterations2 += stats->qpIterations;
            stats->qpIterations = 0;
        }
    }

    return retval;
}


qpOASES::returnValue SQPmethod::postprocessQP_Id( qpOASES::returnValue ret, Matrix &deltaXi, Matrix &lambdaQP,
                                                  qpOASES::SymmetricMatrix *H, double *g, qpOASES::Matrix *A,
                                                  double *lb, double *lu, double *lbA, double *luA )
{
    int iQP, maxIt;
    qpOASES::SolutionAnalysis solAna;
    qpOASES::returnValue retval;
    stats->qpResolve = 0;
    double deltaH = 0.0, cpuTime;

    // If first QP was solved successfully check if we have to convexify and solve again to satisfy assumption (G3*)
    if( ret == qpOASES::SUCCESSFUL_RETURN )
    {
        // Remove some of the bounds. Is it still positive definite?
        if( solAna.checkCurvatureOnStronglyActiveConstraints( qp ) == qpOASES::SUCCESSFUL_RETURN )
            return qpOASES::SUCCESSFUL_RETURN;
    }

    /* If QP was not successfully solved or assumption (G3*) is violated, we have to convexify and re-solve */
    for( iQP=0; iQP<param->maxConvQP; iQP++ )
    {
        // Set convexification parameter
        if( iQP == 0 )
        {
            // If deltaH has not been set before, we have to start somewhere
            if( vars->deltaH( 0 ) == 0.0 )
                deltaH = param->deltaH0;
            else// If we have convexified before, use a fraction of the last used value as first guess
                deltaH = fmax( 10*param->eps, param->kappaMinus * vars->deltaH( 0 ) );
        }
        else
        {
            // If we have no history for deltaH, increase a lot
            if( vars->deltaH( 0 ) == 0.0 )
                deltaH *= param->kappaPlusMax;
            else// If we have some history for deltaH, increase a little
                deltaH *= param->kappaPlus;
        }

        // Modify Hessian
        ///\note maybe each block could be convexified by a different parameter
        H->addToDiag( deltaH );
        //printf("deltaH = %23.16e\n", deltaH );

        #if (MYDEBUGLEVEL >= 3)
        stats->dumpQPCpp( prob, vars, qp );
        #endif

        // Solve QP
        maxIt = 2500;
        cpuTime = 10000.0;

        // Warmstart with the last successfully solved QP
        *qp = qpSave;
        retval = qp->hotstart( H, g, A, lb, lu, lbA, luA, maxIt, &cpuTime );

        if( ret == qpOASES::RET_SETUP_AUXILIARYQP_FAILED )
            stats->qpIterations2 += 1;
        else
            stats->qpIterations2 += maxIt + 1;
        stats->qpResolve++;

        if( retval == qpOASES::SUCCESSFUL_RETURN )
        {
            // Remove some of the bounds. Is it still positive definite?
            if( solAna.checkCurvatureOnStronglyActiveConstraints( qp ) == 0 )
            {
                // Subtract the regularization before exiting
                H->addToDiag( -deltaH );
                vars->deltaH( 0 ) = deltaH;
                return qpOASES::SUCCESSFUL_RETURN;
            }
        }
        //else
            //printf( "postProcessQP(), it %i: mu = %g. qpOASES errormessage \"%s\"\n", iQP, deltaH, qpOASES::getGlobalMessageHandler()->getErrorCodeMessage( retval ) );
    }

    return retval;
}
#endif


/**
 * Check if QP returns a solution and analyze inertia of reduced Hessian
 * (see filter line search Waechter paper, section 4.2)
 */
#if defined QPSOLVER_QPOASES_DENSE || defined QPSOLVER_QPOASES_SPARSE
qpOASES::returnValue SQPmethod::QPLoop( qpOASES::Options opts, qpOASES::returnValue ret, Matrix &deltaXi, Matrix &lambdaQP,
                                        double *g, qpOASES::Matrix *A, double *lb, double *lu, double *lbA, double *luA )
{
    printf("QPLoop not implemented for dense qpOASES\n");
    return qpOASES::SUCCESSFUL_RETURN;
}

qpOASES::returnValue SQPmethod::postprocessQP_Id( qpOASES::returnValue ret, Matrix &deltaXi, Matrix &lambdaQP,
                                                  qpOASES::SymmetricMatrix *H, double *g, qpOASES::Matrix *A,
                                                  double *lb, double *lu, double *lbA, double *luA )
{
    printf("postprocessQP_ID not implemented for dense qpOASES.\n");
    return qpOASES::SUCCESSFUL_RETURN;
}
#endif


/**
 * Set bounds on the step (in the QP), either according
 * to variable bounds in the NLP or according to
 * trust region box radius
 */
void SQPmethod::updateStepBounds( bool soc )
{
    int i;
    int nVar = prob->nVar;
    int nCon = prob->nCon;

    // Bounds on step
    for( i=0; i<nVar; i++ )
    {
        if( prob->bl(i) != param->inf )
            vars->deltaBl( i ) = prob->bl( i ) - vars->xi( i );
        else
            vars->deltaBl( i ) = param->inf;

        if( prob->bu(i) != param->inf )
            vars->deltaBu( i ) = prob->bu( i ) - vars->xi( i );
        else
            vars->deltaBu( i ) = param->inf;
    }

    // Bounds on linearized constraints
    for( i=0; i<nCon; i++ )
    {
        if( prob->bl( nVar+i ) != param->inf )
        {
            vars->deltaBl( nVar+i ) = prob->bl( nVar+i ) - vars->constr( i );
            if( soc ) vars->deltaBl( nVar+i ) += vars->AdeltaXi( i );
        }
        else
            vars->deltaBl( nVar+i ) = param->inf;

        if( prob->bu( nVar+i ) != param->inf )
        {
            vars->deltaBu( nVar+i ) = prob->bu( nVar+i ) - vars->constr( i );
            if( soc ) vars->deltaBu( nVar+i ) += vars->AdeltaXi( i );
        }
        else
            vars->deltaBu( nVar+i ) = param->inf;
    }
}

} // namespace blockSQP


