#include "blocksqp.hpp"
#include "blocksqp_general_purpose.hpp"

namespace blockSQP
{

/**
 * Call external solver qpOASES. The classes and methods are declared in qpOASES.hpp
 * flag = 0: First QP in a major iteration (default)
 * flag = 1: QP for second order correction step
 */
#ifdef QPSOLVER_QPOASES
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
        if( param->hessUpdate == 1 || (param->hessUpdate == 2 && !param->hessDamp) )
            enableQPLoop = true;

    // set options for qpOASES
    opts.enableInertiaCorrection = qpOASES::BT_TRUE;
    opts.enableEqualities = qpOASES::BT_TRUE;
    //opts.initialStatusBounds = qpOASES::ST_UPPER;
    opts.initialStatusBounds = qpOASES::ST_INACTIVE;
    opts.printLevel = qpOASES::PL_NONE;
    opts.numRefinementSteps = 2;
    opts.epsLITests =  2.2204e-08;
    cpuTime = 10000.0;

    if( flag == 0 )
        maxIt = 5000;
    else // don't want to spend too much time for second order correction
        maxIt = 100;

    if( flag == 0 )
    {
#ifdef QPSOLVER_SPARSE
        // Convert Hessian to sparse format
        vars->convertHessian( prob );
        H = new qpOASES::SymSparseMat( prob->nVar, prob->nVar,
                                       vars->hessIndRow, vars->hessIndCol,
                                       vars->hessNz, vars->hessIndLo );

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

#ifdef MYDEBUG
    //stats->dumpQPCpp( prob, vars, qp );
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
        else if( flag == 1 )
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
    int iQP, maxQP = 1, maxIt, hessDampSave;
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
    for( iQP=0; iQP<maxQP; iQP++ )
    {

        mu1 = (iQP+1.0) / maxQP;
        mu2 = 1.0 - mu1;

        /* 1.) New Hessian: (1-mu)*H_SR1 + mu*(scale*I) */
        //stats->itCount--;
        //calcHessianUpdateLimitedMemory( param->hessUpdate, param->hessScaling, mu1 );
        //stats->itCount++;

        /* 2.) New Hessian: (1-mu)*H_SR1 + mu*H_BFGS */
        // fallback: damped BFGS update

        // Limited memory: call update routine again, completely rebuild Hessian
        if( param->hessLimMem )
        {
            stats->itCount--;
            hessDampSave = param->hessDamp;
            param->hessDamp = 1;
            calcHessianUpdateLimitedMemory( param->fallbackUpdate, param->fallbackScaling );
            param->hessDamp = hessDampSave;
            stats->itCount++;

            if( iQP != maxQP-1 )
            {
                // Store BFGS Hessian
                SymMatrix hessSave[vars->nBlocks];
                for( int iBlock=0; iBlock<vars->nBlocks; iBlock++ )
                    hessSave[iBlock] = SymMatrix( vars->hess[iBlock] );

                // SR1 update
                stats->itCount--;
                calcHessianUpdateLimitedMemory( param->hessUpdate, param->hessScaling );
                stats->itCount++;

                // Combine both
                int i, j;
                for( int iBlock=0; iBlock<vars->nBlocks; iBlock++ )
                    for( i=0; i<vars->hess[iBlock].M(); i++ )
                        for( j=i; j<vars->hess[iBlock].N(); j++ )
                            vars->hess[iBlock]( i, j ) = mu2*vars->hess[iBlock]( i, j ) +
                                                         mu1*hessSave[iBlock]( i, j );
            }
        }
        else
        { // Full memory: set hess pointer to hess2, update is automatically maintained

            /// \todo convex combination is not yet implemented for the full memory case!
            vars->hess = vars->hess2;
        }

        // Convert Hessian to sparse format
        vars->convertHessian( prob );
        H = new qpOASES::SymSparseMat( prob->nVar, prob->nVar, vars->hessIndRow,
                                       vars->hessIndCol, vars->hessNz, vars->hessIndLo );

        #ifdef MYDEBUG
        //stats->dumpQPCpp( prob, vars, qp );
        #endif

        // Solve QP: Warmstart with the last successfully solved QP
        maxIt = 3000;
        cpuTime = 10000.0;
        *qp = qpSave;
        if( iQP == maxQP-1 )
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
        if( iQP != maxQP-1 )
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
    int iQP, maxQP = 10, maxIt;
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
    for( iQP=0; iQP<maxQP; iQP++ )
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

        #ifdef MYDEBUG
        //stats->dumpQPCpp( prob, vars, qp );
        #endif

        // Solve QP
        maxIt = 3000;
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
    //qpOASES::returnValue retval;
    //retval = ret;

    //int iQP, iBlock, i, j, k, l, count;
    //double cpuTime;
    //double averageMu;
    //int maxIt;

    //int n = deltaXi.M();
    //int m = lambdaQP.M() - n;
    //int mBig;
    //int dimNull;
    //int nFixed;
    //int *xFixed;
    //double *Z, *ZHZ;
    //double *work;
    //int info, lwork, *ipiv, nNegEig;

    //stats->qpResolve = 0;

    //int strategy = 1;
    //int maxQP = 5;
    //Matrix deltaH;
    //deltaH.Dimension( vars->nBlocks );

    ///* Compute (possible) regularization paramter to convexify QP */
    //if( strategy == 1 )
    //{// Estimate the smallest eigenvalue and use a fraction of this estimate to convexify
        //deltaH.Initialize( 0.0 );
        //for( iBlock=0; iBlock<vars->nBlocks; iBlock++ )
            //deltaH( iBlock ) = -estimateSmallestEigenvalue( vars->hess[iBlock] ) / maxQP;
    //}

    ///// \todo if all mu positive, we don't need to do anything
    //for( iQP=0; iQP<maxQP; iQP++ )
    //{
        //if( retval == qpOASES::SUCCESSFUL_RETURN )
        //{// If QP was solved successfully check if we have to convexify and solve again to satisfy assumption (G3*)

            ///* Read solution */
            //qp->getPrimalSolution( deltaXi.ARRAY() );
            //qp->getDualSolution( lambdaQP.ARRAY() );

            ///* Compute set S_k = {fixed xi AND deltaXi} */
            //nFixed = 0;
            //for( i=0; i<n; i++ )
                //if( ( fabs( vars->xi( i ) - prob->bl( i ) ) < 10*param->eps && fabs( vars->deltaXi( i ) - vars->deltaBl( i ) ) < 10*param->eps  ) ||
                    //( fabs( vars->xi( i ) - prob->bu( i ) ) < 10*param->eps && fabs( vars->deltaXi( i ) - vars->deltaBu( i ) ) < 10*param->eps  ) )
                    //nFixed++;// only if variable and step are at their bounds, the variable is considered fixed

            //count = 0;
            //xFixed = new int[nFixed];
            //for( i=0; i<n; i++ )
                //if( ( fabs( vars->xi( i ) - prob->bl( i ) ) < 10*param->eps && fabs( vars->deltaXi( i ) - vars->deltaBl( i ) ) < 10*param->eps  ) )
                    //xFixed[count++] = i+1;// only if variable and step are at their bounds, the variable is considered fixed
                //else if( ( fabs( vars->xi( i ) - prob->bu( i ) ) < 10*param->eps && fabs( vars->deltaXi( i ) - vars->deltaBu( i ) ) < 10*param->eps  ) )
                    //xFixed[count++] = -(i+1);// only if variable and step are at their bounds, the variable is considered fixed

            //// DEBUG: Print set S_k
            ////for(i=0;i<nFixed;i++)printf("xFixed[%i]=%i\n",i,xFixed[i]);

            ///* Form Matrix AS = ( A**T, I_S ) */
            //mBig = m+nFixed;
            //dimNull = n - mBig;
            //double *AS = new double[n*mBig];
            //for( i=0; i<n; i++ )
            //{
                //for( j=0; j<m; j++ )
                    //AS[i+j*n] = A[i+j*n];
                //for( j=m; j<mBig; j++ )
                    //if( xFixed[j-m] == i+1 )
                        //AS[i+j*n] = 1.0;
                    //else if( xFixed[j-m] == -(i+1) )
                        //AS[i+j*n] = -1.0;
                    //else
                        //AS[i+j*n] = 0.0;
            //}

            //// DEBUG: Print AS
            ////double *AS2 = new double[n*mBig];
            ////for( i=0; i<n*mBig; i++ )
                ////AS2[i] = AS[i];
            ////for( i=0; i<n; i++ )
            ////{
                ////for( j=0; j<mBig; j++ )
                    ////printf( "%23.16e ", AS2[i+j*n] );
                ////printf("\n");
            ////}

            ///* Compute a nullspace basis Z for AS */
            //double *tau = new double[mBig];
            //lwork = mBig*mBig+dimNull*dimNull;
            //work = new double[lwork];
            //info = 0;

            //// QR decomposition of AS**T (AS will be overwritten)
            //dgeqrf( &n, &mBig, AS, &n, tau, work, &lwork, &info );
            //if( info )
                //printf("Error: info = %i in DGEQRF\n", info );

            //// We want the last (n-mBig) columns of the nxn matrix Q
            ///// \note This will probably fail if the problem has inactive inequality constraints, so make sure it doesn't!
            ////printf("n=%i, m=%i, mBig=%i\n", n, m, mBig);
            //Z = new double[n*dimNull];
            //for( i=0; i<n*dimNull; i++ )
                //Z[i] = 0.0;
            //for( i=0; i<dimNull; i++ )
                //Z[i+mBig + i*n] = 1.0;

            //// Z := Q * (0 I)**T
            //dormqr( "l", "n", &n, &dimNull, &mBig, AS, &n, tau, Z, &n,
                    //work, &lwork, &info, strlen("l"), strlen("n") );
            //if( info )
                //printf("Error: info = %i in DORMQR\n", info );
            //delete[] work;
            //delete[] tau;

            //// DEBUG: print AS**T * Z. Should be zero!
            ////double *test = new double[mBig*dimNull];
            ////for( k=0; k<mBig*dimNull; k++ )
                ////test[k] = 0.0;
            ////for( i=0; i<mBig; i++ )
                ////for( j=0; j<dimNull; j++ )
                    ////for( k=0; k<n; k++ )
                        ////test[i+j*mBig] += AS2[k+i*n]*Z[k+j*n];
            ////for( i=0; i<mBig; i++ )
            ////{
                ////for( j=0; j<dimNull; j++ )
                    ////printf("%23.16e ", test[i+j*mBig]);
                ////printf("\n");
            ////}
            ////printf("\n");
            ////delete[] AS2;
            ////delete[] test;

            //delete[] AS;

            ///* Form Z'*H*Z */
            //ZHZ = new double[dimNull*dimNull];
            //for( i=0; i<dimNull; i++ )
                //for( j=0; j<dimNull; j++ )
                //{
                    //ZHZ[i+j*dimNull] = 0.0;
                    //for( k=0; k<n; k++ )
                        //for( l=0; l<n; l++ )
                            //ZHZ[i+j*dimNull] += Z[k+i*n] * H[k+l*n] * Z[l+j*n];
                //}
            //delete[] xFixed;
            //delete[] Z;

            //// DEBUG: print Z'HZ
            ////for( i=0; i<dimNull; i++ )
            ////{
                ////for( j=0; j<dimNull; j++ )
                    ////printf( "%23.16e ", ZHZ[i+j*dimNull] );
                ////printf("\n");
            ////}

            ///* Compute inertia of Z'*H*Z */
            //nNegEig = 0;
            //if( dimNull > 0 )
            //{
                //lwork = dimNull*dimNull+1;
                //work = new double[lwork];
                //ipiv = new int[dimNull];
                //dsytrf( "u", &dimNull, ZHZ, &dimNull, ipiv, work, &lwork, &info, strlen("u") );
                //if( info )
                    //printf("Warning: info = %i in DSYTRF\n", info );
                //delete[] work;
                //delete[] ipiv;
                //for( i=0; i<dimNull; i++ )
                    //if( ZHZ[i+i*dimNull] < -param->eps )
                        //nNegEig++;
                //delete[] ZHZ;
            //}

            //// If all eigenvalues positive, assumption (G3*) from Waechter paper should be satisfied and we're done
            //if( nNegEig == 0 )
            //{
                //for( iBlock=0; iBlock<vars->nBlocks; iBlock++ )
                    //vars->deltaH( iBlock ) = deltaH( iBlock );
                //return retval;
            //}
            //else
            //{// If not form B + deltaH*I and solve again
                //for( iBlock=0; iBlock<vars->nBlocks; iBlock++ )
                //{
                    //// Set convexification parameter
                    //if( strategy == 0 && iQP == 0 )
                    //{
                        //// If deltaH has not been set before
                        //if( vars->deltaH( iBlock ) == 0.0 )
                            //deltaH( iBlock ) = param->deltaH0;
                        //else
                            //deltaH( iBlock ) = fmax( 10*param->eps, param->kappaMinus * vars->deltaH( iBlock ) );
                    //}
                    //else if( strategy == 0 )
                    //{
                        //// If we have no history for deltaH, increase a lot
                        //if( vars->deltaH( iBlock ) == 0.0 )
                            //deltaH( iBlock ) = param->kappaPlusMax * deltaH( iBlock );
                        //else// If we have some history for deltaH, increase a little
                            //deltaH( iBlock ) = param->kappaPlus * deltaH( iBlock );
                    //}

                    //for( i=prob->blockIdx[iBlock]; i<prob->blockIdx[iBlock+1]; i++ )
                        //H[i+i*n] += deltaH( iBlock );
                //}
            //}
        //}
        //else
        //{// If QP is unbounded or reduced Hessian is indefinite, assumption (G3*) is violated in any case and we have to convexify and resolve

            //// Set convexification parameter
            //if( strategy == 0 && iQP == 0 )
            //{
                //// If deltaH has not been set before
                //if( vars->deltaH( iBlock ) == 0.0 )
                    //deltaH( iBlock ) = param->deltaH0;
                //else
                    //deltaH( iBlock ) = fmax( 10*param->eps, param->kappaMinus * vars->deltaH( iBlock ) );
            //}
            //else if( strategy == 0 )
            //{
                //// If we have no history for deltaH, increase a lot
                //if( vars->deltaH( iBlock ) == 0.0 )
                    //deltaH( iBlock ) = param->kappaPlusMax * deltaH( iBlock );
                //else// If we have some history for deltaH, increase a little
                    //deltaH( iBlock ) = param->kappaPlus * deltaH( iBlock );
            //}

            //// Form B + mu*I
            //for( iBlock=0; iBlock<vars->nBlocks; iBlock++ )
                //for( i=prob->blockIdx[iBlock]; i<prob->blockIdx[iBlock+1]; i++ )
                    //H[i+i*n] += deltaH( iBlock );
        //}

//#ifdef MYDEBUG
        //stats->dumpQPCpp( prob, vars, qp );
//#endif

        //// Solve QP
        //maxIt = 5000;
        //cpuTime = 10000.0;

        //// If previous QP could not be solved we cannot use hotstart
        //if( retval == 0 )
            //retval = qp->hotstart( H, g, A, lb, lu, lbA, luA, maxIt, &cpuTime );
        //else
        //{
            //// Use active set from the previous major (successful) iteration, this might be a better guess than from the failed QP
            ///// \note sometimes that does not seem to work well at all, why?
            //double *xOpt = 0;
            //double *yOpt = 0;
            //retval = qp->init( H, g, A, lb, lu, lbA, luA, maxIt, &cpuTime, xOpt, yOpt, &b, &c );
        //}

        //if( retval != qpOASES::SUCCESSFUL_RETURN )
            //printf( "%s\n", qpOASES::getGlobalMessageHandler()->getErrorCodeMessage( retval ) );

        //stats->qpIterations2 += (maxIt+1);
        //stats->qpResolve++;
    //}
    //for( iBlock=0; iBlock<vars->nBlocks; iBlock++ )
        //vars->deltaH( iBlock ) = deltaH( iBlock );

    //return retval;
}
#endif


#ifdef QPSOLVER_QPOPT
/**
 * Generic static function (only visible in this file)
 * for Hessian-vector products assuming Hessian is stored in column-compressed format
 */
static void qpoptHess( int *n, int *ldH, int *jthcol, double *H, double *x, double *Hx,
                       int *iw, int *leniw, double *w, int *lenw )
{
    int i, k;

    // Hessian is stored in column compressed sparse format in H with index information in ldH
    int nnz = ldH[*n];
    int *hessIndCol = ldH;
    int *hessIndRow = ldH - nnz;

    // Compute actually x^T * H
    for( i=0; i<*n; i++ )
    {
        Hx[i] = 0.0;
        // k runs over all elements in one column
        for( k=hessIndCol[i]; k<hessIndCol[i+1]; k++ )
            Hx[i] += H[k] * x[hessIndRow[k]];
    }
}

/**
 * Call external solver QPOPT. The calling routine is declared in qpopt.h
 */
int SQPmethod::solveQP( Matrix &deltaXi, Matrix &lambdaQP )
{
    int n, nclin, i, info, iter, ioptns, ios, ibuffer;
    int *istate;
    double obj, *clambda;
    char* printfile;
    char* buffer;

    double *A, *H;
    int ldA, *ldH;

    int *iw, leniw, lenw;
    double *w;

    // Use QPOPT nomenclature to avoid confusion
    n = prob->nVar;
    nclin = prob->nCon;
    ldA = nclin;
    A = vars->constrJac.ARRAY();
    clambda = lambdaQP.ARRAY();

    // Convert Hessian to sparse format
    vars->convertHessian( prob );

    // ldH should in fact be an int variable, I try to get the Hessian
    // index information into my Hess subroutine with this
    ldH = vars->hessIndCol;
    H = vars->hessNz;

    // Needs to be set for warm starts
    istate = vars->istate;

    // workspace
    leniw = 2*n + 3 + 10000;
    iw = new int[leniw];
    lenw = 2*n*n + 8*n + 5*nclin + 10000;
    w = new double[lenw];

    // Set options for qpopt
    //ioptns = 14;
    //sprintf( printfile, "%s%s", stats->outpath, "qpopt.opt" );
    //fortopen( &ioptns, printfile, "unknown", &ios,
              //strlen(printfile), strlen("unknown") );

    //qpprms( &ioptns, &info );
    //fortclose( &ioptns );

    // Set Options
    strcpy( buffer, "Print Level 0\0" );
    ibuffer = 0;
    qpprm( buffer, 13);
    // Maximum number of iterations
    strcpy( buffer, "It\0" );
    ibuffer = 10000;
    qpprmi( buffer, &ibuffer, 2 );

    qpopt( &n, &nclin, &ldA, ldH, A, vars->deltaBl.ARRAY(), vars->deltaBu.ARRAY(),
            vars->gradObj.ARRAY(), H, qpoptHess, istate, deltaXi.ARRAY(),
            &info, &iter, &obj, vars->AdeltaXi.ARRAY(), clambda,
            iw, &leniw, w, &lenw );

    stats->qpIterations += iter;
    vars->qpObj = obj;

    delete[] iw;
    delete[] w;

    /* 0: Success
     * 1: Maximum number of iterations reached
     * 2: Unbounded
     * 3: Infeasible
     * 4: Other error */
    if( info == 0 || info == 1 )
        return 0;
    else if( info == 4 )
        return 1;
    else if( info == 2 )
        return 2;
    else if( info == 3 )
        return 3;
    else
        return info;
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


