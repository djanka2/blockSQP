#include "blocksqp.hpp"
#include "blocksqp_general_purpose.hpp"

namespace blockSQP
{

/**
 * Initial Hessian: Identity matrix
 */
void SQPmethod::calcInitialHessian()
{
    int iBlock;

    for( iBlock=0; iBlock<vars->nBlocks; iBlock++ )
        //if objective derv is computed exactly, don't set the last block!
        if( !(param->whichSecondDerv == 1 && param->blockHess && iBlock == vars->nBlocks-1) )
            calcInitialHessian( iBlock );
}


/**
 * Initial Hessian for one block: Identity matrix
 */
void SQPmethod::calcInitialHessian( int iBlock )
{
    vars->hess[iBlock].Initialize( 0.0 );

    // Each block is a diagonal matrix
    for( int i=0; i<vars->hess[iBlock].M(); i++ )
        vars->hess[iBlock]( i, i ) = param->iniHessDiag;
}


void SQPmethod::resetHessian()
{
    for( int iBlock=0; iBlock<vars->nBlocks; iBlock++ )
        //if objective derv is computed exactly, don't set the last block!
        if( !(param->whichSecondDerv == 1 && param->blockHess && iBlock == vars->nBlocks - 1) )
            resetHessian( iBlock );
}


void SQPmethod::resetHessian( int iBlock )
{
    Matrix smallDelta, smallGamma;
    int nVarLocal = vars->hess[iBlock].M();

    // smallGamma and smallDelta are either subvectors of gamma and delta
    // or submatrices of gammaMat, deltaMat, i.e. subvectors of gamma and delta from m prev. iterations (for L-BFGS)
    smallGamma.Submatrix( vars->gammaMat, nVarLocal, vars->gammaMat.N(), vars->blockIdx[iBlock], 0 );
    smallDelta.Submatrix( vars->deltaMat, nVarLocal, vars->deltaMat.N(), vars->blockIdx[iBlock], 0 );

    // Remove past information on Lagrangian gradient difference
    smallGamma.Initialize( 0.0 );

    // Remove past information on steps
    smallDelta.Initialize( 0.0 );

    // Remove information on old scalars (used for COL sizing)
    vars->deltaNormOld( iBlock ) = 1.0;
    vars->deltaGammaOld( iBlock ) = 0.0;

    vars->noUpdateCounter[iBlock] = -1;

    calcInitialHessian( iBlock );
}

/**
 * Approximate Hessian by finite differences
 */
int SQPmethod::calcFiniteDiffHessian()
{
    int iVar, jVar, iBlock, nBlocks, infoConstr, infoObj, info, offset, nVarBlock;
    double dummy;
    SymMatrix *dummySym;
    SQPiterate *varsP1, *varsP2;

    const double DELTA = 1.0e-8;

    varsP1 = new SQPiterate( vars );
    varsP2 = new SQPiterate( vars );

    info = 0;
    for( iBlock=0; iBlock<vars->nBlocks; iBlock++ )
    {
        offset = vars->blockIdx[iBlock];
        nVarBlock = vars->blockIdx[iBlock+1] - vars->blockIdx[iBlock];
        vars->hess[iBlock].Initialize( 0.0 );
        for( iVar=vars->blockIdx[iBlock]; iVar<vars->blockIdx[iBlock+1]; iVar++ )
        {
            vars->xi( iVar ) += DELTA;

            // Evaluate objective and constraints for perturbed xi and compute perturbed Lagrange gradient
            #ifdef QPSOLVER_SPARSE
            prob->evaluate( vars->xi, vars->lambda, &dummy, varsP1->constr, varsP1->gradObj, varsP1->jacNz, varsP1->jacIndRow, varsP1->jacIndCol, vars->hess, 1, &info );
            calcLagrangeGradient( vars->lambda, varsP1->gradObj, varsP1->jacNz, varsP1->jacIndRow, varsP1->jacIndCol, varsP1->gradLagrange, 0 );
            #else
            prob->evaluate( vars->xi, vars->lambda, &dummy, varsP1->constr, varsP1->gradObj, varsP1->constrJac, vars->hess, 1, &info );
            calcLagrangeGradient( vars->lambda, varsP1->gradObj, varsP1->constrJac, varsP1->gradLagrange, 0 );
            #endif

            vars->xi( iVar ) -= 2*DELTA;

            // Evaluate objective and constraints for perturbed xi and compute perturbed Lagrange gradient
            #ifdef QPSOLVER_SPARSE
            prob->evaluate( vars->xi, vars->lambda, &dummy, varsP2->constr, varsP2->gradObj, varsP2->jacNz, varsP2->jacIndRow, varsP2->jacIndCol, vars->hess, 1, &info );
            calcLagrangeGradient( vars->lambda, varsP2->gradObj, varsP2->jacNz, varsP2->jacIndRow, varsP2->jacIndCol, varsP2->gradLagrange, 0 );
            #else
            prob->evaluate( vars->xi, vars->lambda, &dummy, varsP2->constr, varsP2->gradObj, varsP2->constrJac, vars->hess, 1, &info );
            calcLagrangeGradient( vars->lambda, varsP2->gradObj, varsP2->constrJac, varsP2->gradLagrange, 0 );
            #endif

            vars->xi( iVar ) += DELTA;

            // Compute central finite difference approximation
            for( jVar=iVar; jVar<vars->blockIdx[iBlock+1]; jVar++ )
                vars->hess[iBlock]( iVar-offset, jVar-offset ) = ( varsP1->gradLagrange( jVar ) - varsP2->gradLagrange( jVar ) ) / (2.0 * DELTA);
        }
    }

    delete varsP1, varsP2;

    return info;
}


void SQPmethod::sizeHessianNocedal( const Matrix &gamma, const Matrix &delta, int iBlock )
{
    int i, j;
    double scale;
    double myEps = 1.0e2 * param->eps;

    scale = adotb( delta, gamma );
    if( scale > 0.0 )
    {
        scale = adotb( gamma, gamma ) / fmax( scale, myEps );
        for( i=0; i<vars->hess[iBlock].M(); i++ )
            for( j=i; j<vars->hess[iBlock].M(); j++ )
                vars->hess[iBlock]( i,j ) *= scale;
    }
    else
        scale = 1.0;

    // statistics: average sizing factor
    stats->averageSizingFactor += scale;
}

void SQPmethod::sizeHessianMean( const Matrix &gamma, const Matrix &delta, int iBlock )
{
    int i, j;
    double scale;
    double myEps = 1.0e2 * param->eps;

    scale = sqrt( adotb( gamma, gamma ) / fmax( adotb( delta, delta ), myEps ) );
    for( i=0; i<vars->hess[iBlock].M(); i++ )
        for( j=i; j<vars->hess[iBlock].M(); j++ )
            vars->hess[iBlock]( i,j ) *= scale;

    // statistics: average sizing factor
    stats->averageSizingFactor += scale;
}

void SQPmethod::sizeHessianOL( const Matrix &gamma, const Matrix &delta, int iBlock )
{
    if( gamma.M() == 0 )
        return;
    int i, j;
    double value1, value2, value3;
    double myEps = 1.0e2 * param->eps;

    value1 = adotb( gamma, delta );
    value2 = adotb( delta, delta );

    /// \todo still has the form from COL sizing with all the safeguards. Should look like Nocedal
    /// and Mean after computations are done! (maybe put them into one function)
    if( value2 > myEps )
        value3 = value1 / value2;
    else
        value3 = 1.0;

    if( value3 > 0.0 && value3 < 1.0 )
    {
        value3 = fmax( param->colEps, value3 );
        for( i=0; i<vars->hess[iBlock].M(); i++ )
            for( j=i; j<vars->hess[iBlock].M(); j++ )
                vars->hess[iBlock]( i,j ) *= value3;
    }
    else
        value3 = 1.0;

    // statistics: average sizing factor
    stats->averageSizingFactor += value3;
}

void SQPmethod::sizeHessianTapia( const Matrix &gamma, const Matrix &delta, int iBlock )
{
    int i, j;
    double theta, scale, myEps = 1.0e2 * param->eps;
    double deltaNorm, deltaGamma, deltaBdelta;

    deltaNorm = adotb( delta, delta );
    deltaGamma = adotb( delta, gamma );
    deltaBdelta = 0.0;
    for( i=0; i<delta.M(); i++ )
        for( j=0; j<delta.M(); j++ )
            deltaBdelta += delta( i ) * vars->hess[iBlock]( i, j ) * delta( j );

    // Centered Oren-Luenberger factor
    if( vars->noUpdateCounter[iBlock] == -1 ) // in the first iteration, this should equal the OL factor
        theta = 1.0;
    else
        theta = fmin( param->colTau1, param->colTau2 * deltaNorm );
    if( deltaNorm > myEps && vars->deltaNormOld(iBlock) > myEps )
        scale = ( (1.0 - theta)*vars->deltaGammaOld(iBlock) / vars->deltaNormOld(iBlock) + theta*deltaGamma / deltaNorm ) /
            ( (1.0 - theta)*vars->deltaGammaOld(iBlock) / vars->deltaNormOld(iBlock) + theta*deltaBdelta / deltaNorm );
        //if( (scale = (1.0 - theta)*vars->deltaGammaOld(iBlock) / vars->deltaNormOld(iBlock) + theta*vars->deltaBdelta(iBlock) / vars->deltaNorm(iBlock)) > myEps )
            //scale = ( (1.0 - theta)*vars->deltaGammaOld(iBlock) / vars->deltaNormOld(iBlock) + theta*vars->deltaGamma(iBlock) / vars->deltaNorm(iBlock) ) / scale;
    else
        scale = 1.0;

    // Size only if factor is between zero and one
    if( scale < 1.0 && scale > 0.0 )
    {
        scale = fmax( param->colEps, scale );
        //printf("Sizing value (Tapia) block %i = %g\n", iBlock, scale );
        for( i=0; i<vars->hess[iBlock].M(); i++ )
            for( j=i; j<vars->hess[iBlock].M(); j++ )
                vars->hess[iBlock]( i,j ) *= scale;

        // statistics: average sizing factor
        stats->averageSizingFactor += scale;
    }
    else
        stats->averageSizingFactor += 1.0;

    // Save deltaNorm and deltaGamma for the next step
    vars->deltaNormOld(iBlock) = deltaNorm;
    vars->deltaGammaOld(iBlock) = deltaGamma;
}

/**
 * Apply BFGS or SR1 update blockwise and size blocks
 */
void SQPmethod::calcHessianUpdate( int updateType, int hessScaling )
{
    int iBlock, nBlocks;
    int nVarLocal;
    Matrix smallGamma, smallDelta;
    bool firstIter;

    //if objective derv is computed exactly, don't set the last block!
    if( param->whichSecondDerv == 1 && param->blockHess )
        nBlocks = vars->nBlocks - 1;
    else
        nBlocks = vars->nBlocks;

    // Statistics: how often is damping active, what is the average COL sizing factor?
    stats->hessDamped = 0;
    stats->averageSizingFactor = 0.0;

    for( iBlock=0; iBlock<nBlocks; iBlock++ )
    {
        nVarLocal = vars->hess[iBlock].M();

        // smallGamma and smallDelta are subvectors of gamma and delta, corresponding to partially separability
        smallGamma.Submatrix( vars->gammaMat, nVarLocal, vars->gammaMat.N(), vars->blockIdx[iBlock], 0 );
        smallDelta.Submatrix( vars->deltaMat, nVarLocal, vars->deltaMat.N(), vars->blockIdx[iBlock], 0 );

        // Is this the first iteration or the first after a Hessian reset?
        firstIter = ( vars->noUpdateCounter[iBlock] == -1 );

        // Sizing before the update
        if( hessScaling == 1 && firstIter )
            sizeHessianNocedal( smallGamma, smallDelta, iBlock );
        else if( hessScaling == 2 && firstIter )
            sizeHessianOL( smallGamma, smallDelta, iBlock );
        else if( hessScaling == 3 && firstIter )
            sizeHessianMean( smallGamma, smallDelta, iBlock );
        else if( hessScaling == 4 )
            sizeHessianTapia( smallGamma, smallDelta, iBlock );

        // Compute the new update
        if( updateType == 1 )
            calcSR1( smallGamma, smallDelta, iBlock );
        else if( updateType == 2 )
            calcBFGS( smallGamma, smallDelta, iBlock );

        // If an update is skipped to often, reset Hessian block
        if( vars->noUpdateCounter[iBlock] > param->maxConsecSkippedUpdates )
            resetHessian( iBlock );
    }

    // statistics: average sizing factor
    stats->averageSizingFactor /= nBlocks;
}


int SQPmethod::sizeHessianByrdLu( const Matrix &gammaMat, const Matrix &deltaMat, int iBlock )
{
    Matrix W, Dinv, N, C, gammai, deltai;
    int i, j, k, m, posDelta, posGamma, posOldest, posNewest;
    int nVar = vars->hess[iBlock].M();

    // Memory structure
    if( stats->itCount > gammaMat.N() )
    {
        m = gammaMat.N();
        posOldest = stats->itCount % m;
        posNewest = (stats->itCount-1) % m;
    }
    else
    {
        m = stats->itCount;
        posOldest = 0;
        posNewest = m-1;
    }

    // Construct \tilde W and D**(-1)
    W.Dimension( m, m ).Initialize( 0.0 );
    Dinv.Dimension( m ).Initialize( 0.0 );
    for( i=0; i<m; i++ )
    {
        // D**(-1) = diag( 1/deltaTdelta_i )
        posDelta = (posOldest+i) % m;
        deltai.Submatrix( deltaMat, nVar, 1, 0, posDelta );
        for( k=0; k<nVar; k++ )
            Dinv( i ) += deltai( k ) * deltai( k );
        Dinv( i ) = sqrt( Dinv( i ) );
        if( Dinv( i ) > param->eps )
            Dinv( i ) = 1.0 / Dinv( i );
        else
        {
            if( i == m-1 ) // don't issue warning every time for the same step
                printf("Warning! deltaTdelta (block %i, step %i) = %23.16e\n", iBlock, i, Dinv( i ));
            Dinv( i ) = 1.0e16;
            return 1; // choose different scaling factor
        }

        // \tilde W
        for( j=i; j<m; j++ )
        {
            posGamma = (posOldest+j) % m;
            gammai.Submatrix( gammaMat, nVar, 1, 0, posGamma );

            for( k=0; k<nVar; k++ )
                W( j, i ) += gammai( k ) * deltai( k );
        }
    }

    // D**(-1) \tilde W D**(-1)
    for( i=0; i<m; i++ )
        for( j=i; j<m; j++ )
            W( j, i ) = Dinv( j ) * W( j, i ) * Dinv( i );

    // Symmetry
    for( i=0; i<m; i++ )
        for( j=i+1; j<m; j++ )
            W( i, j ) = W( j, i );

    // Eigendecomposition
    int lwork = 3*m, info = 0;
    double *work;
    int ldim = W.LDIM();
    N.Dimension( m );
    C = Matrix( W );
    work = new double[lwork];

    dsyev_( "V", "L", &m, C.ARRAY(), &ldim,
            N.ARRAY(), work, &lwork, &info, strlen("V"), strlen("L") );

    if( info != 0 )
        printf( "Warning! DSYEV returned info = %i!\n", info );

    double minEv = 1.0e20;
    for( i=0; i<m; i++ )
        if( N( i ) < minEv ) minEv = N( i );
    //printf("minEv = %23.16e\n", minEv );
    if( minEv > 0.0 )
    {
        Matrix G, Q, NCD, NCDG, gammaT;
        G.Dimension( m, nVar ).Initialize( 0.0 );
        Q.Dimension( m, m ).Initialize( 0.0 );
        NCD.Dimension( m, m ).Initialize( 0.0 );
        NCDG.Dimension( m, m ).Initialize( 0.0 );
        gammaT.Dimension( m, nVar ).Initialize( 0.0 );

        // gammaT = gammaMat**T (reorder columns!)
        for( i=0; i<m; i++ )
            for( j=0; j<nVar; j++ )
            {
                posGamma = (posOldest+i) % m;
                gammaT( i, j ) = gammaMat( j, posGamma );
            }

        // N = N**(-1/2)
        for( i=0; i<m; i++ )
            N( i ) = 1.0 / sqrt( N( i ) );

        //
        // G = [N**(-1/2) * C**T * Dinv * gammaMat**T] * [gammaMat * Dinv * C * N**(-1/2)]
        //

        // NCD = N**(-1/2) * C**T * Dinv
        for( i=0; i<m; i++ )
            for( j=0; j<m; j++ )
                NCD( i, j ) = N( j ) * C( j, i ) * Dinv( i );

        // G = NCD * gammaT
        for( i=0; i<m; i++ )
            for( j=0; j<nVar; j++ )
                for( k=0; k<m; k++ )
                    G( i, j ) += NCD( i, k ) * gammaT( k, j );

        // Q = G * G**T
        for( i=0; i<m; i++ )
            for( j=0; j<m; j++ )
                for( k=0; k<nVar; k++ )
                    Q( i, j ) += G( i, k ) * G( j, k );

        Matrix ev;
        ev.Dimension( m ).Initialize( 0.0 );
        ldim = Q.LDIM();
        dsyev_( "N", "L", &m, Q.ARRAY(), &ldim,
               ev.ARRAY(), work, &lwork, &info, strlen("N"), strlen("L") );

        double maxEig = 0.0;
        for( i=0; i<m; i++ )
            if( ev( i ) > maxEig ) maxEig = ev( i );
        maxEig *= 1.5;

        for( i=0; i<vars->hess[iBlock].M(); i++ )
            for( j=i; j<vars->hess[iBlock].M(); j++ )
                vars->hess[iBlock]( i,j ) *= maxEig;

        printf("scale = %g\n", maxEig );
        stats->averageSizingFactor += maxEig;
    }
    else // choose different scaling factor
    {
        stats->averageSizingFactor += 1.0;
        return 1;
    }

    return 0;
}


void SQPmethod::calcHessianUpdateLimitedMemory( int updateType, int hessScaling, double mu )
{
    int iBlock, nBlocks, nVarLocal;
    Matrix smallGamma, smallDelta;
    Matrix gammai, deltai;
    int i, j, m, pos, posOldest, posNewest;
    int hessDamped, hessSkipped;
    double averageSizingFactor;

    //if objective derv is computed exactly, don't set the last block!
    if( param->whichSecondDerv == 1 && param->blockHess )
        nBlocks = vars->nBlocks - 1;
    else
        nBlocks = vars->nBlocks;

    // Statistics: how often is damping active, what is the average COL sizing factor?
    stats->hessDamped = 0;
    stats->hessSkipped = 0;
    stats->averageSizingFactor = 0.0;

    for( iBlock=0; iBlock<nBlocks; iBlock++ )
    {
        nVarLocal = vars->hess[iBlock].M();

        // smallGamma and smallDelta are submatrices of gammaMat, deltaMat,
        // i.e. subvectors of gamma and delta from m prev. iterations
        smallGamma.Submatrix( vars->gammaMat, nVarLocal, vars->gammaMat.N(), vars->blockIdx[iBlock], 0 );
        smallDelta.Submatrix( vars->deltaMat, nVarLocal, vars->deltaMat.N(), vars->blockIdx[iBlock], 0 );

        // Memory structure
        if( stats->itCount > smallGamma.N() )
        {
            m = smallGamma.N();
            posOldest = stats->itCount % m;
            posNewest = (stats->itCount-1) % m;
        }
        else
        {
            m = stats->itCount;
            posOldest = 0;
            posNewest = m-1;
        }

        // Set B_0 (pretend it's the first step)
        calcInitialHessian( iBlock );
        vars->deltaNormOld( iBlock ) = 1.0;
        vars->deltaGammaOld( iBlock ) = 0.0;
        vars->noUpdateCounter[iBlock] = -1;

        // Size the initial update, but with the most recent delta/gamma-pair
        gammai.Submatrix( smallGamma, nVarLocal, 1, 0, posNewest );
        deltai.Submatrix( smallDelta, nVarLocal, 1, 0, posNewest );
        if( hessScaling == 1 )
            sizeHessianNocedal( gammai, deltai, iBlock );
        else if( hessScaling == 2 )
            sizeHessianOL( gammai, deltai, iBlock );
        else if( hessScaling == 3 )
            sizeHessianMean( gammai, deltai, iBlock );
        else if( hessScaling == 5 )
            sizeHessianByrdLu( smallGamma, smallDelta, iBlock );

        for( i=0; i<m; i++ )
        {
            pos = (posOldest+i) % m;

            // Get new vector from list
            gammai.Submatrix( smallGamma, nVarLocal, 1, 0, pos );
            deltai.Submatrix( smallDelta, nVarLocal, 1, 0, pos );

            // Save statistics, we want to record them only for the most recent update
            averageSizingFactor = stats->averageSizingFactor;
            hessDamped = stats->hessDamped;
            hessSkipped = stats->hessSkipped;

            // Selective sizing before the update
            if( hessScaling == 4 )
                sizeHessianTapia( gammai, deltai, iBlock );

            // Compute the new update
            if( updateType == 1 && vars->updateSequence[pos] == 1 )
                calcSR1( gammai, deltai, iBlock );
            else if( updateType == 1 && vars->updateSequence[pos] == 2 )
                calcBFGS( gammai, deltai, iBlock );
            else if( updateType == 2 )
                calcBFGS( gammai, deltai, iBlock );

            // Count damping statistics only for the most recent update
            if( pos != posNewest )
            {
                stats->hessDamped = hessDamped;
                stats->hessSkipped = hessSkipped;
                if( hessScaling == 4 )
                    stats->averageSizingFactor = averageSizingFactor;
            }
        }

        // Convex combination between SR1 and scaled identity
        if( updateType == 1 && mu > 0.0 )
        {
            const double myEps = 1.0e2 * param->eps;
            const double fac1 = 1.0 - mu;
            double fac2;

            // H = (1-mu)*H_SR1 + mu*(scal*I)
            for( i=0; i<vars->hess[iBlock].M(); i++ )
                for( j=i; j<vars->hess[iBlock].M(); j++ )
                    vars->hess[iBlock]( i,j ) *= fac1;

            fac2 = adotb( deltai, gammai );
            if( fac2 > 0.0 )
                fac2 = adotb( gammai, gammai ) / fmax( fac2, myEps );
            else
                fac2 = 1.0;

            for( i=0; i<vars->hess[iBlock].M(); i++ )
                vars->hess[iBlock]( i,i ) += fac2;
        }

        // If an update is skipped to often, reset Hessian block
        if( vars->noUpdateCounter[iBlock] > param->maxConsecSkippedUpdates )
            resetHessian( iBlock );
    }//blocks
    stats->averageSizingFactor /= nBlocks;
}


void SQPmethod::calcBFGS( const Matrix &gamma, const Matrix &delta, int iBlock )
{
    int i, j, k, dim = gamma.M();
    Matrix Bdelta;
    SymMatrix *B;
    double h1 = 0.0;
    double h2 = 0.0;
    double thetaPowell = 0.0;
    int damped;

    // Work with a local copy of gamma because damping may need to change gamma.
    // Note that vars->gamma needs to remain unchanged! This may be important in a limited memory context:
    // When information is "forgotten", B_i-1 is different and the original gamma might lead to an undamped update with the new B_i-1!
    Matrix gamma2 = gamma;

    B = &vars->hess[iBlock];

    // Bdelta = B*delta (if sizing is enabled, B is the sized B!)
    // h1 = delta^T * B * delta
    // h2 = delta^T * gamma
    Bdelta.Dimension( dim ).Initialize( 0.0 );
    for( i=0; i<dim; i++ )
    {
        for( k=0; k<dim; k++ )
            Bdelta( i ) += (*B)( i,k ) * delta( k );

        h1 += delta( i ) * Bdelta( i );
        h2 += delta( i ) * gamma( i );
    }

    // Damping strategy to maintain pos. def. (Nocedal/Wright p.537; SNOPT paper)
    // Interpolates between current approximation and unmodified BFGS
    damped = 0;
    if( param->hessDamp )
        if( h2 < param->hessDampFac * h1 / vars->alpha && fabs( h1 - h2 ) > 1.0e-12 )
        {// At the first iteration h1 and h2 are equal due to Tapia scaling

            thetaPowell = (1.0 - param->hessDampFac)*h1 / ( h1 - h2 );

            // Redefine gamma and h2 = delta^T * gamma
            h2 = 0.0;
            for( i=0; i<dim; i++ )
            {
                gamma2( i ) = thetaPowell*gamma2( i ) + (1.0 - thetaPowell)*Bdelta( i );
                h2 += delta( i ) * gamma2( i );
            }

            // Also redefine deltaGammaOld for computation of sizing factor in the next iteration
            vars->deltaGammaOld( iBlock ) = h2;

            damped = 1;
        }

    // For statistics: count number of damped blocks
    stats->hessDamped += damped;

    // B_k+1 = B_k - Bdelta * (Bdelta)^T / h1 + gamma * gamma^T / h2
    double myEps = 1.0e1 * param->eps;
    if( fabs( h1 ) < myEps || fabs( h2 ) < myEps )
    {// don't perform update because of bad condition, might introduce negative eigenvalues
        /// \note it is not clear how to choose thresholds though
        //printf("block = %i, h1 = %g, h2 = %g\n", iBlock, h1, h2);
        vars->noUpdateCounter[iBlock]++;
        stats->hessDamped -= damped;
        stats->hessSkipped++;
    }
    else
    {
        for( i=0; i<dim; i++ )
            for( j=i; j<dim; j++ )
                (*B)( i,j ) = (*B)( i,j ) - Bdelta( i ) * Bdelta( j ) / h1
                                          + gamma2( i ) * gamma2( j ) / h2;

        vars->noUpdateCounter[iBlock] = 0;
    }
}


void SQPmethod::calcSR1( const Matrix &gamma, const Matrix &delta, int iBlock )
{
    int i, j, k, dim = gamma.M();
    Matrix gmBdelta;
    SymMatrix *B;
    double myEps = 1.0e1 * param->eps;
    double r = 1.0e-8;
    double h = 0.0;

    B = &vars->hess[iBlock];

    // gmBdelta = gamma - B*delta
    // h = (gamma - B*delta)^T * delta
    gmBdelta.Dimension( dim );
    for( i=0; i<dim; i++ )
    {
        gmBdelta( i ) = gamma( i );
        for( k=0; k<dim; k++ )
            gmBdelta( i ) -= ( (*B)( i,k ) * delta( k ) );

        h += ( gmBdelta( i ) * delta( i ) );
    }

    // B_k+1 = B_k + gmBdelta * gmBdelta^T / h
    if( fabs( h ) < r * l2VectorNorm( delta ) * l2VectorNorm( gmBdelta ) || fabs( h ) < myEps )
    {// Skip update if denominator is too small
        //printf("block %i, h = %23.16e\n", iBlock, h );
        vars->noUpdateCounter[iBlock]++;
        stats->hessSkipped++;
    }
    else
    {
        for( i=0; i<dim; i++ )
            for( j=i; j<dim; j++ )
                (*B)( i,j ) = (*B)( i,j ) + gmBdelta( i ) * gmBdelta( j ) / h;
        vars->noUpdateCounter[iBlock] = 0;
    }
}


/**
 * Set deltaXi and gamma as a column in the matrix containing
 * the m most recent delta and gamma
 */
void SQPmethod::updateDeltaGamma()
{
    int nVar = vars->gammaMat.M();
    int m = vars->gammaMat.N();

    if( m == 1 )
        return;

    vars->deltaXi.Submatrix( vars->deltaMat, nVar, 1, 0, stats->itCount % m );
    vars->gamma.Submatrix( vars->gammaMat, nVar, 1, 0, stats->itCount % m );
}

} // namespace blockSQP

