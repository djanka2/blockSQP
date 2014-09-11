#include "blocksqp_problemspec.hpp"

namespace blockSQP
{

RestorationProblem::RestorationProblem( Problemspec *parentProblem, Matrix xiReference, Matrix constrReference )
{
    int i, iVar, iCon, ieqCnt;
    double myEps = 2.2204e-16;

    parent = parentProblem;

    // Need those to compute initial values and objective
    xiRef = Matrix( xiReference );
    constrRef = Matrix( constrReference );

    // Find out equality constraints
    neq = 0;
    isEqCon = new bool[parent->nCon];
    for( i=0; i<parent->nCon; i++ )
        if( fabs(parent->bu( parent->nVar+i ) - parent->bl( parent->nVar+i )) < myEps )
        {
            isEqCon[i] = 1;
            neq++;
        }
        else
            isEqCon[i] = 0;

    // We need nVar variables plus 2*nCon slack variables
    nVar = parent->nVar + 2*parent->nCon;
    // For each equality constraint, we need one constraint in the restoration problem,
    // for each inequality constraint we need two new constraints
    nCon = parent->nCon + ( parent->nCon - neq );
    // \todo nochmal genauer ueberlegen...
    nnCon = nCon;

    // If we want to use block matrices the variables must be ordered accordingly!
    nBlocks = 1;
    blockIdx = new int[2];
    blockIdx[0] = 0;
    blockIdx[1] = nVar;
    /*
    nBlocks = parent->nBlocks;
    blockIdx = new int[nBlocks+1];
    for( i=0; i<nBlocks+1; i++ )
        blockIdx[i] = parent->blockIdx[i];
    */

    /* Set bounds */

    // slack variables greater than zero
    bl.Dimension( nVar + nCon + 1 ).Initialize( 0.0 );
    bu.Dimension( nVar + nCon + 1 ).Initialize( myInf );

    objLo = 0.0;
    objUp = myInf;

    // bound constraints for original variables
    for( iVar=0; iVar<parent->nVar; iVar++ )
    {
        bl( iVar ) = parent->bl( iVar );
        bu( iVar ) = parent->bu( iVar );
    }

    ieqCnt = 0;
    for( iCon=0; iCon<parent->nCon; iCon++ )
    {
        if( isEqCon[iCon] )
        {// bounds for equality constraints stay the same
            bl( nVar + iCon ) = parent->bl( parent->nVar + iCon );
            bu( nVar + iCon ) = parent->bu( parent->nVar + iCon );
        }
        else
        {// inequalities in original problem are split up in two inequalities
            bl( nVar + iCon ) = parent->bl( parent->nVar + iCon );
            bu( nVar + iCon ) = myInf;

            bl( nVar + parent->nCon + ieqCnt ) = -myInf;
            bu( nVar + parent->nCon + ieqCnt ) = parent->bu( parent->nVar + iCon );

            ieqCnt++;
        }
    }
    if( ieqCnt != parent->nCon - neq ) printf("Error! ieqCnt = %i, parent->nCon - neq = %i\n",
                                              ieqCnt, parent->nCon - neq );
}


void RestorationProblem::evaluate( Matrix xi, Matrix lambda,
                                   double *objval, Matrix &constr,
                                   Matrix &gradObj, double *&jacNz, int *&jacIndRow, int *&jacIndCol,
                                   SymMatrix *&hess, int dmode, int *info )
{
    printf( "Implement sparse version of restoration phase!\n" );
    *info = 1;
}

void RestorationProblem::evaluate( Matrix xi, Matrix lambda,
                                   double *objval, Matrix &constr,
                                   Matrix &gradObj, Matrix &constrJac,
                                   SymMatrix *&hess, int dmode, int *info )
{
    int ieqCnt, iCon, iVar;
    double objOrig;
    Matrix xiOrig, constrJacOrig, constrOrig, gradObjOrig;
    Matrix pSlack, nSlack;
    Matrix lambdaDummy;
    SymMatrix *hessDummy;

    // The first nVar elements of the variable vector correspond to the variables of the original problem
    xiOrig.Submatrix( xi, parent->nVar, 1, 0, 0 );
    pSlack.Submatrix( xi, parent->nCon, 1, parent->nVar, 0 );
    nSlack.Submatrix( xi, parent->nCon, 1, parent->nVar+parent->nCon, 0 );
    gradObjOrig.Dimension( parent->nVar );

    // The first nCon elements of the constraint vector correspond to the constraints of the original problem
    constrOrig.Submatrix( constr, parent->nCon, 1, 0, 0 );
    constrJacOrig.Submatrix( constrJac, parent->nCon, parent->nVar, 0, 0 );

    // Evaluate constraints of the original problem
    parent->evaluate( xiOrig, lambdaDummy,
                      &objOrig, constrOrig,
                      gradObjOrig, constrJacOrig, hessDummy, dmode, info );

    /*
     * From this evaluation, build up constraints of the restoration problem
     * .-------------------------------------------------.
     * |                      {+ p(i) - n(i) <= up(i)    |
     * |   lo(i) <= constr(i)                            | nCon
     * |                      {+ p(i)        <= \infty   |
     * |-------------------------------------------------|
     * | -\infty <= constr(i)         - n(i) <= up(i)    | nCon - neq
     * '-------------------------------------------------'
     */
    ieqCnt = 0;
    for( iCon=0; iCon<parent->nCon; iCon++ )
    {
        if( isEqCon[iCon] )
        {// An equality constraint corresp to one equality constraints with pos and a neg slack
            constr( iCon ) = constrOrig( iCon ) + pSlack( iCon ) - nSlack( iCon );
            constrJac( iCon, parent->nVar + iCon ) = 1.0;
            constrJac( iCon, parent->nVar + parent->nCon + iCon ) = -1.0;
        }
        else
        {// An inequality constraint corresp to two inequality constraints, one with pos, one with neg slack

            // copy constraint value
            constr( parent->nCon + ieqCnt ) = constrOrig( iCon ) - nSlack( iCon );
            // copy derivatives w.r.t. original variables
            for( iVar=0; iVar<parent->nVar; iVar++ )
                constrJac( parent->nCon + ieqCnt, iVar ) = constrJacOrig( iCon, iVar );
            // derivative w.r.t. nSlack( iCon )
            constrJac( parent->nCon + ieqCnt, parent->nVar + parent->nCon + iCon ) = -1.0;

            // Be careful to first copy the original constraint value before adding the slack
            constr( iCon ) = constrOrig( iCon ) + pSlack( iCon );
            // derivative w.r.t. pSlack( iCon )
            constrJac( iCon, parent->nVar + iCon ) = 1.0;
            // derivatives w.r.t. original variables is already in place!

            ieqCnt++;
        }
    }
    if( ieqCnt != parent->nCon - neq ) printf("Error! ieqCnt = %i, parent->nCon - neq = %i\n",
                                              ieqCnt, parent->nCon - neq );

    evalObjective( xi, objval, gradObj, hess, dmode, info );
}


void RestorationProblem::evalObjective( Matrix xi, double *objval,
                                        Matrix &gradObj, SymMatrix *&hess,
                                        int dmode, int *info )
{
    if( dmode < 0 )
        return;

    int i;
    Matrix xiOrig, slack;

    // The first nVar elements of the variable vector correspond to the variables of the original problem
    xiOrig.Submatrix( xi, parent->nVar, 1, 0, 0 );
    slack.Submatrix( xi, 2*parent->nCon, 1, parent->nVar, 0 );

    *objval = 0.0;

    // First part: regularization term
    for( i=0; i<parent->nVar; i++ )
        *objval += diagScale( i ) * pow( xiOrig( i ) - xiRef( i ), 2);

    // Second part: sum of slack variables
    for( i=0; i<nVar - parent->nVar; i++ )
        *objval += rho * slack( i );

    if( dmode > 0 )
    {// compute gradient

        // gradient w.r.t. xi (regularization term)
        for( i=0; i<parent->nVar; i++ )
            gradObj( i ) = zeta * diagScale( i ) * diagScale( i ) * (xiOrig( i ) - xiRef( i ));
        // gradient w.r.t. slack variables
        for( i=parent->nVar; i<nVar; i++ )
            gradObj( i ) = rho;
    }

    if( dmode > 1 )
    {// compute hessian
        ; // n/a
    }

    *info = 0;
}


void RestorationProblem::convertJacobian( Matrix constrJac, double *&jacNz,
                                          int *&jacIndRow, int *&jacIndCol, bool firstCall )
{
    int nnz, i, j, count;

    //if( firstCall )
    if( true )
    {
        // 1st run: Count nonzeros
        nnz = 0;

        for( i=0; i<constrJac.M(); i++ )
            for( j=0; j<constrJac.N(); j++ )
                if( constrJac( i, j ) < myInf )
                    nnz++;

        if( jacNz != NULL ) delete[] jacNz;
        if( jacIndRow != NULL ) delete[] jacIndRow;

        jacNz = new double[nnz];
        jacIndRow = new int[nnz + (nVar+1) + nVar];
        jacIndCol = jacIndRow + nnz;
    }
    else
    {
        /* arrays jacInd* are already allocated! */
        nnz = jacIndCol[nVar];
    }

    // 2nd Run: store matrix entries columnwise in jacNz
    count = 0; // runs over all nonzero elements
    for( j=0; j<constrJac.N(); j++ )
        for( i=0; i<constrJac.M(); i++ )
        {
            jacIndCol[j] = count;
            if( fabs(constrJac(i,j)) < myInf )
            {
                jacNz[count] = constrJac(i,j);
                jacIndRow[count] = i;
                count++;
            }
        }
    jacIndCol[nVar] = count;
}


void RestorationProblem::initialize( Matrix &xi, Matrix &lambda, double *&jacNz, int *&jacIndRow, int *&jacIndCol )
{
    printf( "ToDo: Implement sparse version of restoration phase!\n" );
}


void RestorationProblem::initialize( Matrix &xi, Matrix &lambda, Matrix &constrJac )
{
    int i, iCon, iVar, ieqCnt;
    Matrix xiOrig, pSlack, nSlack, constrJacOrig;

    xiOrig.Submatrix( xi, parent->nVar, 1, 0, 0 );
    constrJacOrig.Submatrix( constrJac, parent->nCon, parent->nVar, 0, 0 );

    // Call setInitialValues of the parent problem to set up linear constraint matrix correctly
    parent->initialize( xiOrig, lambda, constrJacOrig );

    // Duplicate the linear INequality constraints
    ieqCnt = 0;
    for( iCon=0; iCon<parent->nCon; iCon++ )
        if( !isEqCon[iCon] )
        {
            for( iVar=0; iVar<parent->nVar; iVar++ )
                constrJac( parent->nCon + ieqCnt, iVar ) = constrJac( iCon, iVar );
            ieqCnt++;
        }

    // The reference point is the starting value for the restoration phase
    for( i=0; i<xiRef.M(); i++ )
        xiOrig( i ) = xiRef( i );

    // Initialize slack variables such that the constraints are feasible
    pSlack.Submatrix( xi, parent->nCon, 1, parent->nVar, 0 ).Initialize( 0.0 );
    nSlack.Submatrix( xi, parent->nCon, 1, parent->nVar+parent->nCon, 0 ).Initialize( 0.0 );
    for( i=0; i<parent->nCon; i++ )
    {
        // if lower bound is violated
        if( constrRef( i ) < parent->bl( parent->nVar + i ) )
            pSlack( i ) = parent->bl( parent->nVar + i ) - constrRef( i );
        else if( constrRef( i ) > parent->bu( parent->nVar + i ) )// if upper bound is violated
            nSlack( i ) = constrRef( i ) - parent->bu( parent->nVar + i );
    }

    // Set diagonal scaling matrix
    diagScale.Dimension( parent->nVar ).Initialize( 1.0 );
    for( i=0; i<parent->nVar; i++ )
        if( fabs( xi( i ) ) > 1.0 )
            diagScale( i ) = 1.0 / fabs( xiRef( i ) );

    // Regularization factor zeta and rho \todo wie setzen?
    zeta = 1.0e-3;
    rho = 1.0e3;

    lambda.Initialize( 0.0 );
}


void RestorationProblem::printVariables( Matrix xi, Matrix lambda, int verbose )
{
    int k;

    printf("\n<|----- Original Variables -----|>\n");
    for( k=0; k<parent->nVar; k++ )
        printf("%7i: %-30s   %7g <= %10.3g <= %7g   |   mul=%10.3g\n", k+1, parent->varNames[k], bl(k), xi(k), bu(k), lambda(k));
    printf("\n<|----- Slack Variables -----|>\n");
    for( k=parent->nVar; k<nVar; k++ )
        printf("%7i: slack   %7g <= %10.3g <= %7g   |   mul=%10.3g\n", k+1, bl(k), xi(k), bu(k), lambda(k));
}


void RestorationProblem::printConstraints( Matrix constr, Matrix lambda )
{

    /*
     * 31 red foreground
     * 32 green foreground
     * 33 brown foreground
     * 34 blue foreground
     * 35 magenta (purple) foreground
     * 36 cyan (light blue) foreground
     * 37 gray foreground
     */
    int k;
    double myeps = 1.0e-6;

    printf("\n<|----- Constraints -----|>\n");
    for( k=0; k<parent->nCon; k++ )
        printf("%5i: %-30s   %7g <= %10.4g <= %7g   |   mul=%10.3g\n", k+1, parent->conNames[parent->nVar+k], bl(nVar+k), constr(k), bu(nVar+k), lambda(nVar+k));
}


void RestorationProblem::printInfo()
{
    printf("Minimum norm NLP to find a point acceptable to the filter\n");
}

} // namespace blockSQP
