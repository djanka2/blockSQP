/*
 * blockSQP -- Sequential quadratic programming for problems with
 *             block-diagonal Hessian matrix.
 * Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef BLOCKSQP_HPP
#define BLOCKSQP_HPP

#include <omp.h>
#include "blocksqp_defs.hpp"
#include "blocksqp_matrix.hpp"
#include "blocksqp_problemspec.hpp"

namespace blockSQP
{

/**
 * \brief Contains algorithmic options and parameters for SQP method
 */
class SQPoptions
{
    /*
     * Variables
     */
    public:
        int printLevel;
        int printColor;
        int debugLevel;
        double eps;                         ///< values smaller than this are regarded as numerically zero
        double inf;                         ///< values larger than this are regarded as numerically infinity
        double opttol;                      ///< optimality tolerance
        double nlinfeastol;                 ///< nonlinear feasibility tolerance

        /* Algorithmic options */
        int sparseQP;                       ///< which qpOASES variant is used (dense/sparse/Schur)
        int globalization;                  ///< Globalization strategy
        int restoreFeas;                    ///< Use feasibility restoration phase
        int maxLineSearch;                  ///< Maximum number of steps in line search
        int maxConsecReducedSteps;          ///< Maximum number of consecutive reduced steps
        int maxConsecSkippedUpdates;        ///< Maximum number of consecutive skipped updates
        int maxItQP;                        ///< Maximum number of QP iterations per SQP iteration
        int blockHess;                      ///< Blockwise Hessian approximation?
        int hessScaling;                    ///< Scaling strategy for Hessian approximation
        int fallbackScaling;                ///< If indefinite update is used, the type of fallback strategy
        double maxTimeQP;                   ///< Maximum number of time in seconds per QP solve per SQP iteration
        double iniHessDiag;                 ///< Initial Hessian guess: diagonal matrix diag(iniHessDiag)
        double colEps;                      ///< epsilon for COL scaling strategy
        double colTau1;                     ///< tau1 for COL scaling strategy
        double colTau2;                     ///< tau2 for COL scaling strategy
        int hessDamp;                       ///< activate Powell damping for BFGS
        double hessDampFac;                 ///< damping factor for BFGS Powell modification
        int hessUpdate;                     ///< Type of Hessian approximation
        int fallbackUpdate;                 ///< If indefinite update is used, the type of fallback strategy
        int hessLimMem;                     ///< Full or limited memory
        int hessMemsize;                    ///< Memory size for L-BFGS updates
        int whichSecondDerv;                ///< For which block should second derivatives be provided by the user
        bool skipFirstGlobalization;        ///< If set to true, no globalization strategy in first iteration is applied
        int convStrategy;                   ///< Convexification strategy
        int maxConvQP;                      ///< How many additional QPs may be solved for convexification per iteration?

        /* Filter line search parameters */
        int maxSOCiter;                     ///< Maximum number of SOC line search iterations
        double gammaTheta;                  ///< see IPOPT paper
        double gammaF;                      ///< see IPOPT paper
        double kappaSOC;                    ///< see IPOPT paper
        double kappaF;                      ///< see IPOPT paper
        double thetaMax;                    ///< see IPOPT paper
        double thetaMin;                    ///< see IPOPT paper
        double delta;                       ///< see IPOPT paper
        double sTheta;                      ///< see IPOPT paper
        double sF;                          ///< see IPOPT paper
        double kappaMinus;                  ///< see IPOPT paper
        double kappaPlus;                   ///< see IPOPT paper
        double kappaPlusMax;                ///< see IPOPT paper
        double deltaH0;                     ///< see IPOPT paper
        double eta;                         ///< see IPOPT paper

    /*
     * Methods
     */
    public:
        SQPoptions();
        /// Some options cannot be used together. In this case set defaults
        void optionsConsistency();
};


/**
 * \brief Holds all variables that are updated during one SQP iteration
 */
class SQPiterate
{
    /*
     * Variables
     */
    public:
        double obj;                                   ///< objective value
        double qpObj;                                 ///< objective value of last QP subproblem
        double cNorm;                                 ///< constraint violation
        double cNormS;                                ///< scaled constraint violation
        double gradNorm;                              ///< norm of Lagrangian gradient
        double lambdaStepNorm;                        ///< norm of step in dual variables
        double tol;                                   ///< current optimality tolerance

        Matrix xi;                                    ///< variable vector
        Matrix lambda;                                ///< dual variables
        Matrix constr;                                ///< constraint vector

        Matrix constrJac;                             ///< full constraint Jacobian (not used in sparse mode)
        double *jacNz;                                ///< nonzero elements of Jacobian (length)
        int *jacIndRow;                               ///< row indices (length)
        int *jacIndCol;                               ///< indices to first entry of columns (nCols+1)

        Matrix deltaMat;                              ///< last m steps
        Matrix deltaXi;                               ///< alias for current step
        Matrix gradObj;                               ///< gradient of objective
        Matrix gradLagrange;                          ///< gradient of Lagrangian
        Matrix gammaMat;                              ///< Lagrangian gradient differences for last m steps
        Matrix gamma;                                 ///< alias for current Lagrangian gradient

        int nBlocks;                                  ///< number of diagonal blocks in Hessian
        int *blockIdx;                                ///< indices in the variable vector that correspond to diagonal blocks (nBlocks+1)

        SymMatrix *hess;                              ///< [blockwise] pointer to Hessian of the Lagrangian
        SymMatrix *hess1;                             ///< [blockwise] first Hessian approximation
        SymMatrix *hess2;                             ///< [blockwise] second Hessian approximation
        double *hessNz;                               ///< nonzero elements of Hessian (length)
        int *hessIndRow;                              ///< row indices (length)
        int *hessIndCol;                              ///< indices to first entry of columns (nCols+1)
        int *hessIndLo;                               ///< Indices to first entry of lower triangle (including diagonal) (nCols)

        /*
         * Variables for QP solver
         */
        Matrix deltaBl;                               ///< lower bounds for current step
        Matrix deltaBu;                               ///< upper bounds for current step
        Matrix lambdaQP;                              ///< dual variables of QP
        Matrix AdeltaXi;                              ///< product of constraint Jacobian with deltaXi

        /*
         * For modified BFGS updates
         */
        Matrix deltaNorm;                             ///< sTs
        Matrix deltaNormOld;                          ///< (from previous iteration)
        Matrix deltaGamma;                            ///< sTy
        Matrix deltaGammaOld;                         ///< (from previous iteration)
        int *noUpdateCounter;                         ///< count skipped updates for each block
        int *updateSequence;                          ///< if mixed updates are used in limited memory context

        /*
         * Variables for globalization strategy
         */
        int steptype;                                 ///< is current step a restoration step (1)?
        double alpha;                                 ///< stepsize for line search
        int nSOCS;                                    ///< number of second-order correction steps
        int reducedStepCount;                         ///< count number of consecutive reduced steps,
        Matrix deltaH;                                ///< scalars for inertia correction (filter line search w indef Hessian)
        Matrix trialXi;                               ///< new trial iterate (for line search)
        std::set< std::pair<double,double> > *filter; ///< Filter contains pairs (constrVio, objective)

    /*
     * Methods
     */
    public:
        /// Call allocation and initializing routines
        SQPiterate( Problemspec* prob, SQPoptions* param, bool full );
        SQPiterate( const SQPiterate &iter );
        /// Allocate variables that any SQP code needs
        void allocMin( Problemspec* prob );
        /// Allocate space for Jacobian
        void allocJac( Problemspec *prob );
        /// Allocate diagonal block Hessian
        void allocHess( SQPoptions* param );
        /// Convert *hess to column compressed sparse format
        void convertHessian( Problemspec *prob, double eps, SymMatrix *&hess_,
                             double *&hessNz_, int *&hessIndRow_, int *&hessIndCol_, int *&hessIndLo_ );
        /// Convert *hess to double array (dense matrix)
        void convertHessian( Problemspec *prob, double eps, SymMatrix *&hess_ );
        /// Allocate variables specifically needed by vmused SQP method
        void allocAlg( Problemspec* prob, SQPoptions* param );
        /// Set initial filter, objective function, tolerances etc.
        void initIterate( SQPoptions* param );
        ~SQPiterate( void );
};


/**
 * \brief Contains objects and pointers to files that record statistics during an SQP run
 */
class SQPstats
{
    /*
     * Variables
     */
    public:
        int itCount;                 ///< iteration number
        int qpIterations;            ///< number of qp iterations in the current major iteration
        int qpIterations2;           ///< number of qp iterations for solving convexified QPs
        int qpItTotal;               ///< total number of qp iterations
        int qpResolve;               ///< how often has QP to be convexified and resolved?
        int nFunCalls;               ///< number of function calls
        int nDerCalls;               ///< number of derivative calls
        int nRestHeurCalls;          ///< number calls to feasibility restoration heuristic
        int nRestPhaseCalls;         ///< number calls to feasibility restoration phase
        int rejectedSR1;             ///< count how often the SR1 update is rejected
        int hessSkipped;             ///< number of block updates skipped in the current iteration
        int hessDamped;              ///< number of block updates damped in the current iteration
        double averageSizingFactor;  ///< average value (over all blocks) of COL sizing factor
        PATHSTR outpath;             ///< path where log files are stored

        FILE *progressFile;          ///< save stats for each SQP step
        FILE *updateFile;            ///< print update sequence (SR1/BFGS) to file
        FILE *primalVarsFile;        ///< primal variables for every SQP iteration
        FILE *dualVarsFile;          ///< dual variables for every SQP iteration
        FILE *jacFile;               ///< Jacobian of one iteration
        FILE *hessFile;              ///< Hessian of one iteration

    /*
     * Methods
     */
    public:
        SQPstats( PATHSTR myOutpath );
        /// Open output files
        void initStats( SQPoptions *param );
        /// Print Debug information in logfiles
        void printDebug( SQPiterate *vars, SQPoptions *param );
        /// Print current iterate of primal variables to file
        void printPrimalVars( const Matrix &xi );
        /// Print current iterate of dual variables to file
        void printDualVars( const Matrix &lambda );
        /// Print all QP data to files to be read in MATLAB
        void dumpQPMatlab( Problemspec *prob, SQPiterate *vars, int sparseQP );
        void dumpQPCpp( Problemspec *prob, SQPiterate *vars, qpOASES::SQProblem *qp, int sparseQP );
        void printVectorCpp( FILE *outfile, double *vec, int len, char* varname );
        void printVectorCpp( FILE *outfile, int *vec, int len, char* varname );
        void printCppNull( FILE *outfile, char* varname );
        /// Print current (full) Jacobian to Matlab file
        void printJacobian( const Matrix &constrJacFull );
        void printJacobian( int nCon, int nVar, double *jacNz, int *jacIndRow, int *jacIndCol );
        /// Print current (full) Hessian to Matlab file
        void printHessian( int nBlocks, SymMatrix *&hess );
        void printHessian( int nVar, double *hesNz, int *hesIndRow, int *hesIndCol );
        /// Print a sparse Matrix in (column compressed) to a MATLAB readable file
        void printSparseMatlab( FILE *file, int nRow, int nVar, double *nz, int *indRow, int *indCol );
        /// Print one line of output to stdout about the current iteration
        void printProgress( Problemspec *prob, SQPiterate *vars, SQPoptions *param, bool hasConverged );
        /// Must be called before returning from run()
        void finish( SQPoptions *param );
};



/**
 * \brief SQP method for a given problem and set of options
 */
class SQPmethod
{
    /*
     * Variables
     */
    public:
        Problemspec*             prob;        ///< Problem structure (has to provide evaluation routines)
        SQPiterate*              vars;        ///< All SQP variables for this method
        SQPoptions*              param;       ///< Set of algorithmic options and parameters for this method
        SQPstats*                stats;       ///< Statistics object for current SQP run
        qpOASES::SQProblem*      qp;          ///< qpOASES qp object
        qpOASES::SQProblem*      qpSave;      ///< qpOASES qp object

    private:
        bool                     initCalled;  ///< indicates if init() has been called (necessary for run())

    /*
     * Methods
     */
    public:
        /// Construct a method for a given problem and set of algorithmic options
        SQPmethod( Problemspec *problem, SQPoptions *parameters, SQPstats *statistics );
        ~SQPmethod();
        /// Initialization, has to be called before run
        void init();
        /// Main Loop of SQP method
        int run( int maxIt, int warmStart = 0 );
        /// Call after the last call of run, to close output files etc.
        void finish();
        /// Print information about the SQP method
        void printInfo( int printLevel );
        /// Compute gradient of Lagrangian function (dense version)
        void calcLagrangeGradient( const Matrix &lambda, const Matrix &gradObj, const Matrix &constrJacFull, Matrix &gradLagrange, int flag );
        /// Compute gradient of Lagrangian function (sparse version)
        void calcLagrangeGradient( const Matrix &lambda, const Matrix &gradObj, double *jacNz, int *jacIndRow, int *jacIndCol, Matrix &gradLagrange, int flag );
        /// Overloaded function for convenience, uses current variables of SQPiterate vars
        void calcLagrangeGradient( Matrix &gradLagrange, int flag );
        /// Update optimization tolerance (similar to SNOPT) in current iterate
        bool calcOptTol();

        /*
         * Solve QP subproblem
         */
        /// Update the bounds on the current step, i.e. the QP variables
        void updateStepBounds( bool soc );
        /// Solve a QP with QPOPT or qpOASES to obtain a step deltaXi and estimates for the Lagrange multipliers
        int solveQP( Matrix &deltaXi, Matrix &lambdaQP, bool matricesChanged = true );
        void computeNextHessian( int idx, int maxQP );

        /*
         * Globalization Strategy
         */
        /// No globalization strategy
        int fullstep();
        /// Set new primal dual iterate
        void acceptStep( const Matrix &deltaXi, const Matrix &lambdaQP, double alpha, int nSOCS );
        /// Overloaded function for convenience, uses current variables of SQPiterate vars
        void acceptStep( double alpha );
        /// Reduce stepsize if a step is rejected
        void reduceStepsize( double *alpha );
        /// Determine steplength alpha by a filter based line search similar to IPOPT
        int filterLineSearch();
        /// Remove all entries from filter
        void initializeFilter();
        /// Is a pair (cNorm, obj) in the current filter?
        bool pairInFilter( double cNorm, double obj );
        /// Augment current filter by pair (cNorm, obj)
        void augmentFilter( double cNorm, double obj );
        /// Perform a second order correction step (solve QP)
        bool secondOrderCorrection( double cNorm, double cNormTrial, double dfTdeltaXi, bool swCond, int it );
        /// Reduce stepsize if a second order correction step is rejected
        void reduceSOCStepsize( double *alphaSOC );
        /// Start feasibility restoration heuristic
        int feasibilityRestorationHeuristic();
        /// Start feasibility restoration phase (solve NLP)
        int feasibilityRestorationPhase();
        /// Check if full step reduces KKT error
        int kktErrorReduction( );

        /*
         * Hessian Approximation
         */
        /// Set initial Hessian: Identity matrix
        void calcInitialHessian();
        /// [blockwise] Set initial Hessian: Identity matrix
        void calcInitialHessian( int iBlock );
        /// Reset Hessian to identity and remove past information on Lagrange gradient and steps
        void resetHessian();
        /// [blockwise] Reset Hessian to identity and remove past information on Lagrange gradient and steps
        void resetHessian( int iBlock );
        /// Compute current Hessian approximation by finite differences
        int calcFiniteDiffHessian();
        /// Compute full memory Hessian approximations based on update formulas
        void calcHessianUpdate( int updateType, int hessScaling );
        /// Compute limited memory Hessian approximations based on update formulas
        void calcHessianUpdateLimitedMemory( int updateType, int hessScaling );
        /// [blockwise] Compute new approximation for Hessian by SR1 update
        void calcSR1( const Matrix &gamma, const Matrix &delta, int iBlock );
        /// [blockwise] Compute new approximation for Hessian by BFGS update with Powell modification
        void calcBFGS( const Matrix &gamma, const Matrix &delta, int iBlock );
        /// Set pointer to correct step and Lagrange gradient difference in a limited memory context
        void updateDeltaGamma();

        /*
         * Scaling of Hessian Approximation
         */
        /// [blockwise] Update scalars for COL sizing of Hessian approximation
        void updateScalars( const Matrix &gamma, const Matrix &delta, int iBlock );
        /// [blockwise] Size Hessian using SP, OL, or mean sizing factor
        void sizeInitialHessian( const Matrix &gamma, const Matrix &delta, int iBlock, int option );
        /// [blockwise] Size Hessian using the COL scaling factor
        void sizeHessianCOL( const Matrix &gamma, const Matrix &delta, int iBlock );
};

} // namespace blockSQP

#endif
