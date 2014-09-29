/**
 * \file blocksqp_problemspec.hpp
 * \author Dennis Janka
 */

#ifndef BLOCKSQP_PROBLEMSPEC_HPP
#define BLOCKSQP_PROBLEMSPEC_HPP

#include "blocksqp_matrix.hpp"

namespace blockSQP
{

/**
 * \brief Base class for problem specification as required by SQP algorithm
 */
class Problemspec
{
    /*
     * CLASS VARIABLES
     */
    public:
        int         nVar;               ///< number of variables
        int         nCon;               ///< number of constraints
        int         nnCon;              ///< number of nonlinear constraints
        char**      varNames;           ///< names of variables
        char**      conNames;           ///< names of constraints

        double      objLo;              ///< lower bound for objective
        double      objUp;              ///< upper bound for objective
        Matrix      bl;                 ///< lower bounds of variables and constraints
        Matrix      bu;                 ///< upper bounds of variables and constraints

        int         nBlocks;            ///< number of separable blocks of Lagrangian
        int*        blockIdx;           ///< [blockwise] index in the variable vector where a block starts

    protected:
        double      objScale;           ///< scaling factor of objective function
        int         nFunCalls;          ///< number of function calls
        int         nDerCalls;          ///< number of derivative calls

    /*
     * METHODS
     */
    public:
        Problemspec( ) : objScale(0.0), nFunCalls(0), nDerCalls(0){};

        /// Set initial values for xi and lambda, may also set matrix for linear constraints (dense version)
        virtual void initialize( Matrix &xi, Matrix &lambda, Matrix &constrJac ) = 0;

        /// Set initial values for xi and lambda, may also set matrix for linear constraints (sparse version)
        virtual void initialize( Matrix &xi, Matrix &lambda, double *&jacNz, int *&jacIndRow, int *&jacIndCol ) = 0;

        /// Evaluate all problem functions and their derivatives (dense version)
        virtual void evaluate( const Matrix &xi, const Matrix &lambda,
                               double *objval, Matrix &constr,
                               Matrix &gradObj, Matrix &constrJac,
                               SymMatrix *&hess, int dmode, int *info ) = 0;

        /// Evaluate all problem functions and their derivatives (sparse version)
        virtual void evaluate( const Matrix &xi, const Matrix &lambda,
                               double *objval, Matrix &constr,
                               Matrix &gradObj, double *&jacNz, int *&jacIndRow, int *&jacIndCol,
                               SymMatrix *&hess, int dmode, int *info ) = 0;

        /// Short cut if no derivatives are needed
        virtual void evaluate( const Matrix &xi, double *objval, Matrix &constr, int *info );

        /*
         * Optional Methods
         */
        /// Problem specific heuristic to reduce constraint violation
        virtual void reduceConstrVio( Matrix &xi, int *info ){};

        /// Print information about the current problem
        virtual void printInfo(){};
        virtual void printVariables( const Matrix &xi, const Matrix &lambda, int verbose ){};
        virtual void printConstraints( const Matrix &constr, const Matrix &lambda ){};
};


/**
 * \brief If feasibility restoration phase is invoked, create an NLP to minimize constraint violation
 */
class RestorationProblem : public Problemspec
{
    /*
     * CLASS VARIABLES
     */
    public:
        Problemspec *parent;
        Matrix xiRef;
        Matrix constrRef;
        Matrix diagScale;
        int neq;
        bool *isEqCon;

        double zeta;
        double rho;

    /*
     * METHODS
     */
    public:
        RestorationProblem( Problemspec *parent, const Matrix &xiReference, const Matrix &constrReference );

        /// Set initial values for xi and lambda, may also set matrix for linear constraints (dense version)
        virtual void initialize( Matrix &xi, Matrix &lambda, Matrix &constrJac );

        /// Set initial values for xi and lambda, may also set matrix for linear constraints (sparse version)
        virtual void initialize( Matrix &xi, Matrix &lambda, double *&jacNz, int *&jacIndRow, int *&jacIndCol );

        /// Evaluate all problem functions and their derivatives (dense version)
        virtual void evaluate( const Matrix &xi, const Matrix &lambda,
                               double *objval, Matrix &constr,
                               Matrix &gradObj, Matrix &constrJac,
                               SymMatrix *&hess, int dmode, int *info );

        /// Evaluate all problem functions and their derivatives (sparse version)
        virtual void evaluate( const Matrix &xi, const Matrix &lambda,
                               double *objval, Matrix &constr,
                               Matrix &gradObj, double *&jacNz, int *&jacIndRow, int *&jacIndCol,
                               SymMatrix *&hess, int dmode, int *info );

        virtual void convertJacobian( const Matrix &constrJac, double *&jacNz, int *&jacIndRow,
                                      int *&jacIndCol, bool firstCall = 0 );

        virtual void evalObjective( const Matrix &xi, double *objval,
                                    Matrix &gradObj, SymMatrix *&hess,
                                    int dmode, int *info );

        virtual void printInfo();
        virtual void printVariables( const Matrix &xi, const Matrix &lambda, int verbose );
        virtual void printConstraints( const Matrix &constr, const Matrix &lambda );
};

} // namespace blockSQP

#endif

