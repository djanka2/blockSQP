/*
 * blockSQP -- Sequential quadratic programming for problems with
 *             block-diagonal Hessian matrix.
 * Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

/**
 * \file blocksqp_problemspec.hpp
 * \author Dennis Janka
 * \date 2012-2015
 *
 *  Declaration of ProblemSpec class to describes an NLP to be solved by blockSQP.
 */

#ifndef BLOCKSQP_PROBLEMSPEC_HPP
#define BLOCKSQP_PROBLEMSPEC_HPP

#include "blocksqp_defs.hpp"
#include "blocksqp_matrix.hpp"

namespace blockSQP
{

/**
 * \brief Base class for problem specification as required by SQPmethod.
 * \author Dennis Janka
 * \date 2012-2015
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

        double      objLo;              ///< lower bound for objective
        double      objUp;              ///< upper bound for objective
        Matrix      bl;                 ///< lower bounds of variables and constraints
        Matrix      bu;                 ///< upper bounds of variables and constraints

        int         nBlocks;            ///< number of separable blocks of Lagrangian
        int*        blockIdx;           ///< [blockwise] index in the variable vector where a block starts

    /*
     * METHODS
     */
    public:
        Problemspec( ){};
        virtual ~Problemspec( ){};

        /// Set initial values for xi and lambda, may also set matrix for linear constraints (dense version)
        virtual void initialize( Matrix &xi, Matrix &lambda, Matrix &constrJac ){};

        /// Set initial values for xi and lambda, may also set matrix for linear constraints (sparse version)
        virtual void initialize( Matrix &xi, Matrix &lambda, double *&jacNz, int *&jacIndRow, int *&jacIndCol ){};

        /// Evaluate all problem functions and their derivatives (dense version)
        virtual void evaluate( const Matrix &xi, const Matrix &lambda,
                               double *objval, Matrix &constr,
                               Matrix &gradObj, Matrix &constrJac,
                               SymMatrix *&hess, int dmode, int *info ){};

        /// Evaluate all problem functions and their derivatives (sparse version)
        virtual void evaluate( const Matrix &xi, const Matrix &lambda,
                               double *objval, Matrix &constr,
                               Matrix &gradObj, double *&jacNz, int *&jacIndRow, int *&jacIndCol,
                               SymMatrix *&hess, int dmode, int *info ){};

        /// Short cut if no derivatives are needed
        virtual void evaluate( const Matrix &xi, double *objval, Matrix &constr, int *info );

        /*
         * Optional Methods
         */
        /// Problem specific heuristic to reduce constraint violation
        virtual void reduceConstrVio( Matrix &xi, int *info ){};

        /// Print information about the current problem
        virtual void printInfo(){};
};

} // namespace blockSQP

#endif

