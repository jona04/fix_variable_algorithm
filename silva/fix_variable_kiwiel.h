/*
 * Module Name: cont_quad_knapsack
 *
 * Description: Implements the semismooth Newton method described in
 *     the manuscript 
 * 
 *     R. Cominetti, W. F. Mascarenhas, and Paulo J. S. Silva. "A
 *     Newton method for the continuous quadratic knapsack problem". 
 *
 *     It also implements a simple structure to hold the problems.
 *
 * Copyright: Paulo J. S. Silva <pjssilva@gmail.com> 2012.
 */

#ifndef CONT_QUAD_KNAPSACK_H
#define CONT_QUAD_KNAPSACK_H

#include <float.h>

/* 
 * Structure to describe a continuous quadratic knapsack problem in the form
 * 
 * min  1/2 x'Dx - a'x
 * s.t. b'x = r
 *      low <= x <= up,
 *
 * where D is a positive diagonal matrix.
 *
 * IMPORTANT: As in the manuscript the code assumes that b > 0. It is
 * the user responsibility to change the problem to fit this format.
 */
typedef struct cqk_problem {
    unsigned n;           /* Dimension of the problem. */
    double *restrict d;   /* D positive diagonal (D = diag(d)). */
    double *restrict a;   /* The a in the objective funtion definition. */
                          /* If the problem is interpreted as projection */
                          /* in the norm induced by D, the point being */  
                          /* projected is D^-1a. */
    double *restrict b;   /* Slopes (positive) that define the plane. */
    double r;             /* Righthand side of the plane equation. */
    double *restrict low; /* lower and */
    double *restrict up;  /* upper bounds for the constraints. */
} cqk_problem;

/********** Interface of a cqn_structure *********/

/*
 * Function Name: allocate_cqk_problem
 *
 * Description: Allocate the memory necessary hold a CQN problem in
 *     dimension n.
 *
 * Input:
 *     unsigned n: The dimension of the problem.
 * 
 * Output: cqk_problem *p: pointer to the cqk_problem whose memory is
 *     being allocated.
 */
void allocate_cqk_problem(unsigned n, cqk_problem *restrict p);

/*
 * Function Name: free_cqk_problem
 *
 * Description: Free the memory associated to a cqk_problem.
 *
 * Input/Output: 
 *     cqk_problem *p: at entry a cqk_problem with allocated
 *         memory, at output all the memory will be freed and the pointers
 *         in the cqk_problem structure will point to NULL.
 */
void free_cqk_problem(cqk_problem *restrict p);


double initial_multiplier(cqk_problem *restrict p, double *restrict x,
                          double *restrict slopes, double *sum_slopes,
                          unsigned *ind);
                          
/********** Interface for the Newton method **********/

/* Used to indicate an invalid value for a multiplier */
#define INVALIDLAMBDA DBL_MAX


/********** Algorithmic parameters **********/

/* Maximum number of iterations for each method. */
#define MAXITERS 100000

/* Precision required in the stopping criteria, it is set very high to
   ensure that the problems are solved "exactly". */
#define PREC 1.0e-12

/* If a bracket interval becomes very narrow, assume that its mid
   point is the solution. */
#define BRACKETPREC 2.2204460492503131e-16

/********** Interface to define new methods 

            If you are only willing to use the Newton solver you don't
            need to read or understand this.
 **********/

/********** Auxiliary macros **********/

/* Macro to print information while debugging. */
#ifdef DEBUG
#define DEBUG_PRINT(...) printf(__VA_ARGS__);
#else
#define DEBUG_PRINT(...) {};
#endif

/* Macro to compute the maximum and minimum between two values. */
#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

/* Structure to hold a bracket interval */
typedef struct {
    double neg_lambda; /* Largest multiplier to the left of the zero */
    double pos_lambda; /* Smallest multiplier to the right of the zero */
    double negPhi;     /* Value of phi at neg_lambda  */
    double posPhi;     /* Value of phi at pos_lambda    */
} bracket;

#endif
