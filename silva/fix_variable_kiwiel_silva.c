/*
 * Module Name: third_party_methods
 *
 * Description: Our implementation of established state-of-the-art
 *     methods to solve the continuous quadratic knapsack problem.
 *
 * Copyright: Paulo J. S. Silva <pjssilva@gmail.com> 2012.
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "fix_variable_kiwiel.h"


void allocate_cqk_problem(unsigned n, cqk_problem *restrict p) {
    p->n = n;
    p->d = (double *) malloc(p->n*sizeof(double));
    p->a = (double *) malloc(p->n*sizeof(double));
    p->b = (double *) malloc(p->n*sizeof(double));
    p->low = (double *) malloc(p->n*sizeof(double));
    p->up = (double *) malloc(p->n*sizeof(double));
    if (!p->d || !p->a || !p->b || !p->low || !p->up) {
        fprintf(stderr, "Memory allocation error, line %d, file %s\n",
                __LINE__, __FILE__);
        exit(1);
    }
}

void free_cqk_problem(cqk_problem *restrict p) {
    p->n = 0;
    free(p->d);
    p->d = NULL;
    free(p->a);
    p->a = NULL;
    free(p->b);
    p->b = NULL;
    free(p->low);
    p->low = NULL;
    free(p->up);
    p->up = NULL;
}

/********** Variable fixing method (variation of the Britan and Harx
            method) suggested in the paper

   Kiwiel KC. Variable Fixing Algorithms for the Continuous
   Quadratic Knapsack Problem. Journal of Optimization Theory and
   Applications. 2008;136(3):445-458. Available at:
   http://www.springerlink.com/index/10.1007/s10957-007-9317-7.

**********/ 





double initial_multiplier(cqk_problem *restrict p, double *restrict x,
                          double *restrict slopes, double *sum_slopes,
                          unsigned *ind) {

    double r = p->r;
    *sum_slopes = 0.0;
    double sum_cts = 0.0;
    double pre_slope;
    /* Estimate the initial multiplier from a given solution. */
    if (x != NULL) {
        for (unsigned i = 0; i < p->n; ++i) {
            ind[i] = i;
            pre_slope = p->b[i]/p->d[i];
            slopes[i] = pre_slope*p->b[i];
            if (x[i] == p->low[i] || x[i] == p->up[i])
                r -= p->b[i]*x[i];
            else {
                sum_cts += pre_slope*p->a[i];
                *sum_slopes += slopes[i];
            }
        }
        /* If all the variables were biding, a very rare and
           unexpected situation, throw the computed values away and
           try again now ignoring all the bounds. */
        if (*sum_slopes == 0.0) {
            r = p->r;
            for(unsigned i = 0; i < p->n; ++i) {
                pre_slope = p->b[i]/p->d[i];
                sum_cts += pre_slope*p->a[i];
                *sum_slopes += pre_slope*p->b[i];
            }
        }     
    }
    else {
        r = p->r;
        for(unsigned i = 0; i < p->n; ++i) {
            ind[i] = i;
            pre_slope = p->b[i]/p->d[i];
            sum_cts += pre_slope*p->a[i];
            slopes[i] = pre_slope*p->b[i];
            *sum_slopes += slopes[i];
        }
    }        

    return (r - sum_cts)/(*sum_slopes);
}







/*
 * Function Name: feasibility
 *
 * Description: Compute the feasibility measure associated to lambda,
 *    as described in Kiwiel's paper. It also finds the primal
 *    solution associated to the multiplier.
 *
 * Input:
 *    cqk_problem *restrict p: description of the cont. quad. knapsack.
 *    double slopes: slopes associated to each linear component of phi.
 *    double lambda: current multiplier.      
 *    double r: current rhs, after fixing variables.
 *    unsigned n: number of free variables.
 *    unsigned *restrict ind: indices of the free variables.
 *
 * Output:
 *    double *restrict nabla: negative component of the feasibility.
 *    double *restrict delta: positive component of the feasibility.
 *    double *restrict x: problem solution estimated from current 
 *       multiplier.
 */
void feasibility(cqk_problem *restrict p, double *restrict slopes, 
                 double lambda, unsigned n, unsigned *restrict ind, 
                 double *restrict x, double *restrict nabla, 
                 double *restrict delta) {

    *nabla = 0.0;
    *delta = 0.0;
    for (unsigned i = 0; i < n; ++i) {
        /* Actual index after indirection. */
        unsigned ii = ind[i];

        double new_x = (p->a[ii] + lambda*p->b[ii])/p->d[ii];
        if (new_x <= p->low[ii]) {
            *nabla += p->b[ii]*(p->low[ii] - new_x);
            new_x = p->low[ii];
        } else if (new_x >= p->up[ii]) {
            *delta += p->b[ii]*(new_x - p->up[ii]);
            new_x = p->up[ii];
        }
        x[ii] = new_x;
    }

    if (fabs(*nabla - *delta) <= PREC*(*nabla + *delta)) {
        *nabla = *delta;
    }
    DEBUG_PRINT("Feasibility = %.16e\n", *nabla - *delta);
}

/*
 * Function Name: kiwiel_fix
 *
 * Description: Fix variables whenever possible.
 *
 * Input:
 *    double feas: feasibility measure.
 *    cqk_problem *restrict p: description of the cont. quad. knapsack.
 *    double *restrict slopes: slopes associated to each linear component 
 *       of phi.
 *    double *restrict x: primal solution associated to current multiplier.
 *
 * Output:
 *    double *restrict r: updated rhs.
 *    unsigned *restrict n: number of free variables.
 *    unsigned *restrict ind: new indices of the free variables.
 *    double *restrict q: new sum of slopes.
 */
void kiwiel_fix(double feas, cqk_problem *restrict p, double *restrict slopes, 
                double *restrict x, double *restrict r, unsigned *restrict n, 
                unsigned *restrict ind, double *restrict q) {
    
    // printf("q %f \n",*q);
    unsigned len = 0;
    if (feas > 0)
        for (unsigned i = 0; i < *n; ++i) {
            /* Actual index after indirection. */
            unsigned ii = ind[i];

            if (x[ii] <= p->low[ii]) {
                *q -= slopes[ii];
                // printf("slopes[ii] %f \n",slopes[ii]);
                *r -= p->b[ii]*x[ii];
            }
            else {
                ind[len] = ii;
                ++len;
            }
        }
    else
        for (unsigned i = 0; i < *n; ++i) {
            /* Actual index after indirection. */
            unsigned ii = ind[i];

            if (x[ii] >= p->up[ii]) {
                *q -= slopes[ii];
                *r -= p->b[ii]*x[ii];
            }
            else {
                ind[len] = ii;
                ++len;
            }
        }
    DEBUG_PRINT("Fixed %d variables (from %d to %d).\n", 
                *n - len, *n, len);
    *n = len;
}

/*
 * Function Name: kiwiel_var_fix
 *
 * Description: Kiwiel's variation of the variable fixing method for
 *     the cont. quad. knapsack problem. The method is based on 
 *
 *     Kiwiel KC. Variable Fixing Algorithms for the Continuous
 *     Quadratic Knapsack Problem. Journal of Optimization Theory and
 *     Applications. 2008;136(3):445-458. Available at:
 *     http://www.springerlink.com/index/10.1007/s10957-007-9317-7.
 *
 * Input:
 *     cqk_problem *restrict p: the problem description.
 * 
 * Output:
 *     double *restrict x: solution vector.
 * 
 * Return value: the number of iterations required to compute the
 *     solution or -1 in case of failure.
 */
int kiwiel_var_fix(cqk_problem *restrict p, double *restrict x) {
    /* Allocate and initialize working area*/
    unsigned n = p->n;
    unsigned *restrict ind = (unsigned *) malloc(n*sizeof(unsigned));
    double *restrict slopes = (double *) malloc(n*sizeof(double));
    if (!ind || !slopes) {
        fprintf(stderr, "Memory allocation error, line %d, file %s\n",
                __LINE__, __FILE__);
        exit(1);
    }

    /* Initialization */
    unsigned n_iters = 0; /* Number of iterations */
    double 
        nabla = 0.0,      /* We want to solve the equation */
        delta = 1.0,      /* nabla = delta; */
        r = p->r,         /* Current rhs of the linear constraint */
        q,                /* Sum of the slopes of the free variables */
        lambda;           /* Multiplier with simple initial value */

    /* Estimate initial multiplier and initialize the free variable
       list. */
    DEBUG_PRINT("Starting Kiwiel variable fixing method\n");
    lambda = initial_multiplier(p, NULL, slopes, &q, ind);
    // printf("lambda %f \n",lambda);
    /* Start iteration */
    feasibility(p, slopes, lambda, n, ind, x, &nabla, &delta);
    n_iters++;
    DEBUG_PRINT("Initial lambda = %e, feasibility = %e\n",
                lambda, nabla - delta);
    while(nabla != delta) {
        if (n_iters > MAXITERS) break;
        
        kiwiel_fix(nabla - delta, p, slopes, x, &r, &n, ind, &q);
        if (nabla > delta)
            lambda -= nabla/q;
        else
            lambda += delta/q;

        // printf("lambda %f \n",lambda);
        // printf("nabla %f \n",nabla);
        // printf("q %f \n",q);

        feasibility(p, slopes, lambda, n, ind, x, &nabla, &delta);
        // printf("iter %d - lambda = %e, feasibility = %e\n", 
                    // n_iters, lambda, nabla - delta);
        DEBUG_PRINT("iter %d - lambda = %e, feasibility = %e\n", 
                    n_iters, lambda, nabla - delta);
        ++n_iters;
    }

    /* Free working memory */
    free(ind);
    free(slopes);

    /* Return the number of iterations */
    if (nabla == delta) {
        DEBUG_PRINT("Kiwiel projection done!\n\n");
        return n_iters;
    }
    else {
        DEBUG_PRINT("Kiwiel projection failed!\n\n");
        return -1;
    }
}