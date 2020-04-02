#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "fix_variable_kiwiel.h"
#define limit_for 10000



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

void imprime_resultado_x(cqk_problem *restrict p,double *x){
    double total = 0;
    for (int i = 0; i < p->n;i++){
     	total = total + x[i]*p->b[i];
        //  printf("x%d = %f  l = %f   u = %f  \n",i,x[i],p->low[i],p->up[i]);
    }
    printf("r %f\n",p->r);
    printf("total %f \n",total);
}

void calc_initial_lambda(cqk_problem *restrict p,double *lambda,double *qk,int *Ik){

    double soma_ab_by_d = 0;
    double soma_b_by_d = 0;
    // omp_set_num_threads(3);
    #pragma omp parallel
    {   
        #pragma omp for reduction(+:soma_ab_by_d,soma_b_by_d)
            for(unsigned i = 0; i < p->n;i++){
                soma_ab_by_d = soma_ab_by_d + ((p->a[i]*p->b[i])/p->d[i]);
                soma_b_by_d = soma_b_by_d + (pow(p->b[i],2)/p->d[i]);
            }
    }
    *qk = soma_b_by_d;
    *lambda = (soma_ab_by_d - p->r) / soma_b_by_d;
}

void checa_viabilidade(cqk_problem *restrict p,double qk,double *qk_l,double *qk_u,int *Ik,double *x,
    double lambda,double *nabla,double *delta){

    double soma_qk_l = 0;
    double soma_qk_u = 0;
    double tmp_nabla = 0;
    double tmp_delta = 0;
    

    #pragma omp parallel for reduction(+:soma_qk_l,soma_qk_u,tmp_nabla,tmp_delta)
    for(int i = 0; i < p->n;i++){
        if (Ik[i] != -1){
            
            // int id = omp_get_thread_num();

            if (Ik[i] == -2)
                Ik[i] = 1;
            else if (Ik[i] == -3)
                Ik[i] = 1;

            double xk = ((p->a[i] - lambda*p->b[i]) / p->d[i]);

            if (xk < p->low[i]){
                // printf("menor %d\n",i);
                x[i] = p->low[i];
                tmp_nabla += p->b[i] * (p->low[i] - xk);
                Ik[i] = -2;
                soma_qk_l += (pow(p->b[i],2)/p->d[i]);


            }else if (xk > p->up[i]){
                x[i] = p->up[i];
                tmp_delta += p->b[i] * (xk - p->up[i]);
                Ik[i] = -3;

                soma_qk_u += (pow(p->b[i],2)/p->d[i]);
            }else{
                x[i] = xk;
                // Ik[i] = -1;
            }
            
        }
        
    }
    // #pragma omp critical
    *qk_l = qk - soma_qk_l;
    *qk_u = qk - soma_qk_u;

    *nabla = tmp_nabla;
    *delta = tmp_delta;
          
}

void fixar_variaveis(cqk_problem *restrict p,double *qk, int *Ik, double *lambda,double qk_l, double qk_u,double nabla, double delta)
{
    if (nabla > delta){
            #pragma omp parallel for
            for (int i = 0; i < p->n;i++)
            {
                if (Ik[i] == -2 ){
                    Ik[i] = -1; 
                }
            }
            *lambda += (nabla / qk_l) ;
            *qk = qk_l;
        }else{
            #pragma omp parallel for
            for (int i = 0; i < p->n;i++)
            {
                if (Ik[i] == -3 ){
                    Ik[i] = -1;  
                }
            }
            
            *lambda += (-delta / qk_u) ;
            *qk = qk_u;
        }
}
int fix_variable_parallel(cqk_problem *restrict p,double *x){

    // Passo 0 //
    int *Ik = (int *) malloc(p->n*sizeof(int));
    double lambda = 0;
    double qk = 0;
    double qk_l = 0;
    double qk_u = 0;
    // FIM Passo 0 //

    calc_initial_lambda(p,&lambda,&qk,Ik);
    // Passo 2 //
    double nabla = 0;
    double delta = 0;
    
    int steps = 0;

    checa_viabilidade(p,qk,&qk_l,&qk_u,Ik,x,lambda,
            &nabla,&delta);
    // Fim passo 2 //
    while(nabla != delta){

        // Passo 4 //
        fixar_variaveis(p,&qk,Ik,&lambda,qk_l,qk_u,nabla,delta);
        // FIM Passo 4 //

        // Passo 2 //
        nabla = 0;
        delta = 0;
        checa_viabilidade(p,qk,&qk_l,&qk_u,Ik,x,lambda,
            &nabla,&delta);
        //FIM passo 2

        steps++;
        // printf("steps %d r = %f \n",steps,nabla);

    }
    // imprime_resultado_x(p,x);
    return steps;
    return 0;

}