#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "fix_variable_kiwiel.h"
#define NUM_THREADS 2

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

void calc_initial_lambda(cqk_problem *restrict p,double *lambda,double *qk,int *Ik,int *temp_Ikl, int *temp_Iku){

    double soma_ab_by_d = 0;
    double soma_b_by_d = 0;
    Ik[0] = p->n;
    #pragma omp parallel
    {   
        #pragma omp for reduction(+:soma_ab_by_d,soma_b_by_d)
            for(unsigned i = 0; i < p->n;i++){
                soma_ab_by_d = soma_ab_by_d + ((p->a[i]*p->b[i])/p->d[i]);
                soma_b_by_d = soma_b_by_d + (pow(p->b[i],2)/p->d[i]);
                Ik[i+1] = i;
            }
    }
    *qk = soma_b_by_d;
    *lambda = (soma_ab_by_d - p->r) / soma_b_by_d;
}

void checa_viabilidade(int steps,cqk_problem *restrict p,double qk,double *qk_l,double *qk_u,int *Ik,double *x,
    double lambda,int *temp_Ikl,int *temp_Iku,double *soma_delta_l,double *soma_delta_u, int *count_temp_Ikl, int *count_temp_Iku){

    double soma_qk_l = 0;
    double soma_qk_u = 0;
    for(unsigned j = 0; j < Ik[0];j++){
        // if (Ik[i] != -1){
            int i = Ik[j+1];
            double xk = ((p->a[i] - lambda*p->b[i]) / p->d[i]);
            if (xk < p->low[i]){
                // printf("menor %d\n",i);
                x[i] = p->low[i];
                temp_Ikl[*count_temp_Ikl] = i;
                *soma_delta_l += p->b[i] * (p->low[i] - xk);
                *count_temp_Ikl = *count_temp_Ikl+1;
                soma_qk_l += (pow(p->b[i],2)/p->d[i]);
                
            }else if (xk > p->up[i]){
                // printf("maior %d\n",i);
                x[i] = p->up[i];
                temp_Iku[*count_temp_Iku] = i;
                *soma_delta_u += p->b[i] * (xk - p->up[i]);
                *count_temp_Iku = *count_temp_Iku+1;
                soma_qk_u += (pow(p->b[i],2)/p->d[i]);
            }else{
                // printf("igual %d \n",i);
                x[i] = xk;
                // *qk_l = qk - (pow(p->b[i],2)/p->d[i]);
                // Ik[i] = -1;
            }
            
        // }
        
    }
    *qk_l = qk - soma_qk_l;
    *qk_u = qk - soma_qk_u;

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

int fix_variable_parallel(cqk_problem *restrict p,double *x){

    // Passo 0 //
    int *Ik = (int *) malloc(1+ p->n*sizeof(int));
    int *temp_Ikl = (int *) malloc(p->n*sizeof(int));
    int *temp_Iku = (int *) malloc(p->n*sizeof(int));
    double lambda = 0;
    double qk = 0;
    double qk_l = 0;
    double qk_u = 0;
    int tam_Ik = p->n;
    // FIM Passo 0 //

    calc_initial_lambda(p,&lambda,&qk,Ik,temp_Ikl,temp_Iku);
    // Passo 2 //
    double soma_delta_l = 0;
    double soma_delta_u = 0;
    int count_temp_Ikl = 0;
    int count_temp_Ukl = 0;
    
    int steps = 0;

    checa_viabilidade(steps,p,qk,&qk_l,&qk_u,Ik,x,lambda,temp_Ikl,temp_Iku,
            &soma_delta_l,&soma_delta_u,&count_temp_Ikl,&count_temp_Ukl);
    // Fim passo 2 //
    for (steps = 0; steps < 15;steps++){

        // printf("step %d \n", steps);

        // Passo 3 //
        if (soma_delta_l == soma_delta_u){
            // imprime_resultado_x(p,x);
            return steps;
        }
        // FIM Passo 3 //

        // Passo 4 //
        if (soma_delta_l > soma_delta_u){
            int count_temp = 0;
            int count_Ik = 0;
            for (unsigned i = 0; i < Ik[0];i++){
                if (Ik[i+1] != temp_Ikl[count_temp]){
                    Ik[count_Ik+1] = Ik[i+1];
                    count_Ik++;
                }else{
                    count_temp++;
                }
            }
            Ik[0] = count_Ik;

            lambda += (soma_delta_l / qk_l) ;
            qk = qk_l;

            // printf("count_temp %d \n",count_temp);
        }
        if (soma_delta_l < soma_delta_u){
            int count_temp = 0;
            int count_Ik = 0;
            for (unsigned i = 0; i < Ik[0];i++){
                if (Ik[i+1] != temp_Iku[count_temp]){
                    Ik[count_Ik+1] = Ik[i+1];
                    count_Ik++;
                }else{
                    count_temp++;
                }
            }
            Ik[0] = count_Ik;

            lambda += (-soma_delta_u / qk_u) ;
            qk = qk_u;
        }
        // FIM Passo 4 //

        // Passo 2 //
        soma_delta_l = 0;
        soma_delta_u = 0;
        count_temp_Ikl = 0;
        count_temp_Ukl = 0;
        checa_viabilidade(steps,p,qk,&qk_l,&qk_u,Ik,x,lambda,temp_Ikl,temp_Iku,
            &soma_delta_l,&soma_delta_u,&count_temp_Ikl,&count_temp_Ukl);
        //FIM passo 2
    }

    return 0;

}