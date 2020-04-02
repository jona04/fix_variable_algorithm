#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <omp.h>

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

double calc_lambda(cqk_problem *restrict p,int *Ik,int *Uk,int *Lk){
    double soma_bl;
    double soma_bu;
    double rk;
    double soma_ab_by_d;
    double soma_b_by_d;
    // printf("n = %d \n",omp_get_num_threads());
    int i;
    int j;

    #pragma omp parallel
    {   
        #pragma omp for reduction(+:soma_bl)
            for (i = 0;i < p->n; i++){
                if (Lk[i] != -1){
                    soma_bl = soma_bl + (p->low[i]*p->b[i]);
                }
            
            }

        #pragma omp for reduction(+:soma_bu)
            for (j = 0;j < p->n; j++){
                if (Uk[j] != -1){
                    soma_bu = soma_bu + (p->up[j]*p->b[j]);
                }
            }

        #pragma omp barrier
        rk = p->r - soma_bl - soma_bu;

        int k;
  
        #pragma omp for reduction(+:soma_ab_by_d,soma_b_by_d)
        for(k = 0; k < p->n;k++){
            if (Ik[k] != -1){
                // printf("Ik = %d\n",i);
                soma_ab_by_d = soma_ab_by_d + ((p->a[k]*p->b[k])/p->d[k]);
                soma_b_by_d = soma_b_by_d + (pow(p->b[k],2)/p->d[k]);
            }
        }
    
    }
    return (soma_ab_by_d - rk) / soma_b_by_d;
}

void checa_viabilidade(cqk_problem *restrict p,int *Ik,double *x,double lambda,int *temp_Ikl,int *temp_Iku,double *soma_delta_l,double *soma_delta_u, int *count_temp_Ikl, int *count_temp_Ukl){

    for(unsigned i = 0; i < p->n;i++){
        if (Ik[i] != -1){
            double xk = ((p->a[i] - lambda*p->b[i]) / p->d[i]);
            if (xk < p->low[i]){
                x[i] = p->low[i];
                temp_Ikl[*count_temp_Ikl] = i;
                *soma_delta_l = *soma_delta_l + p->b[i] * (p->low[i] - xk);
                *count_temp_Ikl = *count_temp_Ikl+1;

            }else if (xk > p->up[i]){
                x[i] = p->up[i];
                temp_Iku[*count_temp_Ukl] = i;
                *soma_delta_u = *soma_delta_u + p->b[i] * (xk - p->up[i]);
                *count_temp_Ukl = *count_temp_Ukl+1;
            }else{
                // printf("igual \n");
                x[i] = xk;
                // Ik[i] = -1;
            }
            
        }
        
    }
    
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
    int *Ik = (int *) malloc(p->n*sizeof(int));
    int *Uk = (int *) malloc(p->n*sizeof(int));
    int *Lk = (int *) malloc(p->n*sizeof(int));
    int *temp_Ikl = (int *) malloc(p->n*sizeof(int));
    int *temp_Iku = (int *) malloc(p->n*sizeof(int));

    double lambda = 0;
    
    for (unsigned i = 0; i<p->n;i++){
        Uk[i] = -1;
        Lk[i] = -1;
        temp_Ikl[i] = -1;
        temp_Iku[i] = -1;
    }
    int k = 1;

    // FIM Passo 0 //
    double total_passo1 = 0;
    for (unsigned steps = 0; steps < 20;steps++){

        printf("step %d \n", steps);
        
        // Passo 1 // 
        double start_omp = omp_get_wtime();
        
        lambda = calc_lambda(p,Ik,Uk,Lk);
        double end_omp = omp_get_wtime();
        total_passo1 += end_omp - start_omp;
        printf("Passo 1: %f \n", omp_get_wtime()-start_omp);
       
        // FIM Passo 1 //

        // Passo 2 //
        double soma_delta_l = 0;
        double soma_delta_u = 0;
        int count_temp_Ikl = 0;
        int count_temp_Ukl = 0;


        checa_viabilidade(p,Ik,x,lambda,temp_Ikl,temp_Iku,
            &soma_delta_l,&soma_delta_u,&count_temp_Ikl,&count_temp_Ukl);

        // FIM Passo 2 //

        // Passo 3 //
        if (soma_delta_l == soma_delta_u){
            imprime_resultado_x(p,x);
            printf("\npasso 1 %f \n", total_passo1);
            // printf("\npasso 2 %f \n", cpu_time_used2);
            // printf("\npasso 4 %f \n", cpu_time_used3);
            return steps;
        }
        // FIM Passo 3 //

        // Passo 4 //
        if (soma_delta_l > soma_delta_u){
            for (unsigned i = 0; i < count_temp_Ikl;i++){
                Ik[temp_Ikl[i]] = -1;
                Lk[temp_Ikl[i]] = 0;
                // printf("- %d / %d - ",temp_Ikl[i],Lk[temp_Ikl[i]]);
            }
        }
        if (soma_delta_l < soma_delta_u){
            for (unsigned i = 0; i < count_temp_Ukl;i++){
                Ik[temp_Iku[i]] = -1;
                Uk[temp_Iku[i]] = 0;
                // printf(" - %d - ",temp_Iku[i]);
            }
        }

        // FIM Passo 4 //

    }

    return 0;

}