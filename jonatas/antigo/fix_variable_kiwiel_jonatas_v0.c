#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

void calc_initial_lambda(cqk_problem *restrict p,double *lambda,double *qk,double *pk,int *Ik,int *Uk,int *Lk,int *temp_Ikl, int *temp_Iku){

    // printf("step %d rk %f\n",steps,rk);
    double soma_ab_by_d = 0;
    double soma_b_by_d = 0;
    
    for(unsigned i = 0; i < p->n;i++){
        Uk[i] = -1;
        Lk[i] = -1;
        temp_Ikl[i] = -1;
        temp_Iku[i] = -1;
        if (Ik[i] != -1){
            // printf("Ik = %d\n",i);
            soma_ab_by_d = soma_ab_by_d + ((p->a[i]*p->b[i])/p->d[i]);
            soma_b_by_d = soma_b_by_d + (pow(p->b[i],2)/p->d[i]);
        }
    }
    *qk = soma_b_by_d;
    *pk = soma_ab_by_d;
    *lambda = (soma_ab_by_d - p->r) / soma_b_by_d;
}

void calc_lambda(cqk_problem *restrict p,double *lambda,int *Ik,int *Uk,int *Lk,int tam_lk,int tam_uk){
    double soma_bl = 0;
    for (unsigned i = 0;i < tam_lk; i++){
        if (Lk[i] != -1){
            soma_bl = soma_bl + (p->low[Lk[i]]*p->b[Lk[i]]);
        }
    }
    double soma_bu = 0;
    for (unsigned i = 0;i < tam_uk; i++){
        if (Uk[i] != -1){
            soma_bu = soma_bu + (p->up[Uk[i]]*p->b[Uk[i]]);
        }
    }
    double rk = p->r - soma_bl - soma_bu;
    // printf("step %d rk %f\n",steps,rk);
    double soma_ab_by_d = 0;
    double soma_b_by_d = 0;
    
    for(unsigned i = 0; i < p->n;i++){
        if (Ik[i] != -1){
            // printf("Ik = %d\n",i);
            soma_ab_by_d = soma_ab_by_d + ((p->a[i]*p->b[i])/p->d[i]);
            soma_b_by_d = soma_b_by_d + (pow(p->b[i],2)/p->d[i]);
        }
    }
    *lambda = (soma_ab_by_d - rk) / soma_b_by_d;
}

void checa_viabilidade(cqk_problem *restrict p,double qk,double *qk_l,double *qk_u,int *Ik,double *x,
    double lambda,int *temp_Ikl,int *temp_Iku,double *soma_delta_l,double *soma_delta_u, int *count_temp_Ikl, int *count_temp_Ukl){
    
    double soma_qk_l = 0;
    double soma_qk_u = 0;
    
    for(unsigned i = 0; i < p->n;i++){
        if (Ik[i] != -1){
            
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
                temp_Iku[*count_temp_Ukl] = i;
                *soma_delta_u += p->b[i] * (xk - p->up[i]);
                *count_temp_Ukl = *count_temp_Ukl+1;
                soma_qk_u += (pow(p->b[i],2)/p->d[i]);
            }else{
                // printf("igual %d \n",i);
                x[i] = xk;
                // *qk_l = qk - (pow(p->b[i],2)/p->d[i]);
                // Ik[i] = -1;
            }
            
        }
        
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

int fix_variable(cqk_problem *restrict p,double *x){



    // Passo 0 //
    int *Ik = (int *) malloc(p->n*sizeof(int));
    int *Uk = (int *) malloc(p->n*sizeof(int));
    int *Lk = (int *) malloc(p->n*sizeof(int));
    int *temp_Ikl = (int *) malloc(p->n*sizeof(int));
    int *temp_Iku = (int *) malloc(p->n*sizeof(int));
    int tam_lk = 0;
    int tam_uk = 0;
    double lambda = 0;
    double lambda_teste = 0;
    double qk = 0;
    double pk = 0;
    double rk = 0;
    double qk_l = 0;
    double qk_u = 0;

    int k = 1;

    // FIM Passo 0 //

    calc_initial_lambda(p,&lambda,&qk,&pk,Ik,Uk,Lk,temp_Ikl,temp_Iku);
    rk = p->r;
    // Passo 2 //
        double soma_delta_l = 0;
        double soma_delta_u = 0;
        int count_temp_Ikl = 0;
        int count_temp_Ukl = 0;

    checa_viabilidade(p,qk,&qk_l,&qk_u,Ik,x,lambda,temp_Ikl,temp_Iku,
            &soma_delta_l,&soma_delta_u,&count_temp_Ikl,&count_temp_Ukl);
    // printf("count_temp_Ikl = %d\n",count_temp_Ikl);
    // Fim passo 2 //
    for (unsigned steps = 0; steps < 10;steps++){

        // Passo 3 //
        if (soma_delta_l == soma_delta_u){
            // imprime_resultado_x(p,x);
            return steps;
        }
        // FIM Passo 3 //

        // Passo 4 //
        if (soma_delta_l > soma_delta_u){
            for (unsigned i = 0; i < count_temp_Ikl;i++){
                Ik[temp_Ikl[i]] = -1;
                Lk[tam_lk+i] = temp_Ikl[i];
            }
            tam_lk = tam_lk + count_temp_Ikl;
            lambda += (soma_delta_l / qk_l) ;
            qk = qk_l;
            // printf("fix = %d \n",count_temp_Ikl);
        }
        if (soma_delta_l < soma_delta_u){
            // printf("delta u \n");
            for (unsigned i = 0; i < count_temp_Ukl;i++){
                Ik[temp_Iku[i]] = -1;
                
                Uk[tam_uk+i] = temp_Iku[i];
            }
            tam_uk = tam_uk + count_temp_Ukl;
            lambda += (-soma_delta_u / qk_u) ;
            qk = qk_u;
            // printf("fix = %d \n",count_temp_Ukl);
        }
        // FIM Passo 4 //

        // Passo 1 // 
        // calc_lambda(p,&lambda,Ik,Uk,Lk,tam_lk,tam_uk);
        // printf("lambdas %f - %f \n",lambda,lambda_teste);
        // FIM Passo 1 //

        // Passo 2 //
        soma_delta_l = 0;
        soma_delta_u = 0;
        count_temp_Ikl = 0;
        count_temp_Ukl = 0;
        checa_viabilidade(p,qk,&qk_l,&qk_u,Ik,x,lambda,temp_Ikl,temp_Iku,
            &soma_delta_l,&soma_delta_u,&count_temp_Ikl,&count_temp_Ukl);
        //FIM passo 2
    }

    return 0;

}