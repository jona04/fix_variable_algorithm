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

void calc_initial_lambda(cqk_problem *restrict p,double *lambda,double *qk,int *Ik,int *temp_Ikl, int *temp_Iku){

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
    double lambda,int *temp_Ikl,int *temp_Iku,double *soma_delta_l,double *soma_delta_u, int *count_temp_Ikl, int *count_temp_Iku){

    double soma_qk_l = 0;
    double soma_qk_u = 0;
    double tmp_soma_delta_l = 0;
    double tmp_soma_delta_u = 0;
    
    int total_ikl = 0;
    int total_iku = 0;

    //n thread 4
    static int temp_Ikl_1[10000000];
     int tmp_count_temp_Ikl_1 = 0;
    static int temp_Ikl_2[10000000];
     int tmp_count_temp_Ikl_2 = 0;
    static int temp_Ikl_3[10000000];
     int tmp_count_temp_Ikl_3 = 0;
    static int temp_Ikl_4[10000000];
     int tmp_count_temp_Ikl_4 = 0;

    static int temp_Iku_1[10000000];
     int tmp_count_temp_Iku_1 = 0;
    static int temp_Iku_2[10000000];
     int tmp_count_temp_Iku_2 = 0;
    static int temp_Iku_3[10000000];
     int tmp_count_temp_Iku_3 = 0;
    static int temp_Iku_4[10000000];
    int tmp_count_temp_Iku_4 = 0;

    #pragma omp parallel for reduction(+:soma_qk_l,soma_qk_u,tmp_soma_delta_l,tmp_soma_delta_u)
    for(int i = 0; i < p->n;i++){
        if (Ik[i] != -1){
            
            int id = omp_get_thread_num();

            double xk = ((p->a[i] - lambda*p->b[i]) / p->d[i]);

            if (xk < p->low[i]){
                // printf("menor %d\n",i);
                x[i] = p->low[i];
                tmp_soma_delta_l += p->b[i] * (p->low[i] - xk);

                if (id == 0)
                {
                    temp_Ikl_1[tmp_count_temp_Ikl_1] = i;
                    tmp_count_temp_Ikl_1++;
                }else if (id == 1)
                {
                    temp_Ikl_2[tmp_count_temp_Ikl_2] = i;
                    tmp_count_temp_Ikl_2++;
                }else if (id == 2)
                {
                    temp_Ikl_3[tmp_count_temp_Ikl_3] = i;
                    tmp_count_temp_Ikl_3++;
                }else
                {
                    temp_Ikl_4[tmp_count_temp_Ikl_4] = i;
                    tmp_count_temp_Ikl_4++;
                }

                soma_qk_l += (pow(p->b[i],2)/p->d[i]);

            }else if (xk > p->up[i]){
                x[i] = p->up[i];
                tmp_soma_delta_u += p->b[i] * (xk - p->up[i]);

                if (id == 0)
                {
                    temp_Iku_1[tmp_count_temp_Iku_1] = i;
                    tmp_count_temp_Iku_1++;
                }else if (id == 1)
                {
                    temp_Iku_2[tmp_count_temp_Iku_2] = i;
                    tmp_count_temp_Iku_2++;
                }else if (id == 2)
                {
                    temp_Iku_3[tmp_count_temp_Iku_3] = i;
                    tmp_count_temp_Iku_3++;
                }else
                {
                    temp_Iku_4[tmp_count_temp_Iku_4] = i;
                    tmp_count_temp_Iku_4++;
                }

                soma_qk_u += (pow(p->b[i],2)/p->d[i]);

            }else{
                x[i] = xk;
            }
            
        }
        
    }
    // #pragma omp critical
    *qk_l = qk - soma_qk_l;
    *qk_u = qk - soma_qk_u;

    *soma_delta_u = tmp_soma_delta_u;
    *soma_delta_l = tmp_soma_delta_l;

    total_ikl = tmp_count_temp_Ikl_1+tmp_count_temp_Ikl_2+tmp_count_temp_Ikl_3+tmp_count_temp_Ikl_4;
    total_iku = tmp_count_temp_Iku_1+tmp_count_temp_Iku_2+tmp_count_temp_Iku_3+tmp_count_temp_Iku_4;
    
    // printf("total_ikl %d\n",total_ikl);
    // printf("total_iku %d\n",total_iku);
    // if(total_ikl > 10000){
    // #pragma omp parallel 
    // {
        // ------ Ikl --------
        #pragma omp parallel for if (tmp_count_temp_Ikl_1 > limit_for)
        for(int i = 0; i<tmp_count_temp_Ikl_1;i++){
            temp_Ikl[i] = temp_Ikl_1[i];
        }
        #pragma omp parallel for if (tmp_count_temp_Ikl_2 > limit_for)
        for(int i = 0; i<tmp_count_temp_Ikl_2;i++){
            temp_Ikl[i+tmp_count_temp_Ikl_1] = temp_Ikl_2[i];
        }
        #pragma omp parallel for if (tmp_count_temp_Ikl_3 > limit_for)
        for(int i = 0; i<tmp_count_temp_Ikl_3;i++){
            temp_Ikl[i+tmp_count_temp_Ikl_1+tmp_count_temp_Ikl_2] = temp_Ikl_3[i];
        }
        #pragma omp parallel for if (tmp_count_temp_Ikl_4 > limit_for)
        for(int i = 0; i<tmp_count_temp_Ikl_4;i++){
            temp_Ikl[i+tmp_count_temp_Ikl_1+tmp_count_temp_Ikl_2+tmp_count_temp_Ikl_3] = temp_Ikl_4[i];
        }
    

        // ----------- Iku ------------
        #pragma omp parallel for if (tmp_count_temp_Iku_1 > limit_for)
        for(int i = 0; i<tmp_count_temp_Iku_1;i++){
            temp_Iku[i] = temp_Iku_1[i];
        }
        #pragma omp parallel for if (tmp_count_temp_Iku_2 > limit_for)
        for(int i = 0; i<tmp_count_temp_Iku_2;i++){
            temp_Iku[i+tmp_count_temp_Iku_1] = temp_Iku_2[i];
        }
        #pragma omp parallel for if (tmp_count_temp_Iku_3 > limit_for)
        for(int i = 0; i<tmp_count_temp_Iku_3;i++){
            temp_Iku[i+tmp_count_temp_Iku_1+tmp_count_temp_Iku_2] = temp_Iku_3[i];
        }
        #pragma omp parallel for if (tmp_count_temp_Iku_4 > limit_for)
        for(int i = 0; i<tmp_count_temp_Iku_4;i++){
            temp_Iku[i+tmp_count_temp_Iku_1+tmp_count_temp_Iku_2+tmp_count_temp_Iku_3] = temp_Iku_4[i];
        }
    // }
          
    *count_temp_Ikl = total_ikl;
    *count_temp_Iku = total_iku;

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
    int *temp_Ikl = (int *) malloc(p->n*sizeof(int));
    int *temp_Iku = (int *) malloc(p->n*sizeof(int));
    double lambda = 0;
    double qk = 0;
    double qk_l = 0;
    double qk_u = 0;
    // FIM Passo 0 //

    calc_initial_lambda(p,&lambda,&qk,Ik,temp_Ikl,temp_Iku);
    // Passo 2 //
    double soma_delta_l = 0;
    double soma_delta_u = 0;
    int count_temp_Ikl = 0;
    int count_temp_Ukl = 0;
    
    int steps = 0;

    checa_viabilidade(p,qk,&qk_l,&qk_u,Ik,x,lambda,temp_Ikl,temp_Iku,
            &soma_delta_l,&soma_delta_u,&count_temp_Ikl,&count_temp_Ukl);
    // Fim passo 2 //
    while(soma_delta_l != soma_delta_u){

        // Passo 4 //
        if (soma_delta_l > soma_delta_u){
            #pragma omp parallel for if (count_temp_Ikl > limit_for)
                for (unsigned i = 0; i < count_temp_Ikl;i++){
                    // printf("ordem %d\n",temp_Ikl[i]);
                    Ik[temp_Ikl[i]] = -1;
                }
            lambda += (soma_delta_l / qk_l) ;
            qk = qk_l;
        }else{
            #pragma omp parallel for if (count_temp_Ukl > limit_for)
                for (unsigned i = 0; i < count_temp_Ukl;i++){
                    Ik[temp_Iku[i]] = -1;
                }
            lambda += (-soma_delta_u / qk_u) ;
            qk = qk_u;
        }
        // FIM Passo 4 //

        // Passo 2 //
        soma_delta_l = 0;
        soma_delta_u = 0;
        count_temp_Ikl = 0;
        count_temp_Ukl = 0;
        checa_viabilidade(p,qk,&qk_l,&qk_u,Ik,x,lambda,temp_Ikl,temp_Iku,
            &soma_delta_l,&soma_delta_u,&count_temp_Ikl,&count_temp_Ukl);
        //FIM passo 2

        steps++;
    }
    // imprime_resultado_x(p,x);
    return steps;
    return 0;

}