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

int fix_variable(cqk_problem *restrict p,double *x){

    printf("r %f \n",p->r);


    // Passo 0 //
    int *Ik = (int *) malloc(p->n*sizeof(int));
    int *Uk = (int *) malloc(p->n*sizeof(int));
    int *Lk = (int *) malloc(p->n*sizeof(int));
    int *temp_Ikl = (int *) malloc(p->n*sizeof(int));
    int *temp_Iku = (int *) malloc(p->n*sizeof(int));

    for (int i = 0; i<p->n;i++){
        Ik[i] = 0;
        Uk[i] = -1;
        Lk[i] = -1;
        temp_Ikl[i] = -1;
        temp_Iku[i] = -1;
        x[i] = 0;
    }
    int k = 1;
    
    // FIM Passo 0 //

    //variavel utilizada apenas para contar quantidade de variaveis fixadas
    int count_fixes = 0;

    for (int steps = 0; steps < 10;steps++){

        printf("\n ------------------ #### ----------------- \n");
        printf("step %d \n", steps);

        // Passo 1 // 
        double soma_bl = 0;

        

        for (int i = 0;i < p->n; i++){
            if (Lk[i] != -1){
                soma_bl = soma_bl + (p->low[i]*p->b[i]);
            }
        }
        double soma_bu = 0;
        for (int i = 0;i < p->n; i++){
            if (Uk[i] != -1){
                soma_bu = soma_bu + (p->up[i]*p->b[i]);
            }
        }
        double rk = p->r - soma_bl - soma_bu;
        // printf("step %d rk %f\n",steps,rk);
        double soma_ab_by_d = 0;
        double soma_b_by_d = 0;
        
        for(int i = 0; i < p->n;i++){
            if (Ik[i] != -1){
                // printf("Ik = %d\n",i);
                soma_ab_by_d = soma_ab_by_d + ((p->a[i]*p->b[i])/p->d[i]);
                soma_b_by_d = soma_b_by_d + (pow(p->b[i],2)/p->d[i]);
            }
        }
        
        double tk = (soma_ab_by_d - rk) / soma_b_by_d;
        // printf("step %d tk %f\n",steps,tk);
        double soma_delta_l = 0;
        double soma_delta_u = 0;
        
        // FIM Passo 1 //

        // Passo 2 //
        int count_temp_Ikl = 0;
        int count_temp_Ukl = 0;

        for(int i = 0; i < p->n;i++){
            if (Ik[i] != -1){

                double xk = ((p->a[i] - tk*p->b[i]) / p->d[i]);
                // printf("%d xk %f\n",i,xk);
                if (xk < p->low[i]){
                    // printf("menor %d - %d \n",i,count_temp_Ikl);
                    x[i] = p->low[i];
                    temp_Ikl[count_temp_Ikl] = i;
                    soma_delta_l = soma_delta_l + p->b[i] * (p->low[i] - xk);
                    count_temp_Ikl++;

                }else if (xk > p->up[i]){
                    // printf("maior \n");
                    x[i] = p->up[i];
                    temp_Iku[count_temp_Ukl] = i;
                    soma_delta_u = soma_delta_u + p->b[i] * (xk - p->up[i]);
                    count_temp_Ukl++;
                }else{
                    // printf("igual \n");
                    x[i] = xk;
                    // Ik[i] = -1;
                }
            }
            
        }

        // FIM Passo 2 //

        // Passo 3 //
        if (soma_delta_l == soma_delta_u){
            // double total = 0;
            // for (int i = 0; i < p->n;i++){
            //  	total = total + x[i]*p->b[i];
            //      printf("x%d = %f  l = %f   u = %f  \n",i,x[i],p->low[i],p->up[i]);
            // }
            // printf("fixies = %d \n",count_fixes);
            // printf("r %f\n",p->r);
            // printf("total %f \n",total);
            break;
        }
        // FIM Passo 3 //

        // Passo 4 //
        if (soma_delta_l > soma_delta_u){
            printf("Qunatidade indices menor = %d\n",count_temp_Ikl);
            for (int i = 0; i < count_temp_Ikl;i++){
                Ik[temp_Ikl[i]] = -1;
                Lk[temp_Ikl[i]] = 0;
                // printf("- %d / %d - ",temp_Ikl[i],Lk[temp_Ikl[i]]);
            }
            count_fixes = count_fixes + count_temp_Ikl;
        }
        if (soma_delta_l < soma_delta_u){
            printf("Qunatidade indices maior = %d\n",count_temp_Ukl);
            for (int i = 0; i < count_temp_Ukl;i++){
                Ik[temp_Iku[i]] = -1;
                Uk[temp_Iku[i]] = 0;
                // printf(" - %d - ",temp_Iku[i]);
            }
            count_fixes = count_fixes + count_temp_Ukl;
        }
        // FIM Passo 4 //
    }

    return 0;

}




int fix_variable_new(cqk_problem *restrict p,double *x){

    printf("r %f \n",p->r);


    // Passo 0 //
    int *Ik = (int *) malloc(p->n*sizeof(int));
    int *sobraIk = (int *) malloc(p->n*sizeof(int));
    int *Uk = (int *) malloc(p->n*sizeof(int));
    int *Lk = (int *) malloc(p->n*sizeof(int));
    int *temp_Ikl = (int *) malloc(p->n*sizeof(int));
    int *temp_Iku = (int *) malloc(p->n*sizeof(int));

    for (int i = 0; i<p->n;i++){
        Ik[i] = 0;
        sobraIk[i] = 0;
        Uk[i] = -1;
        Lk[i] = -1;
        temp_Ikl[i] = -1;
        temp_Iku[i] = -1;
        x[i] = 0;
    }
    int k = 1;
    
    // FIM Passo 0 //

    //variavel utilizada apenas para contar quantidade de variaveis fixadas
    int count_fixes = 0;

    //tamanho das sobras Ik
    int tam_sobraIk = p->n;

    for (int steps = 0; steps < 10;steps++){

        printf("\n ------------------ #### ----------------- \n");
        printf("step %d \n", steps);

        // Passo 1 // 
        double soma_bl = 0;

        

        for (int i = 0;i < p->n; i++){
            if (Lk[i] != -1){
                soma_bl = soma_bl + (p->low[i]*p->b[i]);
            }
        }
        double soma_bu = 0;
        for (int i = 0;i < p->n; i++){
            if (Uk[i] != -1){
                soma_bu = soma_bu + (p->up[i]*p->b[i]);
            }
        }
        double rk = p->r - soma_bl - soma_bu;
        // printf("step %d rk %f\n",steps,rk);
        double soma_ab_by_d = 0;
        double soma_b_by_d = 0;
        
        // for(int i = 0; i < p->n;i++){
        //     if (Ik[i] != -1){
        //         soma_ab_by_d = soma_ab_by_d + ((p->a[i]*p->b[i])/p->d[i]);
        //         soma_b_by_d = soma_b_by_d + (pow(p->b[i],2)/p->d[i]);
        //     }
        // }

        if (steps==0)
        {
            for(int i = 0; i < p->n;i++)
            {
                if (Ik[i] != -1){
                    soma_ab_by_d = soma_ab_by_d + ((p->a[i]*p->b[i])/p->d[i]);
                    soma_b_by_d = soma_b_by_d + (pow(p->b[i],2)/p->d[i]);
                }
            }
        }else
        {
            for(int j = 0; j < tam_sobraIk;j++)
            {
                int i = sobraIk[j];
                // printf("Ik = %d\n",i);
                soma_ab_by_d = soma_ab_by_d + ((p->a[i]*p->b[i])/p->d[i]);
                soma_b_by_d = soma_b_by_d + (pow(p->b[i],2)/p->d[i]);

            }
        }


        double tk = (soma_ab_by_d - rk) / soma_b_by_d;
        // printf("step %d tk %f\n",steps,tk);
        double soma_delta_l = 0;
        double soma_delta_u = 0;
        
        // FIM Passo 1 //

        // Passo 2 //
        int count_temp_Ikl = 0;
        int count_temp_Ukl = 0;


        if (steps==0)
        {
            for(int i = 0; i < p->n;i++){
                if (Ik[i] != -1){

                    double xk = ((p->a[i] - tk*p->b[i]) / p->d[i]);
                    // printf("%d xk %f\n",i,xk);
                    if (xk < p->low[i]){
                        x[i] = p->low[i];
                        temp_Ikl[count_temp_Ikl] = i;
                        soma_delta_l = soma_delta_l + p->b[i] * (p->low[i] - xk);
                        count_temp_Ikl++;

                    }else if (xk > p->up[i]){
                        // printf("maior \n");
                        x[i] = p->up[i];
                        temp_Iku[count_temp_Ukl] = i;
                        soma_delta_u = soma_delta_u + p->b[i] * (xk - p->up[i]);
                        count_temp_Ukl++;
                    }else{
                        // printf("igual \n");
                        x[i] = xk;
                        // Ik[i] = -1;
                    }
                }
                
            }
        }else
        {
            for(int j = 0; j < tam_sobraIk;j++)
            {
                int i = sobraIk[j];
                double xk = ((p->a[i] - tk*p->b[i]) / p->d[i]);
                    // printf("%d xk %f\n",i,xk);
                if (xk < p->low[i]){
                    // printf("menor %d - %d \n",i,count_temp_Ikl);
                    x[i] = p->low[i];
                    temp_Ikl[count_temp_Ikl] = i;
                    soma_delta_l = soma_delta_l + p->b[i] * (p->low[i] - xk);
                    count_temp_Ikl++;

                }else if (xk > p->up[i]){
                    x[i] = p->up[i];
                    temp_Iku[count_temp_Ukl] = i;
                    soma_delta_u = soma_delta_u + p->b[i] * (xk - p->up[i]);
                    count_temp_Ukl++;
                }else{
                    x[i] = xk;
                    // Ik[i] = -1;
                }
            }
                
            
        }

        // FIM Passo 2 //

        // Passo 3 //
        if (soma_delta_l == soma_delta_u){
            // double total = 0;
            // for (int i = 0; i < p->n;i++){
            //  	total = total + x[i]*p->b[i];
            //      printf("x%d = %f  l = %f   u = %f  \n",i,x[i],p->low[i],p->up[i]);
            // }
            // printf("fixies = %d \n",count_fixes);
            // printf("r %f\n",p->r);
            // printf("total %f \n",total);
            break;
        }
        // FIM Passo 3 //

        // Passo 4 //
        if (soma_delta_l > soma_delta_u){
            printf("Qunatidade indices menor = %d\n",count_temp_Ikl);
            int count = 0;
            int count_sobra = 0;
            
            if(steps == 0){
                for (int i = 0; i < p->n;i++){
                
                    // printf("Ik = %d\n",i);

                    if(temp_Ikl[count] == i){
                        Ik[temp_Ikl[count]] = -1;
                        Lk[temp_Ikl[count]] = 0;
                        // printf("- %d / %d ",temp_Ikl[count],Lk[temp_Ikl[count]]);
                        count++;
                    }else{
                        sobraIk[count_sobra] = i;
                        // printf("\nsobraIk %d \n", sobraIk[count_sobra]);
                        count_sobra++;
                    }
                }
                tam_sobraIk = tam_sobraIk - count;
            }else{
                for(int j = 0; j < tam_sobraIk;j++)
                {
                    int i = sobraIk[j];
                    // printf("Ik = %d\n",i);

                    if(temp_Ikl[count] == i){
                        Ik[temp_Ikl[count]] = -1;
                        Lk[temp_Ikl[count]] = 0;
                        // printf("- %d / %d ",temp_Ikl[count],Lk[temp_Ikl[count]]);
                        count++;
                    }else{
                        sobraIk[count_sobra] = i;
                        // printf("\nsobraIk %d \n", sobraIk[count_sobra]);
                        count_sobra++;
                    }
                }
                tam_sobraIk = tam_sobraIk - count;
            }
            count_fixes = count_fixes + count_temp_Ikl;
            // printf("\nsobra %d \n",tam_sobraIk);
        }



        if (soma_delta_l < soma_delta_u){
            printf("Qunatidade indices maior = %d\n",count_temp_Ukl);
            int count = 0;
            int count_sobra = 0;
            
            if(steps == 0){
                for (int i = 0; i < p->n;i++){
                    if(temp_Iku[count] == i){
                        Ik[temp_Iku[count]] = -1;
                        Uk[temp_Iku[count]] = 0;
                        // printf(" - %d - ",temp_Iku[i]);
                        count++;
                    }else{
                        sobraIk[count_sobra] = i;
                        // printf("\nsobraIk %d \n", sobraIk[count_sobra]);
                        count_sobra++;
                    }
                }
                tam_sobraIk = tam_sobraIk - count;
            }else{
                for(int j = 0; j < tam_sobraIk;j++)
                {
                    int i = sobraIk[j];
                    // printf("Ik = %d\n",i);

                    if(temp_Iku[count] == i){
                        Ik[temp_Iku[count]] = -1;
                        Uk[temp_Iku[count]] = 0;
                        // printf("- %d / %d ",temp_Ikl[count],Lk[temp_Ikl[count]]);
                        count++;
                    }else{
                        sobraIk[count_sobra] = i;
                        // printf("\nsobraIk %d \n", sobraIk[count_sobra]);
                        count_sobra++;
                    }
                }
                tam_sobraIk = tam_sobraIk - count;
            }
            count_fixes = count_fixes + count_temp_Ukl;
        }
        
        
        
        
        
        // FIM Passo 4 //
    }

    return 0;

}