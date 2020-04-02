#include <stdio.h> 
#include <stdlib.h> 
#include <time.h> 

#include "fix_variable_kiwiel_jonatas.c"

/* Possible benchmark types */
typedef enum {
    correlated,
    weakly_correlated,
    uncorrelated,
    flow,
    invalid
} problem_type;


double r2()
{
    return (double)rand() / (double)RAND_MAX ;
}

void printRandoms(int lower, int upper,  
                             int count) 
{ 
    int i; 
    for (i = 0; i < count; i++) { 
        double num = (rand() % 
           (upper - lower + 1)) + lower; 
        printf("%f ", num); 
    } 
} 

double generatorNumber(int lower, int upper) 
{ 
    return (double)(rand() % (upper - lower + 1)) + lower; 
} 


void generate_cqk_problem(problem_type type, cqk_problem *restrict p) {
    double
        temp,
        lb = 0.0, /* Receives p->b^t p->low. */
        ub = 0.0; /* Receives p->b^t p->up. */
    int lower = 10, upper = 25; 
    /* Generate the problem */

    // Use current time as  
    // seed for random generator 
    srand(time(0)); 

    for (unsigned i = 0; i < p->n; ++i) {


        p->b[i] = 10 + 14*r2();
        p->low[i] = 1 + 14*r2();
        p->up[i] = 1 + 14*r2();
        if (p->low[i] > p->up[i]) {
            temp = p->low[i];
            p->low[i] = p->up[i];
            p->up[i] = temp;
        }
        lb += p->b[i]*p->low[i];
        ub += p->b[i]*p->up[i];
        switch (type) {

        case uncorrelated:
            p->d[i] = 10 + 14*r2();
            p->a[i] = (10 + 14*r2());
            break;

        case weakly_correlated:
            p->d[i] = p->b[i] - generatorNumber(10,25);
            p->a[i] = p->b[i] - generatorNumber(10,25);
            break;

        case correlated:
        default:
            p->d[i] = p->b[i] + 5.0;
            p->a[i] = (p->b[i] + 5.0);
        }
    }
    p->r = lb + (ub - lb)*0.1;
}

void open_file(cqk_problem*restrict p,int n){
    FILE *myFile;
    myFile = fopen("instance_test.txt", "r");
    int len_file = (5*n)+1;
    double numberArray[len_file];


    for (int i = 0; i < len_file; i++)
    {
        fscanf(myFile, "%lf", &numberArray[i]);
    }
    // for (int i = 0; i < len_file; i++)
    // {
    //     // printf("%f ", numberArray[i]);
    // }
    for (int i = 0; i < len_file-1; i++)
    {
        if(i < n)
        {
            p->d[i] = numberArray[i];
        }else if (i < n*2)
        {
            p->a[i-n] = numberArray[i];
        }else if (i < n*3)
        {
            p->b[i-(n*2)] = numberArray[i];
        }else if (i < n*4)
        {
            p->up[i-(n*3)] = numberArray[i];
        }else{
            p->low[i-(n*4)] =numberArray[i];
        }
    }
    p->r = numberArray[len_file-1];
}

void run_a_test(cqk_problem *restrict p, double *restrict x) {
    
    clock_t start1, end1;
    double cpu_time_used1;


    double max = 0;
    double min = 999999;
    int repeticao = 1;
    double duracao = 0.0;
    for (int i = 0;i<repeticao;i++){
        int status = -1;
        start1 = clock();
        status = fix_variable(p, x);
        end1 = clock();
        cpu_time_used1 = ((double) (end1 - start1)) / CLOCKS_PER_SEC;
        duracao += cpu_time_used1;
        if (cpu_time_used1 > max){
            max = cpu_time_used1;

        }
        if (cpu_time_used1 < min){
            min = cpu_time_used1;
        }
        printf("Time Total: %f \n", cpu_time_used1);
        printf("result = %d\n",status);
    }
    printf("media time %f \n",duracao/repeticao);
    printf("min time %f \n",min);
    printf("max time %f \n",max);

    //Executa algoritmo JONATAS
    // status = fix_variable(p, x);
    //Executa algoritmo SILVA
    // status = kiwiel_var_fix(p, x);
}

int main(int argc, char *argv[]) {

    problem_type ptype; /* Code for the problem type */
    unsigned n;         /* Problem dimension */ 
    n = 2000000;
    ptype = correlated;

    /* Problem data */
    cqk_problem p;
    allocate_cqk_problem(n, &p);

    /* Solution vector Simplex */
    double *x = (double *) malloc(p.n*sizeof(double));
    generate_cqk_problem(ptype, &p);
    // open_file(&p,n);

    
    run_a_test(&p,x);
    


    free_cqk_problem(&p);
    return 0;

}