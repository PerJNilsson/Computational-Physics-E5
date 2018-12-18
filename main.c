#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>

#define PI 3.14159265358979323846
int main () {
        double rnd;
        int N = 10;
        double *grid_pointer = malloc(sizeof(double*)*N*N);
        double **grid = malloc(sizeof(double)*N);
        for (int i=0, j=0; i<N; i++, j+=N){
                grid[i] =grid_pointer + j;
        }
        double * b = malloc(sizeof(double)*N);
        double * u = malloc(sizeof(double)*N);

        // Set up the A matrix
        grid[0][0] = 0; grid[N-1][N-1]=1;
        for (int i=1; i<N-1; i++){
                for (int j=1; j<N-1; j++){
                        if (i==j){
                                grid[i][j-1] = -1;
                                grid[i][j]   =  2;
                                grid[i][j+1] = -1;
                        }
                }
        }
         // GSL INITIALIZATION
        const gsl_rng_type *T;
        gsl_rng *q;
        // Initializations
        gsl_rng_env_setup();
        T = gsl_rng_default;
        q = gsl_rng_alloc(T);
        gsl_rng_set(q,time(NULL));
  
        for (size_t i=1; i<N-1; i++) {
                rnd = gsl_rng_uniform(q);
                u[i] = (i/(double) N)*(i/(double) N)*4*PI*rnd;
        }
        
          for (int i=0; i<N; i++){
                for (int j=0; j<N; j++){
                        printf("%.0f",grid[i][j]);
                }
                
                printf("\n");
          }



        return 0;
}
