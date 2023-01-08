/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include "atlas_headers.h"


int main(void)
{
    int ierr;
    // nbpoints -> nb lignes de A,                             la -> nombre de colonnes de A
    int nbpoints, la;
    // ku -> nb de surdiagonales de A,                         kl -> nombre de sousdiagonales de A
    // kv -> nb sousdiagonales en plus pour le stockage de LU, lab -> ld de AB ie nombre de lignes de AB
    int ku, kl, kv, lab;
    int *ipiv;
    int info;
    int NRHS;
    // T0 -> solution initiale,                                T1 -> solution finale
    double T0, T1;
    double *RHS, *EX_SOL, *X;
    double *AB;

    double elapsed, relres;
    struct timespec before, after;

    NRHS=1;
    nbpoints=1000000;
    la=nbpoints-2;
    T0=-5.0;
    T1=5.0;

    printf("--------- Poisson 1D ---------\n\n");
    RHS=(double *) malloc(sizeof(double)*la);
    EX_SOL=(double *) malloc(sizeof(double)*la);
    X=(double *) malloc(sizeof(double)*la);

    // TODO : you have to implement those functions
    set_grid_points_1D(X, &la);
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

    /* write_vec(RHS, &la, "RHS.dat"); */
    /* write_vec(EX_SOL, &la, "EX_SOL.dat"); */
    /* write_vec(X, &la, "X_grid.dat"); */

    kv=1; 
    ku=1;
    kl=1;
    lab=kv+kl+ku+1;

    AB = (double *) malloc(sizeof(double)*lab*la);

    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

    /* write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat"); */

    printf("Solution with LAPACK\n");
    /* LU Factorization */
    info=0;
    ipiv = (int *) calloc(la, sizeof(int));


    clock_gettime(CLOCK_MONOTONIC_RAW, &before);


    dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

    /* LU for tridiagonal matrix  (can replace dgbtrf_) */
    // ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

    /* write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU1.dat"); */ 

    /* Solution (Triangular) */
    if (info==0){
        dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
        if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
    }else{
        printf("\n INFO = %d\n",info);
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &after);

    elapsed = (double) after.tv_sec + (double) after.tv_nsec /1000000000 - ((double) before.tv_sec + (double) before.tv_nsec /1000000000);
    // RHS now contains the solution x

    /* It can also be solved with dgbsv */
    // TODO : use dgbsv

    write_xy(X, RHS, &la, "SOL.dat");

    /* Relative forward error */
    // TODO : Compute relative norm of the residual
    cblas_dscal(la, -1.0, RHS, 1);
    cblas_daxpy(la, 1.0, EX_SOL, 1, RHS, 1);
    relres = cblas_dnrm2(la, RHS, 1) / cblas_dnrm2(la, EX_SOL, 1);

    printf("\nThe relative forward error for dgbrf+dgbrs is : %1.3e\n\
            The execution time is : %4.6f s\n",relres, elapsed);

    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);

    clock_gettime(CLOCK_MONOTONIC_RAW, &before);

    dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);

    clock_gettime(CLOCK_MONOTONIC_RAW, &after);

    elapsed = (double) after.tv_sec + (double) after.tv_nsec /1000000000 - ((double) before.tv_sec + (double) before.tv_nsec /1000000000);

    cblas_dscal(la, -1.0, RHS, 1);
    cblas_daxpy(la, 1.0, EX_SOL, 1, RHS, 1);
    relres = cblas_dnrm2(la, RHS, 1) / cblas_dnrm2(la, EX_SOL, 1);

    printf("\nThe relative forward error for dgbsv is : %1.3e\n\
            The execution time is : %4.6f s\n",relres, elapsed);

    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);

    clock_gettime(CLOCK_MONOTONIC_RAW, &before);

    ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

    /* write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat"); */ 

    /* Solution (Triangular) */
    if (info==0){
        dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
        if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
    }else{
        printf("\n INFO = %d\n",info);
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &after);

    elapsed = (double) after.tv_sec + (double) after.tv_nsec /1000000000 - ((double) before.tv_sec + (double) before.tv_nsec /1000000000);

    cblas_dscal(la, -1.0, RHS, 1);
    cblas_daxpy(la, 1.0, EX_SOL, 1, RHS, 1);
    relres = cblas_dnrm2(la, RHS, 1) / cblas_dnrm2(la, EX_SOL, 1);

    printf("\nThe relative forward error for dgbtrftridiag is : %1.3e\n\
            The execution time is : %4.6f s\n",relres, elapsed);

    // Be cautious exponential augmentation => expensive in memory
    const size_t nbIter = 26;
    int *tabTaille = aligned_alloc(32, nbIter * sizeof(int));
    double *elapsed_dgbtrf = aligned_alloc(32, nbIter * sizeof(double)); 
    double *elapsed_dgbsv = aligned_alloc(32, nbIter * sizeof(double)); 
    double *elapsed_dgbtrftridiag = aligned_alloc(32, nbIter * sizeof(double)); 

    nbpoints = 2;

    /* Evaluation of performance */
    for(size_t i = 0; i  < nbIter; i++) {

        nbpoints *= 2;
        la = nbpoints - 2;

        RHS = (double *) realloc(RHS, sizeof(double)*la);
        AB  = (double *) realloc(AB,  sizeof(double)*lab*la);
        ipiv = (int *) realloc(ipiv, la*sizeof(int));

        tabTaille[i] = nbpoints;


        /* dgbtrf + dgtrs measure */

        set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
        set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);

        clock_gettime(CLOCK_MONOTONIC_RAW, &before);

        dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

        if (info==0){
            dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
            if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
        }else{
            printf("\n INFO = %d\n",info);
        }

        clock_gettime(CLOCK_MONOTONIC_RAW, &after);

        elapsed_dgbtrf[i] = (double) after.tv_sec + (double) after.tv_nsec /1000000000 - ((double) before.tv_sec + (double) before.tv_nsec /1000000000);


        /* dgbsv  measure */

        set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
        set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);

        clock_gettime(CLOCK_MONOTONIC_RAW, &before);

        dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);

        clock_gettime(CLOCK_MONOTONIC_RAW, &after);

        elapsed_dgbsv[i] = (double) after.tv_sec + (double) after.tv_nsec /1000000000 - ((double) before.tv_sec + (double) before.tv_nsec /1000000000);
        clock_gettime(CLOCK_MONOTONIC_RAW, &before);


        /* dfbtrftridiag measure */

        set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
        set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);

        ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

        if (info==0){
            dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
            if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
        }else{
            printf("\n INFO = %d\n",info);
        }

        clock_gettime(CLOCK_MONOTONIC_RAW, &after);

        elapsed_dgbtrftridiag[i] = (double) after.tv_sec + (double) after.tv_nsec /1000000000 - ((double) before.tv_sec + (double) before.tv_nsec /1000000000);
        printf("%d\n", i);
    }


    int fd = open("TIME_DIRECT.dat", O_CREAT | O_WRONLY | O_TRUNC, 0666);

    for(size_t i = 0; i < nbIter; i++) {
        dprintf(fd, "%5d; %4.6lf; %4.6lf; %4.6lf\n", 
                tabTaille[i],
                elapsed_dgbtrf[i],
                elapsed_dgbsv[i],
                elapsed_dgbtrftridiag[i]);
    }

    /* close(fd); */
    free(tabTaille);
    free(elapsed_dgbtrf);
    free(elapsed_dgbsv);
    free(elapsed_dgbtrftridiag);

    free(ipiv);
    free(RHS);
    free(EX_SOL);
    free(X);
    free(AB);
    printf("\n--------- End -----------\n");
}
