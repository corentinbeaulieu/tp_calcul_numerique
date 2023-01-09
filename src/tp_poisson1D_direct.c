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
    // nbpoints -> nb rows of A,                          la -> nb columns of A 
    int nbpoints, la;
    // ku -> nb upper diags of A,                         kl -> nb lower diags of A
    // kv -> nb additional diags for LU computation,      lab -> ld of AB ie nb cols here 
    int ku, kl, kv, lab;
    int *ipiv;
    int NRHS;
    // T0 -> initial solution,                            T1 -> final solution
    double T0, T1;
    double *RHS, *EX_SOL, *X;
    double *AB;

    double elapsed, relres;
    struct timespec before, after;

    NRHS=1;
    nbpoints=12;
    la=nbpoints-2;

    /* Dirichlet Boundary conditions */
    T0=5.0;
    T1=20.0;

    printf("--------- Poisson 1D ---------\n\n");
    RHS    = (double *) malloc(sizeof(double)*la);
    if(RHS == NULL) {
        perror("RHS ");
        return ERR_ALLOC;
    }
    EX_SOL = (double *) malloc(sizeof(double)*la);
    if(EX_SOL == NULL) {
        perror("EX_SOL ");
        return ERR_ALLOC;
    }
    X      = (double *) malloc(sizeof(double)*la);
    if(X == NULL) {
        perror("X ");
        return ERR_ALLOC;
    }

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
    if(AB == NULL) {
        perror("AB ");
        return ERR_ALLOC;
    }

    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

    /* write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat"); */

    printf("Solution with LAPACK\n");
    /* LU Factorization */
    ipiv = (int *) calloc(la, sizeof(int));
    if(ipiv == NULL) {
        perror("ipiv ");
        return ERR_ALLOC;
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &before);
    dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &ierr);

    /* LU for tridiagonal matrix  (can replace dgbtrf_) */
    // ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

    /* write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU1.dat"); */ 

    /* Solution (Triangular) */
    if (ierr==0){
        dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &ierr);
        if (ierr!=0){
            fprintf(stderr,"\n INFO DGBTRS = %d\n",ierr);
        }
    }
    else{
        fprintf(stderr,"\n INFO DGBTRF = %d\n",ierr);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &after);

    elapsed = (double) after.tv_sec + (double) after.tv_nsec /1000000000 - ((double) before.tv_sec + (double) before.tv_nsec /1000000000);

    write_xy(X, RHS, &la, "SOL_dgbtrf.dat");

    /* Relative forward error */
    cblas_dscal(la, -1.0, RHS, 1);
    cblas_daxpy(la, 1.0, EX_SOL, 1, RHS, 1);
    relres = cblas_dnrm2(la, RHS, 1) / cblas_dnrm2(la, EX_SOL, 1);

    printf("\nThe relative forward error for dgbrf+dgbrs is : %1.3e\n\
            The execution time is : %4.6f s\n",relres, elapsed);

    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);

    clock_gettime(CLOCK_MONOTONIC_RAW, &before);
    dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &ierr);
    if(ierr != 0) {
        fprintf(stderr, "INFO DGBSV = %d\n", ierr);
    }
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
    dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &ierr);
    /* write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat"); */ 

    /* Solution (Triangular) */
    if (ierr==0){
        dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &ierr);
        if (ierr!=0){
            fprintf(stderr, "\n INFO DGBTRS = %d\n", ierr);
        }
    }
    else{
        fprintf(stderr, "\n INFO DGBTRFTRIDIAG = %d\n",ierr);
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &after);

    elapsed = (double) after.tv_sec + (double) after.tv_nsec /1000000000 - ((double) before.tv_sec + (double) before.tv_nsec /1000000000);

    cblas_dscal(la, -1.0, RHS, 1);
    cblas_daxpy(la, 1.0, EX_SOL, 1, RHS, 1);
    relres = cblas_dnrm2(la, RHS, 1) / cblas_dnrm2(la, EX_SOL, 1);

    printf("\nThe relative forward error for dgbtrftridiag is : %1.3e\n\
            The execution time is : %4.6f s\n",relres, elapsed);


    /* Evaluation of performance */

    printf("\n\n ** Performance Analysis **\n");

    // This parameter can be changed
    // Be cautious exponential => expensive in memory
    const size_t nbIter = 20;

    int *tabTaille = aligned_alloc(32, nbIter * sizeof(int));
    if(tabTaille == NULL) {
        perror("tabTaille ");
        return ERR_ALLOC;
    }
    double *elapsed_dgbtrf = aligned_alloc(32, nbIter * sizeof(double)); 
    if(elapsed_dgbtrf == NULL) {
        perror("elapsed_dgbtrf ");
        return ERR_ALLOC;
    }
    double *elapsed_dgbsv = aligned_alloc(32, nbIter * sizeof(double)); 
    if(elapsed_dgbsv == NULL) {
        perror("elapsed_dgbsv ");
        return ERR_ALLOC;
    }
    double *elapsed_dgbtrftridiag = aligned_alloc(32, nbIter * sizeof(double)); 
    if(elapsed_dgbtrftridiag == NULL) {
        perror("elapsed_dgbtrftridiag ");
        return ERR_ALLOC;
    }

    nbpoints = 2;

    for(size_t i = 0; i  < nbIter; i++) {

        nbpoints *= 2;
        la = nbpoints - 2;

        RHS = (double *) realloc(RHS, sizeof(double)*la);
        if(RHS == NULL) {
            perror("RHS ");
            return ERR_ALLOC;
        }
        AB  = (double *) realloc(AB,  sizeof(double)*lab*la);
        if(AB == NULL) {
            perror("AB ");
            return ERR_ALLOC;
        }
        ipiv = (int *) realloc(ipiv, la*sizeof(int));
        if(ipiv == NULL) {
            perror("ipiv ");
            return ERR_ALLOC;
        }

        tabTaille[i] = nbpoints;


        /* dgbtrf + dgtrs measure */

        set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
        set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);

        clock_gettime(CLOCK_MONOTONIC_RAW, &before);

        dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &ierr);

        if (ierr==0){
            dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &ierr);
            if (ierr!=0){
                fprintf(stderr, "\n INFO DGBTRS = %d\n",ierr);}
        }
        else{
            fprintf(stderr, "\n INFO DGBTRF = %d\n",ierr);
        }

        clock_gettime(CLOCK_MONOTONIC_RAW, &after);

        elapsed_dgbtrf[i] = (double) after.tv_sec + (double) after.tv_nsec /1000000000 - ((double) before.tv_sec + (double) before.tv_nsec /1000000000);


        /* dgbsv  measure */

        set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
        set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);

        clock_gettime(CLOCK_MONOTONIC_RAW, &before);

        dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &ierr);
        if (ierr!=0){
            fprintf(stderr, "\n INFO DGBSV = %d\n", ierr);
        }

        clock_gettime(CLOCK_MONOTONIC_RAW, &after);

        elapsed_dgbsv[i] = (double) after.tv_sec + (double) after.tv_nsec /1000000000 - ((double) before.tv_sec + (double) before.tv_nsec /1000000000);
        clock_gettime(CLOCK_MONOTONIC_RAW, &before);


        /* dfbtrftridiag measure */

        set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
        set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);

        dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &ierr);

        if (ierr==0){
            dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &ierr);
            if (ierr!=0){
                fprintf(stderr, "\n INFO DGBTRS = %d\n",ierr);}
        }
        else{
            fprintf(stderr, "\n INFO DGBTRFTRIDIAG = %d\n",ierr);
        }

        clock_gettime(CLOCK_MONOTONIC_RAW, &after);

        elapsed_dgbtrftridiag[i] = (double) after.tv_sec + (double) after.tv_nsec /1000000000 - ((double) before.tv_sec + (double) before.tv_nsec /1000000000);
    }


    int fd = open("TIME_DIRECT.dat", O_CREAT | O_WRONLY | O_TRUNC, 0666);
    if(fd == -1) {
        perror("open ");
    }

    else {
        for(size_t i = 0; i < nbIter; i++) {
            dprintf(fd, "%9d; %4.6lf; %4.6lf; %4.6lf\n", 
                    tabTaille[i],
                    elapsed_dgbtrf[i],
                    elapsed_dgbsv[i],
                    elapsed_dgbtrftridiag[i]);
        }

        close(fd); 
    }

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
