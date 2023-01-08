/******************************************/
/* tp2_poisson1D_iter.c                   */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include "atlas_headers.h"

int main(void)
{
    int ierr;
    int nbpoints, la;
    int ku, kl, lab, kv;
    int *ipiv = NULL;
    double T0, T1;
    double *RHS, *SOL, *EX_SOL, *X;
    double *AB;
    double *MB;

    double relres;

    double opt_alpha;

    /* Size of the problem */
    nbpoints=12;
    la=nbpoints-2;

    /* Dirichlet Boundary conditions */
    T0=5.0;
    T1=20.0;

    printf("--------- Poisson 1D ---------\n\n");
    RHS=(double *) malloc(sizeof(double)*la);
    SOL=(double *) calloc(la, sizeof(double)); 
    EX_SOL=(double *) malloc(sizeof(double)*la);
    X=(double *) malloc(sizeof(double)*la);

    /* Setup the Poisson 1D problem */
    /* General Band Storage */
    set_grid_points_1D(X, &la);
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

    /* uncomment the following to check analytical solution */
    // write_vec(RHS, &la, "RHS.dat");
    // write_vec(EX_SOL, &la, "EX_SOL.dat");
    // write_vec(X, &la, "X_grid.dat");

    kv=0;
    ku=1;
    kl=1;
    lab=kv+kl+ku+1;

    AB = (double *) malloc(sizeof(double)*lab*la);
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

    /* uncomment the following to check matrix A */
    // write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

    /********************************************/
    /* Solution (Richardson with optimal alpha) */

    /* Computation of optimum alpha */
    opt_alpha = richardson_alpha_opt(&la);
    printf("Optimal alpha for simple Richardson iteration is : %.3lf\n",opt_alpha); 

    // uncomment to compute and write AB eigvalues
    /* double *eigval = aligned_alloc(32, la*sizeof(double));
     * eig_poisson1D(eigval, &la);
     * write_vec(eigval, &la, "EIGVAL.dat");
     *
     * free(eigval); */

    /* Solve */
    double tol=1e-4;
    int maxit=200000;
    double *resvec;
    int nbite=0;

    resvec = (double *) calloc(maxit + 1, sizeof(double));

    /* Solve with Richardson alpha */
    richardson_alpha(AB, RHS, SOL, &opt_alpha, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
    printf("\nNumber of iteration of Richardson is   : %5d\n", nbite);

    /* uncomment to Write richardson_alpha solution */
    // write_vec(SOL, &la, "SOL_RICH.dat");

    /* Write convergence history */
    // The increment of nbite is to write final backward error as well  
    // We have nbite + 1 res [0; nbite]
    nbite++;
    write_vec(resvec, &nbite, "RESVEC_RICH.dat");

    /* Compute forward error */
    cblas_dscal(la, -1.0, SOL, 1);
    cblas_daxpy(la, 1.0, EX_SOL, 1, SOL, 1);
    relres = cblas_dnrm2(la, SOL, 1) / cblas_dnrm2(la, EX_SOL, 1);
    printf("The relative forward  error for Richardson is    : %1.3e\nThe relative backward error for Richardson is    : %1.3e\n",
            relres,
            resvec[nbite-1]);

    /********************************************************/
    /*** Convergence analysis ***/

    const int nbAlpha = 100;

    // The method diverge for alpha >= 0.52 at least 
    const double h = (double) 0.52/nbAlpha;
    int *convergence = aligned_alloc(32, nbAlpha * sizeof(int));
    double *tabAlpha = aligned_alloc(32, nbAlpha * sizeof(double));
    double alpha = 0.0;

    /* Compute the solution for different alpha */
    for(int i = 0; i < nbAlpha; i++) {
        nbite = 0;

        set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
        memset(SOL, 0, la * sizeof(double));
        memset(resvec, 0xff, (maxit + 1) * sizeof(double));

        alpha = i*h;
        tabAlpha[i] = alpha;

        richardson_alpha(AB, RHS, SOL, &alpha, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);

        convergence[i] = nbite;
    }

    /* Write convergence analysis results */
    FILE *file = fopen("ALPHA.dat", "w");
    for(int i = 0; i < nbAlpha; i++) {

        fprintf(file, "%.16f; %5d\n", tabAlpha[i], convergence[i]);

    }

    fclose(file);
    free(tabAlpha);
    free(convergence);


    /********************************************************/
    /*** Richardson General Tridiag ***/

    /* get MB (:=M, D for Jacobi, (D-E) for Gauss-seidel) */
    kv = 1;
    ku = 1;
    kl = 1;


    MB = (double *) malloc(lab * la * sizeof(double));

    /** Jacobi **/

    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    memset(SOL, 0, la * sizeof(double));
    memset(resvec, 0xff, (maxit + 1) * sizeof(double));

    extract_MB_jacobi_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);

    // uncomment to Write D 
    // write_GB_operator_colMajor_poisson1D(MB, &lab, &la, "D.dat");

    richardson_MB(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);

    printf("\nNumber of iteration of Jacobi is       : %5d\n", nbite);

    // uncomment to Write Jacobi solution */
    // write_vec(SOL, &la, "SOL_JACOBI.dat");

    /* Write convergence history */
    // The increment of nbite is to write final backward error as well  
    // We have nbite + 1 res [0; nbite]
    nbite++;
    write_vec(resvec, &nbite, "RESVEC_JACOBI.dat");

    /* Compute and print forward and backward errors */
    cblas_dscal(la, -1.0, SOL, 1);
    cblas_daxpy(la, 1.0, EX_SOL, 1, SOL, 1);
    relres = cblas_dnrm2(la, SOL, 1) / cblas_dnrm2(la, EX_SOL, 1);
    printf("The relative forward  error for Jacobi is        : %1.3e\nThe relative backward error for Jacobi is        : %1.3e\n",
            relres,
            resvec[nbite-1]);

    /** Gauß-Seidel **/

    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    memset(SOL, 0, la * sizeof(double));
    memset(resvec, 0xff, (maxit + 1) * sizeof(double));

    extract_MB_gauss_seidel_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);

    // uncomment to Write DE 
    // write_GB_operator_colMajor_poisson1D(MB, &lab, &la, "DE.dat");

    richardson_MB(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);

    printf("\nNumber of iteration of Gauss-Seidel is : %5d\n", nbite);

    // uncomment to Write Jacobi solution */
    // write_vec(SOL, &la, "SOL_GAUSS_SEIDEL.dat");

    /* Write convergence history */
    // The increment of nbite is to write final backward error as well  
    // We have nbite + 1 res [0; nbite]
    nbite++;
    write_vec(resvec, &nbite, "RESVEC_GAUSS_SEIDEL.dat");

    /* Compute and print forward and backward errors */
    cblas_dscal(la, -1.0, SOL, 1);
    cblas_daxpy(la, 1.0, EX_SOL, 1, SOL, 1);
    relres = cblas_dnrm2(la, SOL, 1) / cblas_dnrm2(la, EX_SOL, 1);
    printf("The relative forward  error for Gauss-Seidel is  : %1.3e\nThe relative backward error for Gauss-Seidel is  : %1.3e\n",
            relres,
            resvec[nbite-1]);

    /* Evaluation of performance */
    // Be cautious exponential augmentation => expensive in memory
    const int nbMeasure = 12;
    int *tabTaille = (int *) aligned_alloc(32, nbMeasure * sizeof(int));
    double *elapsed_alpha = (double *) aligned_alloc(32, nbMeasure * sizeof(double));
    double *elapsed_jacobi = (double *) aligned_alloc(32, nbMeasure * sizeof(double));
    double *elapsed_gauss_seidel = (double *) aligned_alloc(32, nbMeasure * sizeof(double));

    int *nbIte_alpha = (int *) aligned_alloc(32, nbMeasure * sizeof(int));
    int *nbIte_jacobi = (int *) aligned_alloc(32, nbMeasure * sizeof(int));
    int *nbIte_gauss_seidel = (int *) aligned_alloc(32, nbMeasure * sizeof(int));

    struct timespec before, after;

    nbpoints = 2;

    for(size_t i = 0; i  < nbMeasure; i++) {

        nbpoints *= 2;
        la = nbpoints - 2;
        tabTaille[i] = nbpoints;

        RHS = (double *) realloc(RHS, sizeof(double) * la);
        SOL = (double *) realloc(SOL, sizeof(double) * la); 
        AB  = (double *) realloc(AB,  sizeof(double) * lab * la);
        MB  = (double *) realloc(MB,  sizeof(double) * lab * la);

        kv = 0;

        set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
        set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
        /** Richardson Alpha opt **/

        memset(SOL, 0, la * sizeof(double));
        memset(resvec, 0, (maxit + 1) * sizeof(double));
        nbite = 0;

        clock_gettime(CLOCK_MONOTONIC_RAW, &before);
        richardson_alpha(AB, RHS, SOL, &opt_alpha, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
        clock_gettime(CLOCK_MONOTONIC_RAW, &after);

        elapsed_alpha[i] = (double) after.tv_sec + (double) after.tv_nsec /1000000000 - ((double) before.tv_sec + (double) before.tv_nsec /1000000000);

        nbIte_alpha[i] = nbite;

        /** Jacobi **/
        kv = 1;

        memset(SOL, 0, la * sizeof(double));
        memset(resvec, 0xff, (maxit + 1) * sizeof(double));
        nbite = 0;

        clock_gettime(CLOCK_MONOTONIC_RAW, &before);
        extract_MB_jacobi_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);

        richardson_MB(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
        clock_gettime(CLOCK_MONOTONIC_RAW, &after);

        elapsed_jacobi[i] = (double) after.tv_sec + (double) after.tv_nsec /1000000000 - ((double) before.tv_sec + (double) before.tv_nsec /1000000000);

        nbIte_jacobi[i] = nbite;
        
       /** Gauss-Seidel **/

        memset(SOL, 0, la * sizeof(double));
        memset(resvec, 0xff, (maxit + 1) * sizeof(double));
        nbite = 0;

        clock_gettime(CLOCK_MONOTONIC_RAW, &before);
        extract_MB_gauss_seidel_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);

        richardson_MB(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
        clock_gettime(CLOCK_MONOTONIC_RAW, &after);
        elapsed_gauss_seidel[i] = (double) after.tv_sec + (double) after.tv_nsec /1000000000 - ((double) before.tv_sec + (double) before.tv_nsec /1000000000);

        nbIte_gauss_seidel[i] = nbite;
    }

    /* Write performance analysis results */
    file = fopen("TIME_ITER.dat", "w");
    for(int i = 0; i < nbMeasure; i++) {

        fprintf(file, "%8d; %4.6lf; %4.6lf; %4.6lf\n",
                tabTaille[i],
                elapsed_alpha[i],
                elapsed_jacobi[i],
                elapsed_gauss_seidel[i]);
    }

    fclose(file);
    free(elapsed_alpha);
    free(elapsed_jacobi);
    free(elapsed_gauss_seidel);


    file = fopen("ITER.dat", "w");
    for(int i = 0; i < nbMeasure; i++) {

        fprintf(file, "%10d; %5d; %5d; %5d\n",
                tabTaille[i],
                nbIte_alpha[i],
                nbIte_jacobi[i],
                nbIte_gauss_seidel[i]);
    }

    fclose(file);
    free(tabTaille);
    free(nbIte_alpha);
    free(nbIte_jacobi);
    free(nbIte_gauss_seidel);

    free(ipiv);
    free(resvec);
    free(RHS);
    free(SOL);
    free(EX_SOL);
    free(X);
    free(AB);
    free(MB);
    printf("\n--------- End -----------\n");
}
