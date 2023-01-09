/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include "atlas_headers.h"


// Creates a Poisson1D matrix in GB col major of size lab*la
void set_GB_operator_colMajor_poisson1D ( double* AB, int *lab, int *la, int *kv ) {

    for(int i = 0; i < *la; i++) {
        for(int j = 0; j < *(lab); j++) {

            if (j < *kv) {
                AB[i * *lab + j] = 0;
            }
            else if(j == *kv + 1) {
                AB[i * *lab + j] = 2;
            }
            else {
                AB[i * *lab + j] = -1;
            }
        }
    }

    AB[*kv] = 0;
    AB[(*lab) * (*la) -1] = 0;

}


// Creates an identity matrix in GB col major of size lab*la
void set_GB_operator_colMajor_poisson1D_Id ( double* AB, int *lab, int *la, int *kv ) {

    for(int i = 0; i < *la; i++) {
        for(int j = 0; j < *lab; j++) {
            if(j == *kv + 1) {
                AB[i * *lab + j] = 1;
            }
            else {
                AB[i * *lab + j] = 0;
            }
        }
    }
}


// Creates a RHS for a simple Poisson1D problem
void set_dense_RHS_DBC_1D ( double* RHS, int* la, double* BC0, double* BC1 ) {

    // We have g(i) = 0 forall i
    RHS[0] = *BC0;
    for(int i = 1; i < *la - 1; i++) {
        RHS[i] = 0;
    }
    RHS[*la - 1] = *BC1;

}  


// Resolves a Poisson1D problem in GB with the analytical solution
void set_analytical_solution_DBC_1D ( double* EX_SOL, double* X, int* la, double* BC0, double* BC1 ) {

    // The analytical solution is T(x) = T0 + x(T1-T0) if g(i) = 0 forall i as considered here 
    for(int i = 0; i < *la; i++) {

        EX_SOL[i] = *BC0 + X[i] * (*BC1 - *BC0);

    }
}  


// Creates a uniform discretisation of the interval (0;1)
void set_grid_points_1D ( double* X, int* la ) {

    // X(i) = i/la = ih with h a constant step
    // We have n-1 intervalles between n points and la = npoints - 2
    const double h = 1.0/ (double) (*la + 1); 

    for(int i = 1; i < *la + 1; i++) {

        X[i-1] = i*h;

    }
    // To test, we must have forall i, 0 < X(i) < 1 
    // More specifically, X(0) - h = 0 and X(la) + h = 1
}


void write_GB_operator_rowMajor_poisson1D ( double* AB, int* lab, int* la, char* filename ) {

    FILE * file;
    int ii,jj;
    file = fopen(filename, "w");
    //Numbering from 1 to la
    if (file != NULL){
        for (ii=0;ii<(*lab);ii++){
            for (jj=0;jj<(*la);jj++){
                fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
            }
            fprintf(file,"\n");
        }
        fclose(file);
    }
    else{
        perror(filename);
    }
}


void write_GB_operator_colMajor_poisson1D ( double* AB, int* lab, int* la, char* filename ) {

    FILE * file;
    int ii,jj;
    file = fopen(filename, "w");
    //Numbering from 1 to la
    if (file != NULL){
        for (ii=0;ii<(*la);ii++){
            for (jj=0;jj<(*lab);jj++){
                fprintf(file,"%02.6lf\t",AB[ii*(*lab)+jj]);
            }
            fprintf(file,"\n");
        }
        fclose(file);
    }
    else{
        perror(filename);
    }
}


void write_GB2AIJ_operator_poisson1D ( double* AB, int* la, char* filename ) {

    FILE * file;
    int jj;
    file = fopen(filename, "w");
    //Numbering from 1 to la
    if (file != NULL){
        for (jj=1;jj<(*la);jj++){
            fprintf(file,"%d\t%d\t%lf\n",jj,jj+1,AB[(*la)+jj]);
        }
        for (jj=0;jj<(*la);jj++){
            fprintf(file,"%d\t%d\t%lf\n",jj+1,jj+1,AB[2*(*la)+jj]);
        }
        for (jj=0;jj<(*la)-1;jj++){
            fprintf(file,"%d\t%d\t%lf\n",jj+2,jj+1,AB[3*(*la)+jj]);
        }
        fclose(file);
    }
    else{
        perror(filename);
    }
}


void write_vec ( double* vec, int* la, char* filename ) {

    int jj;
    FILE * file;
    file = fopen(filename, "w");
    // Numbering from 1 to la
    if (file != NULL){
        for (jj=0;jj<(*la);jj++){
            fprintf(file,"%.16lf\n",vec[jj]);
        }
        fclose(file);
    }
    else{
        perror(filename);
    } 
}  


void write_xy ( double* vec, double* x, int* la, char* filename ) {

    int jj;
    FILE * file;
    file = fopen(filename, "w");
    // Numbering from 1 to la
    if (file != NULL){
        for (jj=0;jj<(*la);jj++){
            fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
        }
        fclose(file);
    }
    else{
        perror(filename);
    } 
}  


int indexABCol(int i, int j, int *lab){
    return j * (*lab) * i;
}

// Resolves a linear system in GB where A is tridiagonal
int dgbtrftridiag ( int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info ) {

    // Error gestion inspired by lapack dgbsv's one
    if(*la < 0) {
        *info = -1;
    }
    else if(*n < 0) {
        *info = -2;
    }
    else if(*kl != 1) {
        *info = -3;
    }
    else if(*ku != 1) {
        *info = -4; 
    }
    else if(AB == NULL) {
        *info = -5;
    }
    else if(*lab < 3) {
        *info = -6;
    }
    else if(ipiv == NULL) {
        *info = -7;
    }
    else{

        // The algorithm is detailed in the report
        AB[*lab - 1] /= AB[*lab-2];
        for(int i = 1; i < *n; i++) {

            AB[(i+1) * (*lab) - 2] += AB[i *     (*lab) - 1];

            if(AB[(i+1) * (*lab) -2] == 0) {
                *info = i;
                break;
            }

            AB[(i+1) * (*lab) - 1] /= AB[(i+1) * (*lab) - 2];
            if(i == *n -1) *info = 0;

        }
    }

    return *info;
}


// Computes the eigvalues of a Poisson1d matrix of size la*la
void eig_poisson1D ( double* eigval, int *la ) {

    for(int i = 1; i <= *la; i++) {

        eigval[i] = 2.0 - 2.0*cos(i * M_PI / (*la +1));
    }
}


// Computes the maximum eigval of a Poisson1D of size la*la
double eigmax_poisson1D ( int *la ) {

    return 2.0 - 2.0*cos(*la * M_PI / (*la +1));
}


// Computes the minimum eigval of a Poisson1D of size la*la
double eigmin_poisson1D ( int *la ) {

    return 2.0 - 2.0*cos( M_PI / (*la +1));
}


// Computes optimum alpha for our problem of size la*la
double richardson_alpha_opt ( int *la ) {

    return 2.0 / (eigmax_poisson1D(la) + eigmin_poisson1D(la));
}


// Resolves a linear system with Richardson alpha method
int richardson_alpha ( double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite ) {

    int ret = 0;

    if      (AB == NULL)  {
        ret = -1;
    }
    else if (RHS == NULL) {
        ret = -2;
    }
    else if (X == NULL)   {
        ret = -3;
    }
    else if (*lab < *ku + *kl + 1) {
        ret = -5;
    }
    else if (*la <= 0) {
        ret = -6;
    }
    else if (*ku < 0) {
        ret = -7;
    }
    else if (*kl < 0) {
        ret = -8; 
    }
    else if (resvec == NULL) {
        ret = -11;
    }

    else {
        *nbite = 0;

        double *tmp = (double*) aligned_alloc(32, *la * sizeof(double));
        if(tmp == NULL) {
            perror("alloc richardson");
            exit(ERR_ALLOC);
        }

        cblas_dcopy(*la, RHS, 1, tmp, 1);

        // Compute first res
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, tmp, 1);
        resvec[*nbite] = cblas_dnrm2(*la, tmp, 1) / cblas_dnrm2(*la, RHS, 1); 

        while(resvec[*nbite] > *tol && *nbite < *maxit) {

            // If the method doesn't converge return an error
            if(*nbite > 1 && resvec[*nbite] > resvec[*nbite - 1]) {
                *nbite = -1;
                break;
            } 

            // Compute an iteration
            cblas_daxpy(*la, *alpha_rich, tmp, 1, X, 1); 

            // Compute the res
            (*nbite)++;

            cblas_dcopy(*la, RHS, 1, tmp, 1);

            cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, tmp, 1);
            resvec[*nbite] = cblas_dnrm2(*la, tmp, 1) / cblas_dnrm2(*la, RHS, 1);
        }

        free(tmp);
    }

    return ret;
}


// Computes and returns MB = D = A+E+F
int extract_MB_jacobi_tridiag ( double *AB, double *MB, int *lab, int *la, int *ku, int*kl, int *kv ) {

    int ret = 0;

    // Error handling
    if      (AB == NULL)  {
        ret = -1;
    }
    else if (MB == NULL) {
        ret = -2;
    }
    else if (*lab != *kv + *ku + *kl + 1 || *lab < 1) {
        ret = -3;
    }
    else if (*la <= 0) {
        ret = -4;
    }
    else if (*ku < 0) {
        ret = -5;
    }
    else if (*kl < 0) {
        ret = -6; 
    }
    else if (*kv < 0) {
        ret = -7; 
    }

    else {

        for(int i = 0; i < *la; i++) {
            for(int j = 0; j < *lab; j++) {

                if(j == *kv + *kl) {
                    MB[i * *lab + j] = AB[i * *lab + j];
                }
                else {
                    MB[i * *lab + j] = 0.0;
                }
            }
        }
    }

    return ret;
}


// Computes and returns MB = D-E = A+F
int extract_MB_gauss_seidel_tridiag ( double *AB, double *MB, int *lab, int *la, int *ku, int*kl, int *kv ) {

    int ret = 0;

    // Error handling
    if      (AB == NULL)  {
        ret = -1;
    }
    else if (MB == NULL) {
        ret = -2;
    }
    else if (*lab != *kv + *ku + *kl + 1 || *lab < 1) {
        ret = -3;
    }
    else if (*la <= 0) {
        ret = -4;
    }
    else if (*ku < 0) {
        ret = -5;
    }
    else if (*kl < 0) {
        ret = -6; 
    }
    else if (*kv < 0) {
        ret = -7; 
    }

    else {

        for(int i = 0; i < *la; i++) {
            for(int j = 0; j < *lab; j++) {

                if(j >= *kv + *kl) {

                    MB[i * *lab + j] = AB[i * *lab + j];

                }
                else {

                    MB[i * *lab + j] = 0.0;
                }

            }
        }
    }
    return ret;
}


// Resolves a linear system in GB with Richardson general method
int richardson_MB ( double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite ) {
    int ret = 0;

    // Error handling
    if      (AB == NULL)  {
        ret = -1;
    }
    else if (RHS == NULL) {
        ret = -2;
    }
    else if (X == NULL)   {
        ret = -3;
    }
    else if (MB == NULL) {
        ret = -4;
    }
    else if (*lab < *ku + *kl + 1 || *lab < 1) {
        ret = -5;
    }
    else if (*la <= 0) {
        ret = -6;
    }
    else if (*ku < 0) {
        ret = -7;
    }
    else if (*kl < 0) {
        ret = -8; 
    }
    else if (resvec == NULL) {
        ret = -11;
    }

    else {

        int *ipiv = (int *) calloc(*la, sizeof(int));
        if(ipiv == NULL) {
            perror("alloc richardson");
            exit(ERR_ALLOC);
        }

        int info = 0;
        const int kuMB = *ku - 1;

        double *MBtr = aligned_alloc(32, *la * *lab *sizeof(double));
        if(MBtr == NULL) {
            perror("alloc richardson");
            exit(ERR_ALLOC);
        }

        // Compute MBtr = L*U = MB 
        cblas_dcopy(*la * *lab, MB, 1, MBtr, 1);

        dgbtrf_(la, la, kl, &kuMB, MBtr, lab, ipiv, &info); 
        if(info != 0) {
            perror("dgbtrf richarson");
            exit(info);
        }

        // MB := M - A in General Band
        cblas_daxpy(*la * *lab, -1.0, AB, 1, MB, 1);

        *nbite = 0;

        double *tmp = aligned_alloc(32, *la * sizeof(double));
        if(tmp == NULL) {
            perror("alloc richardson");
            exit(ERR_ALLOC);
        }

        cblas_dcopy(*la, RHS, 1, tmp, 1);

        // Compute first res
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, tmp, 1);

        resvec[*nbite] = cblas_dnrm2(*la, tmp, 1) / cblas_dnrm2(*la, RHS, 1); 

        int NRHS = 1;

        while(resvec[*nbite] > *tol && *nbite < *maxit) {

            // If the method doesn't converge return an error
            if(*nbite > 1 && resvec[*nbite] >= resvec[*nbite - 1]) {
                *nbite = -1;
                break;
            } 

            // Compute b + (M-A) * x^(i) the RHS of our linear system
            cblas_dcopy(*la, RHS, 1, tmp, 1);

            cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *kl, 1.0, MB, *lab, X, 1, 1.0, tmp, 1);

            // Resolve M x^(i+1) = b + (M-A) * x^(i)
            cblas_dcopy(*la, tmp, 1, X, 1);

            dgbtrs_("N", la, kl, &kuMB, &NRHS, MBtr, lab, ipiv, X, la, &info); 
            if(info != 0) {
                perror("dgbtrf richarson");
                exit(info);
            }

            // Compute the res
            (*nbite)++;

            cblas_dcopy(*la, RHS, 1, tmp, 1);

            cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, tmp, 1);
            resvec[*nbite] = cblas_dnrm2(*la, tmp, 1) / cblas_dnrm2(*la, RHS, 1);
        }

        free(MBtr);
        free(tmp);
        free(ipiv);
    }

    return ret;
}


