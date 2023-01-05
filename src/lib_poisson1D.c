/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){


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

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){


  for(int i = 0; i < *la; i++) {
    for(int j = 0; j < *lab; j++) {
      if(j == 1) {
        AB[i * *lab + j] = 1;
      }
      else {
        AB[i * *lab + j] = 0;
      }
    }
  }
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  // On a posé g(i) = 0
  RHS[0] = *BC0;
  for(int i = 1; i < *la - 1; i++) {
    RHS[i] = 0;
  }
  RHS[*la - 1] = *BC1;
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  // La solution analytique est T(x) = T0 + x(T1-T0) si il n'y a pas de source de chaleur comme considéré
  for(int i = 0; i < *la; i++) {

    EX_SOL[i] = *BC0 + X[i] * (*BC1 - *BC0);

  }
}  

void set_grid_points_1D(double* X, int* la){
  // X(i) = i/la = ih avec h un pas constant
  
  const double h = 1.0/ (double) (*la + 1); // On a n-1 intervalles entre n points et la = npoints - 2

  for(int i = 1; i < *la + 1; i++) {

    X[i-1] = i*h;

  }
  // Pour tester, on doit avoir pour tout i 0 < X(i) < 1 (strictement car on connait déjà les limites)
  // Plus spécifiquement, X(0) - h = 0 et X(la) + h = 1
}

void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
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

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*lab);ii++){
      for (jj=0;jj<(*la);jj++){
	fprintf(file,"%02.6lf\t",AB[jj*(*lab)+ii]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename){
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

void write_vec(double* vec, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%1.6lf\n",vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void write_xy(double* vec, double* x, int* la, char* filename){
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
  return 0;
}
int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){

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
  

void eig_poisson1D(double* eigval, int *la){
}

double eigmax_poisson1D(int *la){
  return 0;
}

double eigmin_poisson1D(int *la){
  return 0;
}

double richardson_alpha_opt(int *la){
  return 0;
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){

}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){

}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){

}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){

}


