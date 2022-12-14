/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include "atlas_headers.h"

int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr;
  int jj;
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
  double **AAB;
  double *AB;

  double temp, relres;

  NRHS=1;
  nbpoints=102;
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
  dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

  /* LU for tridiagonal matrix  (can replace dgbtrf_) */
  // ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);

  /* write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat"); */
  
  /* Solution (Triangular) */
  if (info==0){
    dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
    if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
  }else{
    printf("\n INFO = %d\n",info);
  }
  // RHS now contains the solution x

  /* It can also be solved with dgbsv */
  // TODO : use dgbsv
  // dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);

  write_xy(X, RHS, &la, "SOL.dat");

  /* Relative forward error */
  // TODO : Compute relative norm of the residual
  cblas_dscal(la, -1.0, RHS, 1);
  cblas_daxpy(la, 1.0, EX_SOL, 1, RHS, 1);
  relres = cblas_dnrm2(la, RHS, 1) / cblas_dnrm2(la, EX_SOL, 1);
  
  printf("\nThe relative forward error for dgbrf+dgbrs = %e\n",relres);

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);

  cblas_dscal(la, -1.0, RHS, 1);
  cblas_daxpy(la, 1.0, EX_SOL, 1, RHS, 1);
  relres = cblas_dnrm2(la, RHS, 1) / cblas_dnrm2(la, EX_SOL, 1);

  printf("\nThe relative forward error for dgbsv = %e\n",relres);

  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  printf("\n\n--------- End -----------\n");
}
