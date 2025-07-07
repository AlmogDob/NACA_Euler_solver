/* example input file:
NACA
0012

ni
51

nj
26

num_points_on_airfoil
31

delta_y
0.02

XSF
1.15

YSF
1.15

r
0.001

omega
0.1
*/
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>

#ifndef MAXDIR
    #define MAXDIR 1000
#endif
#ifndef MAXWORD
    #define MAXWORD 2000
#endif
#ifndef PI
    #ifndef __USE_MISC
    #define __USE_MISC
    #endif
    #include <math.h>
    #define PI M_PI
#endif
#ifndef dprintSTRING
    #define dprintSTRING(expr) printf(#expr " = %s\n", expr)
#endif
#ifndef dprintINT
    #define dprintINT(expr) printf(#expr " = %d\n", expr)
#endif
#ifndef dprintF
    #define dprintF(expr) printf(#expr " = %g\n", expr)
#endif
#ifndef dprintD
    #define dprintD(expr) printf(#expr " = %g\n", expr)
#endif

typedef struct {
    double x;
    double y;
} Vec2;

int create_mesh(double **x_mat, double **y_mat, int NACA, int ni, int nj, int num_points_on_airfoil, double delta_y, double XSF, double YSF, double r, double omega, char *output_dir);
void output_solution(char *file, double *data, int i_max, int j_max);
// int offset2d(int i, int j, int ni);
#define offset2d(i, j, ni) (j) * (ni) + (i)
void initialize(double *x_vals_mat, double *y_vals_mat, double *alpha_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, double *phi_vals_mat, int i_TEL, int i_TEU, int i_LE, double delta_x, double delta_y, double XSF, double YSF, int i_max, int j_max, int i_min, int j_min, int NACA, double phi_valuse, double psi_valuse);
void set_grid_boundaries(double *x_vals_mat, double *y_vals_mat, int i_TEL, int i_TEU, int i_LE, double delta_x, double delta_y, double XSF, double YSF, int i_max, int j_max, int NACA);
void airfoil(double *x_value, double *y_value, double x, int i, int NACA, int i_LE);
void interpulat_mat(double *mat, char diraction, int i_max, int j_max, int i_min, int j_min);
double first_deriv(double *mat, char diraction, int i_min, int i_max, int j_min, int j_max, int i, int j);
double second_deriv(double *mat, char diraction, int i_min, int i_max, int j_min, int j_max, int i, int j);
void alpha_beta_gama(double *alpha_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *x_vals_mat, double *y_vals_mat, int i_max, int i_min, int j_max, int j_min);
void psi_phi(double *psi_vals_mat, double *phi_vals_mat, double *x_vals_mat, double *y_vals_mat, int i_min, int i_max, int j_min, int j_max, double phi_valuse, double psi_valuse);
void copy_mat(double *dst, double *src, int i_max, int j_max);
void copy_row_to_mat(double *dst, double *src, int row_num, int i_max);
void copy_col_to_mat(double *dst, double *src, int col_num, int i_max, int j_max);
double L_x(double *x_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, int i_max, int j_max, int i, int j);
double L_y(double *y_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, int i_max, int j_max, int i, int j);
int sweep1(double *fx_vals_mat, double *fy_vals_mat, double *x_vals_mat, double *y_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, double *A, double *B, double *C,  double *D, double *temp_row, int i_max, int j_max, double r, double omega);
int sweep2(double *Cx_vals_mat, double *Cy_vals_mat, double *fx_vals_mat, double *fy_vals_mat, double *gama_vals_mat, double *A, double *B, double *C, double *D, double *temp_row, int i_max, int j_max, double r);
void LHS_sweep1(double *A, double *B, double *C, double *alpha_vals_mat, int j, int i_max, double r);
void LHS_sweep2(double *A, double *B, double *C, double *gama_vals_mat, int i, int i_max, int j_max, double r);
void RHS_sweep1_x(double *D, double *x_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, int j, int i_max, int j_max, double r, double omega);
void RHS_sweep2_x(double *D,double *fx_vals_mat, int i, int i_max, int j_max);
void RHS_sweep1_y(double *D, double *y_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, int j, int i_max, int j_max, double r, double omega);
void RHS_sweep2_y(double *D,double *fy_vals_mat, int i, int i_max, int j_max);
void BC_sweep1(double *A, double *B, double *C, double *D, int i_max);
void BC_sweep2(double *A, double *B, double *C, double *D, int j_max);
int tridiag(double *a, double *b, double *c, double *d, double *u, int is, int ie);
Vec2 step(double *Cx_vals_mat, double *Cy_vals_mat, double *fx_vals_mat, double *fy_vals_mat, double *x_vals_mat_current, double *x_vals_mat_next, double *y_vals_mat_current, double *y_vals_mat_next, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, double *A, double *B, double *C, double *D, double *temp_row, int i_max, int i_min, int j_max, int j_min, double r, double omega);
void mat_print_to_file(FILE *fp, double *data, int i_max, int j_max);
void mat_print(double *data, int i_max, int j_max);
double calculate_max_L_x(double *x_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, int i_max, int j_max);
double calculate_max_L_y(double *y_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, int i_max, int j_max);