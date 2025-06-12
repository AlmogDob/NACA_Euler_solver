/* example input file:
NACA
0012

i_max
50

j_max
25

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

phi
-1

psi
-1
*/
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#define MAXDIR 100
#define MAXWORD 100
#define PI 3.14159265359
#define dprintSTRING(expr) printf(#expr " = %s\n", expr)    /* macro for easy debuging*/
#define dprintINT(expr) printf(#expr " = %d\n", expr)
#define dprintF(expr) printf(#expr " = %g\n", expr)
#define dprintD(expr) printf(#expr " = %g\n", expr)

typedef struct {
    double x;
    double y;
} Vec2;

int create_mesh(double **x_mat, double **y_mat, int NACA, int ni, int nj, int num_points_on_airfoil, double delta_y, double XSF, double YSF, double r, double omega);
void output_solution(char *dir, double *data, int i_max, int j_max);
int offset2d(int i, int j, int ni);
void initialize(double *x_vals_mat, double *y_vals_mat, double *alpha_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat,double *phi_vals_mat);
void set_grid_boundaries(double *x_vals_mat, double *y_vals_mat);
void airfoil(double *x_value, double *y_value, double x, int i);
void interpulat_mat(double *x_vals_mat, char diraction);
double first_deriv(double *mat, char diraction, int i, int j);
double second_deriv(double *mat, char diraction, int i, int j);
void alpha_beta_gama(double *alpha_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *x_vals_mat, double *y_vals_mat);
void psi_phi(double *psi_vals_mat, double *phi_vals_mat, double *x_vals_mat, double *y_vals_mat);
void copy_mat(double *dst, double *src);
void copy_row_to_mat(double *dst, double *src, int row_num);
void copy_col_to_mat(double *dst, double *src, int col_num);
double L_x(double *x_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, int i, int j);
double L_y(double *y_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, int i, int j);
int sweep1(double *fx_vals_mat, double *fy_vals_mat, double *x_vals_mat, double *y_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, double *A, double *B, double *C,  double *D, double *temp_row);
int sweep2(double *Cx_vals_mat, double *Cy_vals_mat, double *fx_vals_mat, double *fy_vals_mat, double *gama_vals_mat, double *A, double *B, double *C, double *D, double *temp_row);
void LHS_sweep1(double *A, double *B, double *C, double *alpha_vals_mat, int j);
void LHS_sweep2(double *A, double *B, double *C, double *alpha_vals_mat, int i);
void RHS_sweep1_x(double *D, double *x_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, int j);
void RHS_sweep2_x(double *D,double *fx_vals_mat, int j);
void RHS_sweep1_y(double *D, double *y_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, int j);
void RHS_sweep2_y(double *D,double *fy_vals_mat, int i);
void BC_sweep1(double *A, double *B, double *C, double *D);
void BC_sweep2(double *A, double *B, double *C, double *D);
int tridiag(double *a, double *b, double *c, double *d, double *u, int is, int ie);
Vec2 step(double *Cx_vals_mat, double *Cy_vals_mat, double *fx_vals_mat, double *fy_vals_mat, double *x_vals_mat_current, double *x_vals_mat_next, double *y_vals_mat_current, double *y_vals_mat_next, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, double *A, double *B, double *C, double *D, double *temp_row);
void mat_print_to_file(FILE *fp, double *data);
void mat_print(double *data);
double calculate_max_L_x(double *x_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat);
double calculate_max_L_y(double *y_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat);