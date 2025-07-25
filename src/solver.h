#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>

#ifndef MAXDIR
    #define MAXDIR 1000
#endif
#ifndef MAXWORD
    #define MAXWORD 3000
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

int solver(const char *output_dir, double **rho_2Dmat, double **u_2Dmat, double **v_2Dmat, double **e_2Dmat, double *x_vals_mat, double *y_vals_mat, double *final_S_norm, int ni, int nj, int num_points_on_airfoil, const double Mach_inf, const double angle_of_attack_deg, const double density, const double environment_pressure, const double delta_t, const double Gamma, const double epse, const double max_iteration);
void read_mat_from_file(FILE *fp, double *des, int ni, int nj);
void output_solution_solver(const char *output_dir, double *current_Q, double *U_mat, double *V_mat, double *x_vals_mat, double *y_vals_mat, double *dxi_dx_mat, double *dxi_dy_mat, double *deta_dx_mat, double *deta_dy_mat, int ni, int nj, int i_TEL, int i_LE, int i_TEU);
// int offset2d_solver(int i, int j, int ni, int nj);
#define offset2d_solver(i, j, ni, nj)  (j) * (ni) + (i)
// int offset3d(int i, int j, int k, int ni, int nj);
#define offset3d(i, j, k, ni, nj) ((k) * (nj) + (j)) * (ni) + (i)
void print_mat2D(double *data, int ni, int nj);
void print_layer_of_mat3D(double *data, int layer, int ni, int nj);
double first_deriv_solver(double *mat, char diraction, int i, int j, int ni, int nj);
double calculate_one_over_jacobian_at_a_point(double *x_vals_mat, double *y_vals_mat, int i, int j, int ni, int nj);
void contravariant_velocities(double *U, double *V, double *dxi_dx_mat, double *dxi_dy_mat, double *deta_dx_mat, double *deta_dy_mat, double *Q, int i, int j, int ni, int nj);
void calculate_u_and_v(double *u, double *v, double *Q, int i, int j, int ni, int nj);
double calculate_p(double energy, double rho, double u, double v, const double Gamma);
double calculate_energy(double p, double u, double v, double rho, const double Gamma);
void calculate_E_hat_at_a_point(double *E0, double *E1, double *E2, double *E3, double *J_vals_mat, double *dxi_dx_mat, double *dxi_dy_mat, double *deta_dx_mat, double *deta_dy_mat, double *Q, int i, int j, int ni, int nj, const double Gamma);
void calculate_F_hat_at_a_point(double *F0, double *F1, double *F2, double *F3, double *J_vals_mat, double *dxi_dx_mat, double *dxi_dy_mat, double *deta_dx_mat, double *deta_dy_mat, double *Q, int i, int j, int ni, int nj, const double Gamma);
void initialize_flow_field(double *Q, int ni, int nj, const double Mach_inf, const double angle_of_attack_rad, const double environment_pressure, const double density, const double Gamma);
void matrices_coeffic_and_Jacobian(double *J_vals_mat, double *dxi_dx_mat, double *dxi_dy_mat, double *deta_dx_mat, double *deta_dy_mat, double *x_vals_mat, double *y_vals_mat, int ni, int nj);
void initialize_solver(double *Q, double *J_vals_mat, double *dxi_dx_mat, double *dxi_dy_mat, double *deta_dx_mat, double *deta_dy_mat, double *x_vals_mat, double *y_vals_mat, int ni, int nj, const double Mach_inf, const double angle_of_attack_rad, const double environment_pressure, const double density, const double Gamma);
void RHS(double *S, double *W, double *Q, double *J_vals_mat, double *dxi_dx_mat, double *dxi_dy_mat, double *deta_dx_mat, double *deta_dy_mat, double *s2, double *rspec, double *qv, double *dd, int ni, int nj, int max_ni_nj, const double Mach_inf, const double delta_t, const double Gamma, const double epse);
void advance_Q(double *next_Q, double *current_Q ,double *S, double *J_vals_mat, int ni, int nj);
void copy_3Dmat_to_3Dmat(double *dst, double *src, int ni, int nj);
int smooth(double *q, double *s, double *jac, double *xx, double *xy, double *yx, double *yy, int id, int jd, double *s2, double *rspec, double *qv, double *dd, double epse, double gamma, double fsmach, double dt);
void apply_BC(double *Q, double *J_vals_mat, double *dxi_dx_mat, double *dxi_dy_mat, double *deta_dx_mat, double *deta_dy_mat, int ni, int nj, int i_TEL, int i_LE, int i_TEU, const double Gamma);
void output_mat2D_to_file(FILE *fp, double *data, int ni, int nj);
void output_layer_of_mat3D_to_file(FILE *fp, double *data, int layer, int ni, int nj);
void calculate_A_hat_j_const(double *dst, double *Q, double *dxi_dx_mat, double *dxi_dy_mat, int i, int j, int ni, int nj, const double Gamma);
void calculate_B_hat_i_const(double *dst, double *Q, double *deta_dx_mat, double *deta_dy_mat, int i, int j, int ni, int nj, const double Gamma);
int smoothx(double *q, double *xx, double *xy, int id, int jd, double *a, double *b, double *c, int j,double *jac, double *drr, double *drp, double *rspec, double *qv, double *dd, double epsi, double gamma, double fsmach, double dt);
int smoothy(double *q, double *yx, double *yy, int id, int jd, double *a, double *b, double *c, int i,double *jac, double *drr, double *drp, double *rspec, double *qv, double *dd, double epsi, double gamma, double fsmach, double dt);
void LHSX(double *A, double *B, double *C, double *Q, double *dxi_dx_mat, double *dxi_dy_mat, int j, int ni, int nj, const double Gamma, const double delta_t);
void LHSY(double *A, double *B, double *C, double *Q, double *deta_dx_mat, double *deta_dy_mat, int i, int ni, int nj, int max_ni_nj , const double Gamma, const double delta_t);
int btri4s(double *a, double *b, double *c, double *f, int kd, int ks, int ke);
double calculate_S_norm(double *S, int ni, int nj);
double step_solver(double *A, double *B, double *C, double *D, double *current_Q, double *S, double *W, double *J_vals_mat, double *dxi_dx_mat, double *dxi_dy_mat, double *deta_dx_mat, double *deta_dy_mat, double *s2, double *drr, double *drp, double *rspec, double *qv, double *dd, int ni, int nj, int max_ni_nj, const double Mach_inf, const double delta_t, const double Gamma, const double epse, const double epsi);
double fix_delta_t(double *current_Q, double *x_vals_mat, double *y_vals_mat, int ni, int nj);
