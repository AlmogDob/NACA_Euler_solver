#include "mesher.h"

/* return 0 on success */
int create_mesh(double **x_mat, double **y_mat, int NACA, int ni, int nj, int num_points_on_airfoil, double delta_y, double XSF, double YSF, double r, double omega, char *output_dir)
{
    double psi_valuse = -1, phi_valuse = -1;
    int i_max = ni-1;
    int j_max = nj-1;
    int i_min = 0;
    int j_min = 0;

    /* declarations */
    char temp_word[MAXWORD];
    int i_index, j_index;
    double *x_vals_mat_init, *y_vals_mat_init, *x_vals_mat_current, *y_vals_mat_current, *x_vals_mat_next, *y_vals_mat_next, *alpha_vals_mat, *beta_vals_mat, *gama_vals_mat, *psi_vals_mat, *phi_vals_mat, *fx_vals_mat, *fy_vals_mat, *Cx_vals_mat, *Cy_vals_mat;
    /* matrix diaganosl for different sweeps */
    double *A, *B, *C, *D, *temp_row;
    Vec2 result, first_result;


    strcpy(temp_word, output_dir);
    strcat(temp_word, "/Ls_valuse.txt");
    FILE *Ls_fp = fopen(temp_word, "wt");
    if (!Ls_fp) {
        fprintf(stderr, "%s:%d: [ERROR] unable to open file '%s'", __FILE__, __LINE__, temp_word);
        exit(1);
    }

    /*------------------------------------------------------------*/

    int i_LE = i_max / 2;
    int i_TEL = i_LE - num_points_on_airfoil / 2 + 1;
    int i_TEU = i_LE + num_points_on_airfoil / 2 - 1;

    double delta_x = 1.0/(i_LE - i_TEL);

    // dprintINT(i_LE);
    // dprintINT(i_TEL);
    // dprintINT(i_TEU);


    /*------------------------------------------------------------*/
    
    /* Memory allocation */
    x_vals_mat_init = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            x_vals_mat_init[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    y_vals_mat_init = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            y_vals_mat_init[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    x_vals_mat_current = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            x_vals_mat_current[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    y_vals_mat_current = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            y_vals_mat_current[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    x_vals_mat_next = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            x_vals_mat_next[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    y_vals_mat_next = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            y_vals_mat_next[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    alpha_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            alpha_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    beta_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            beta_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    gama_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            gama_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    psi_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            psi_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    phi_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            phi_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    fx_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            fx_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    fy_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            fy_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    Cx_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            Cx_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }
    Cy_vals_mat = (double *)malloc(sizeof(double) * (i_max + 1) * (j_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the matrix with zeros */
        for (j_index = 0; j_index < j_max+1; j_index++) {
            Cy_vals_mat[offset2d(i_index, j_index, i_max+1)] = 0;
        }
    }

    A = (double *)malloc(sizeof(double) * (i_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the array with zeros */
        A[i_index] = 0;
    }
    B = (double *)malloc(sizeof(double) * (i_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the array with zeros */
        B[i_index] = 0;
    }
    C = (double *)malloc(sizeof(double) * (i_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the array with zeros */
        C[i_index] = 0;
    }
    D = (double *)malloc(sizeof(double) * (i_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the array with zeros */
        D[i_index] = 0;
    }
    temp_row = (double *)malloc(sizeof(double) * (i_max + 1));
    for (i_index = 0; i_index < i_max+1; i_index++) {   /* filling the array with zeros */
        temp_row[i_index] = 0;
    }
    
    /*------------------------------------------------------------*/
    
    initialize(x_vals_mat_init, y_vals_mat_init, alpha_vals_mat, beta_vals_mat, gama_vals_mat, psi_vals_mat, phi_vals_mat, i_TEL, i_TEU, i_LE, delta_x, delta_y, XSF, YSF, i_max, j_max, i_min, j_min, NACA, phi_valuse, psi_valuse);
    
    copy_mat(x_vals_mat_current, x_vals_mat_init, i_max, j_max);
    copy_mat(x_vals_mat_next, x_vals_mat_init, i_max, j_max);
    copy_mat(y_vals_mat_current, y_vals_mat_init, i_max, j_max);
    copy_mat(y_vals_mat_next, y_vals_mat_init, i_max, j_max);

    for (i_index = 0; i_index < 1e5; i_index++) {
        result = step(Cx_vals_mat, Cy_vals_mat, fx_vals_mat, fy_vals_mat, x_vals_mat_current, x_vals_mat_next, y_vals_mat_current, y_vals_mat_next, alpha_vals_mat, phi_vals_mat, beta_vals_mat, gama_vals_mat, psi_vals_mat, A, B, C, D, temp_row, i_max, i_min, j_max, j_min, r, omega);
        if (i_index == 0) {
            first_result = result;
        }
        if (result.x == 1 && result.y == 0) {
            fprintf(stderr, "ERROR: Step - sweep 1\n");
            exit(1);
        }
        if (result.x == 2 && result.y == 0) {
            fprintf(stderr, "ERROR: Step - sweep 2\n");
            exit(2);
        }

        copy_mat(x_vals_mat_current, x_vals_mat_next, i_max, j_max);
        copy_mat(y_vals_mat_current, y_vals_mat_next, i_max, j_max);
        
        /* printing Lx and Ly */
        if (!((i_index+1) % 1000)) {
            printf("%4d. Lx_max: %0.10f, Ly_max: %0.10f\n",i_index+1, result.x, result.y);
        }
        fprintf(Ls_fp, "%g, %g\n", result.x, result.y);

        /* checking convergenc */
        if (log10(fabs(first_result.x/result.x)) > 5 && log10(fabs(first_result.y/result.y))) {
            break;
        }
    }

    strcpy(temp_word, output_dir);
    strcat(temp_word, "/x_mat_init.txt");
    output_solution_mesher(temp_word, x_vals_mat_init, i_max, j_max);

    strcpy(temp_word, output_dir);
    strcat(temp_word, "/y_mat_init.txt");
    output_solution_mesher(temp_word, y_vals_mat_init, i_max, j_max);

    strcpy(temp_word, output_dir);
    strcat(temp_word, "/x_mat.txt");
    output_solution_mesher(temp_word, x_vals_mat_next, i_max, j_max);

    strcpy(temp_word, output_dir);
    strcat(temp_word, "/y_mat.txt");
    output_solution_mesher(temp_word, y_vals_mat_next, i_max, j_max);

    *x_mat = x_vals_mat_next;
    *y_mat = y_vals_mat_next;

    /*------------------------------------------------------------*/
    /* freeing the memory */
    free(x_vals_mat_init);
    free(y_vals_mat_init);
    free(x_vals_mat_current);
    free(y_vals_mat_current);
    // free(x_vals_mat_next);
    // free(y_vals_mat_next);
    free(alpha_vals_mat);
    free(beta_vals_mat);
    free(gama_vals_mat);
    free(psi_vals_mat);
    free(phi_vals_mat);
    free(fx_vals_mat);
    free(fy_vals_mat);
    free(Cx_vals_mat);
    free(Cy_vals_mat);
    free(A);
    free(B);
    free(C);
    free(D);
    free(temp_row);

    return 0;
}

/* output data;
argument list:
dir - the directory of the output file.
data - the solution vector */
void output_solution_mesher(char *file, double *data, int i_max, int j_max)
{
    FILE *fp = fopen(file, "wt");
    if (!fp) {
        fprintf(stderr, "%s:%d: [ERROR] unable to open file '%s'", __FILE__, __LINE__, file);
        exit(1);
    }
    int i, j;
    
    for (j = 0; j < j_max+1; j++) {
        for (i = 0; i < i_max+1; i++) {
            fprintf(fp, "%g ", data[offset2d(i, j, i_max+1)]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

/* converts a 2D index into 1D index
argument list:
i - x position
j - y position
ni - stride */
// int offset2d(int i, int j, int ni)
// {
//     return j * ni + i;
// }

/* set inital valsuse of the mash points 
argument list:
x_vals_mat - 1D array of the x valuse 
y_vals_mat - 1D array of the y valuse
alpha_vals_mat - 1D array for the alpha valus
beta_vals_mat - 1D array for the beta valus
gama_vals_mat - 1D array for the gama valus */
void initialize(double *x_vals_mat, double *y_vals_mat, double *alpha_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, double *phi_vals_mat, int i_TEL, int i_TEU, int i_LE, double delta_x, double delta_y, double XSF, double YSF, int i_max, int j_max, int i_min, int j_min, int NACA, double phi_valuse, double psi_valuse)
{
    set_grid_boundaries(x_vals_mat, y_vals_mat, i_TEL, i_TEU, i_LE, delta_x, delta_y, XSF, YSF, i_max, j_max, NACA);
    interpulat_mat(x_vals_mat, 'j', i_max, j_max, i_min, j_min);
    interpulat_mat(y_vals_mat, 'j', i_max, j_max, i_min, j_min);
    alpha_beta_gama(alpha_vals_mat, beta_vals_mat, gama_vals_mat, x_vals_mat, y_vals_mat, i_max, i_min, j_max, j_min);
    psi_phi(psi_vals_mat, phi_vals_mat, x_vals_mat, y_vals_mat, i_min, i_max, j_min, j_max, phi_valuse, psi_valuse);
}

/* set the mash boundaries coorditates 
argument list: x_vals_mat - 1D array of the x valuse 
y_vals_mat - 1D array of the y valuse */
void set_grid_boundaries(double *x_vals_mat, double *y_vals_mat, int i_TEL, int i_TEU, int i_LE, double delta_x, double delta_y, double XSF, double YSF, int i_max, int j_max, int NACA)
{
    int i_index, i_min = 0, j_index, j_min = 0, index = 0, num_points_befor_circle,
    num_of_outer_segments, num_of_top_outer_segments;

    double x, x_temp, x_i_minos_2, x_i_minos_1, y, y_j_minos_2, y_j_minos_1, y_imax_jmax, x_imax_jmax,
    delta_theta, R, theta = 0, length_top_outer, segment_length, current_x_val = 0;

    /* airfoil */
    for (i_index = i_TEL, j_index = j_min; i_index < i_TEU+1; i_index++) { 
        x_temp = 1 - cos(0.5*PI*(i_LE-i_index)*delta_x);
        airfoil(&x, &y, x_temp, i_index, NACA, i_LE);
        x_vals_mat[offset2d(i_index, j_index, i_max+1)] = x;
        y_vals_mat[offset2d(i_index, j_index, i_max+1)] = y;
    }
    /* Wake */
    for (i_index = i_TEU + 1, j_index = j_min; i_index < i_max+1; i_index++) {
        x_i_minos_1 = x_vals_mat[offset2d(i_index-1, j_index, i_max+1)];  
        x_i_minos_2 = x_vals_mat[offset2d(i_index-2, j_index, i_max+1)];  
        x_vals_mat[offset2d(i_index, j_index, i_max+1)] = x_i_minos_1 + (x_i_minos_1 - x_i_minos_2) * XSF;
    }
    for (i_index = i_min, j_index = j_min; i_index < i_TEL; i_index++) {
        x_vals_mat[offset2d(i_index, j_index, i_max+1)] = x_vals_mat[offset2d(i_max-i_index, j_index, i_max+1)];
    }
    /* Exit Boundary */
    y_vals_mat[offset2d(i_max, j_min+1, i_max+1)] = delta_y;
    for (i_index = i_max, j_index = j_min+2; j_index < j_max+1; j_index++) {
        y_j_minos_1 = y_vals_mat[offset2d(i_index, j_index-1, i_max+1)];  
        y_j_minos_2 = y_vals_mat[offset2d(i_index, j_index-2, i_max+1)];  
        y_vals_mat[offset2d(i_max, j_index, i_max+1)] = y_j_minos_1 + (y_j_minos_1 - y_j_minos_2) * YSF;
    }
    for (i_index = i_max, j_index = j_min+1; j_index < j_max+1; j_index++) {
        x_vals_mat[offset2d(i_max, j_index, i_max+1)] = x_vals_mat[offset2d(i_max, j_min, i_max+1)];
        y_vals_mat[offset2d(i_min, j_index, i_max+1)] = -y_vals_mat[offset2d(i_max, j_index, i_max+1)];
        x_vals_mat[offset2d(i_min, j_index, i_max+1)] = x_vals_mat[offset2d(i_min, j_min, i_max+1)];
    }
    /* Outer boundary */
    y_imax_jmax = y_vals_mat[offset2d(i_max, j_max, i_max+1)];
    x_imax_jmax = x_vals_mat[offset2d(i_max, j_max, i_max+1)];
    R = y_imax_jmax;
    /*test*/
    // dprintD(x_imax_jmax);
    /*test*/

    num_of_outer_segments = i_max;
    num_of_top_outer_segments = num_of_outer_segments/2;
    length_top_outer = x_imax_jmax + 0.5*PI*R; /* length of stright part and quarter of the circle */
    segment_length = length_top_outer/num_of_top_outer_segments;

    /* the stright line part */
    for (num_points_befor_circle = 0;
         num_points_befor_circle < num_of_top_outer_segments + 1;
         num_points_befor_circle++) {
            current_x_val = x_imax_jmax - num_points_befor_circle*segment_length; 
            if (current_x_val < 0) {
                break;
            }
            x_vals_mat[offset2d(i_max-num_points_befor_circle, j_max, i_max+1)] = current_x_val;
            y_vals_mat[offset2d(i_max-num_points_befor_circle, j_max, i_max+1)] = y_imax_jmax;
    }

    theta = PI/2 + atan(x_vals_mat[offset2d(i_max-num_points_befor_circle+1, j_max, i_max+1)] / R);
    delta_theta = theta / (num_of_top_outer_segments - num_points_befor_circle + 1);

    /* the quarter circle part */
    for (index = 0; index < num_of_top_outer_segments - num_points_befor_circle + 1; index++) {
        x_vals_mat[offset2d(i_max/2 + index, j_max, i_max+1)] = -R*cos(delta_theta*index);
        y_vals_mat[offset2d(i_max/2 + index, j_max, i_max+1)] = R*sin(delta_theta*index);
    }

    /* coping to the bottom side */
    for (index = 1; index < i_max/2 + 1; index++) {
        x_vals_mat[offset2d(i_max/2 - index, j_max, i_max+1)] = x_vals_mat[offset2d(i_max/2 + index, j_max, i_max+1)];
        y_vals_mat[offset2d(i_max/2 - index, j_max, i_max+1)] = -y_vals_mat[offset2d(i_max/2 + index, j_max, i_max+1)];
    }
}

/* returns the shape of the airfoil as a function of x
argument list: 
i - index along the airfoil */
void airfoil(double *x_value, double *y_value, double x, int i, int NACA, int i_LE)
{
    double m = ((NACA % 10000)/1000) / (double)100;
    double p = ((NACA % 1000)/100) / (double)10;
    double t = (NACA % 100) / (double)100;

    double y_t = 5 * t * (0.2969 * sqrt(x) - 0.1260 * x - 0.3516 * x * x + 0.2843 * x * x * x - 0.1036 * x * x * x * x);
    double y_c, dy_c__dx;
    if (x <= p) {
        y_c      = m / p / p * (2 * p * x - x * x);
        dy_c__dx = m / p / p * (p - x);
    } else {
        y_c      = m / (1 - p) / (1 - p) * ((1 - 2 * p) + 2 * p * x - x * x);
        dy_c__dx = 2 * m / (1 - p) / (1 - p) * (p - x);
    }
    double theta = atan(dy_c__dx);
    double x_U = x - y_t * sin(theta);
    double y_U = y_c + y_t * cos(theta);
    double x_L = x + y_t * sin(theta);
    double y_L = y_c - y_t * cos(theta);
    double x_C = x;
    double y_C = y_c;

    if (i <= i_LE) {
        *x_value = x_L;
        *y_value = y_L;
    } else if (i > i_LE) {
        *x_value = x_U;
        *y_value = y_U;
    // } else {
    //     *x_value = 0; 
    //     *y_value = 0;
    }
    (void)y_C;
    (void)x_C;
}

/* fill the matrix in a way of interpulation the boundarys
argument list:
mat - 1D array of valsuse
diraction - i direction or j direction */
void interpulat_mat(double *mat, char diraction, int i_max, int j_max, int i_min, int j_min)
{
    int i, j; 
    double max, min;

    assert(diraction == 'j' || diraction == 'i'); /* different directions are not implemented */

    if (diraction == 'j') {
        for (i = 1; i < i_max; i++) {
            for (j = 1; j < j_max; j++) {
                max = mat[offset2d(i, j_max, i_max+1)];
                min = mat[offset2d(i, j_min, i_max+1)];
                mat[offset2d(i, j, i_max+1)] = (max - min)/(j_max) * j + min;   /* liniar interpulation */
            }
        }
    }
    if (diraction == 'i') {
        for (j = 1; j < j_max; j++) {
            for (i = 1; i < i_max; i++) {
                max = mat[offset2d(i_max, j, i_max+1)];
                min = mat[offset2d(i_min, j, i_max+1)];
                mat[offset2d(i, j, i_max+1)] = (max - min)/(i_max) * i + min;   /* liniar interpulation */
            }
        }
    }
}

/* return the second order first derivitive from the valuse in the mat matrix
argument list:
mat - 1D array of valuse
diraction - i or j
i, j - the points coordinates */
double first_deriv(double *mat, char diraction, int i_min, int i_max, int j_min, int j_max, int i, int j)
{
    if (diraction == 'j') {
        if (j == j_min || j == j_max) {
            return 0;
        }
        return (mat[offset2d(i, j+1, i_max+1)] - mat[offset2d(i, j-1, i_max+1)]) / (2); /* second order first derivitive */
    }
    if (diraction == 'i') {
        if (i == i_min || i == i_max) {
            return 0;
        }
        return (mat[offset2d(i+1, j, i_max+1)] - mat[offset2d(i-1, j, i_max+1)]) / (2); /* second order first derivitive */
    }
    return NAN;
}

/* return the second order second derivitive from the valuse in the mat matrix
argument list:
mat - 1D array of valuse
diraction - i or j
i, j - the points coordinates */
double second_deriv(double *mat, char diraction, int i_min, int i_max, int j_min, int j_max, int i, int j)
{
    if (diraction == 'j') {
        if (j == j_min || j == j_max) {
            return 0;
        }
        return (mat[offset2d(i, j+1, i_max+1)] -2*mat[offset2d(i, j, i_max+1)] + mat[offset2d(i, j-1, i_max+1)]) / (1); /* second order second derivitive */
    }
    if (diraction == 'i') {
        if (i == i_min || i == i_max) {
            return 0;
        }
        return (mat[offset2d(i+1, j, i_max+1)] -2*mat[offset2d(i, j, i_max+1)] + mat[offset2d(i-1, j, i_max+1)]) / (1); /* second order second derivitive */
    }
    return NAN;
}

/* fills the alpha and beta and gama matrices
argument list:
alpha_vals_mat - 1D array for the alpha valus
beta_vals_mat - 1D array for the beta valus
gama_vals_mat - 1D array for the gama valus
x_vals_mat - 1D array of the x valus
y_vals_mat - 1D array of the y valus */
void alpha_beta_gama(double *alpha_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *x_vals_mat, double *y_vals_mat, int i_max, int i_min, int j_max, int j_min)
{
    /* according to equation 3 */
    int i, j;
    double Dx_Deta, Dy_Deta, Dx_Dxai, Dy_Dxai;

    for (i = 0; i < i_max + 1; i++) {
        for (j = 0; j < j_max + 1; j++) {
            Dx_Deta = first_deriv(x_vals_mat, 'j', i_min, i_max, j_min, j_max, i, j);
            Dy_Deta = first_deriv(y_vals_mat, 'j', i_min, i_max, j_min, j_max, i, j);
            Dx_Dxai = first_deriv(x_vals_mat, 'i', i_min, i_max, j_min, j_max, i, j);
            Dy_Dxai = first_deriv(y_vals_mat, 'i', i_min, i_max, j_min, j_max, i, j);
            alpha_vals_mat[offset2d(i, j, i_max+1)] = Dx_Deta*Dx_Deta + Dy_Deta*Dy_Deta;
            beta_vals_mat[offset2d(i, j, i_max+1)] = Dx_Dxai*Dx_Deta + Dy_Dxai*Dy_Deta;
            gama_vals_mat[offset2d(i, j, i_max+1)] = Dx_Dxai*Dx_Dxai + Dy_Dxai*Dy_Dxai;
        }
    }
}

/* creat and fill the psi and phi matrices
argument list:
psi_vals_mat - 1D array for the psi valus
phi_vals_mat - 1D array for the phi valus
x_vals_mat - 1D array of the x valus
y_vals_mat - 1D array of the y valus */
void psi_phi(double *psi_vals_mat, double *phi_vals_mat, double *x_vals_mat, double *y_vals_mat, int i_min, int i_max, int j_min, int j_max, double phi_valuse, double psi_valuse)
{
    /* according to equation 4 and 5 */
    int i, j;
    double Dx_Deta_min, Dy_Deta_min, Dx_Dxai_min, Dy_Dxai_min, Dx_Deta_max, Dy_Deta_max, Dx_Dxai_max, Dy_Dxai_max,
    Dx_Deta_Deta_min, Dx_Deta_Deta_max, Dy_Deta_Deta_min, Dy_Deta_Deta_max, Dx_Dxai_Dxai_min, Dx_Dxai_Dxai_max,
    Dy_Dxai_Dxai_min, Dy_Dxai_Dxai_max;

    /* eq 4 */
    for (j = 0; j < j_max+1; j++) {
        if (psi_valuse == -1) {
            Dx_Deta_min = first_deriv(x_vals_mat, 'j', i_min, i_max, j_min, j_max, i_min, j);
            Dy_Deta_min = first_deriv(y_vals_mat, 'j', i_min, i_max, j_min, j_max, i_min, j);
            Dx_Deta_max = first_deriv(x_vals_mat, 'j', i_min, i_max, j_min, j_max, i_max, j);
            Dy_Deta_max = first_deriv(y_vals_mat, 'j', i_min, i_max, j_min, j_max, i_max, j);
            Dx_Deta_Deta_min = second_deriv(x_vals_mat, 'j', i_min, i_max, j_min, j_max, i_min, j);
            Dy_Deta_Deta_min = second_deriv(y_vals_mat, 'j', i_min, i_max, j_min, j_max, i_min, j);
            Dx_Deta_Deta_max = second_deriv(x_vals_mat, 'j', i_min, i_max, j_min, j_max, i_max, j);
            Dy_Deta_Deta_max = second_deriv(y_vals_mat, 'j', i_min, i_max, j_min, j_max, i_max, j);

            if (fabs(Dy_Deta_min) > fabs(Dx_Deta_min)) {
                psi_vals_mat[offset2d(i_min, j, i_max+1)] = - Dy_Deta_Deta_min / Dy_Deta_min;
            }
            if (fabs(Dy_Deta_min) < fabs(Dx_Deta_min)) {
                psi_vals_mat[offset2d(i_min, j, i_max+1)] = - Dx_Deta_Deta_min / Dx_Deta_min;
            }
            if (fabs(Dy_Deta_max) > fabs(Dx_Deta_max)) {
                psi_vals_mat[offset2d(i_max, j, i_max+1)] = - Dy_Deta_Deta_max / Dy_Deta_max;
            }
            if (fabs(Dy_Deta_max) < fabs(Dx_Deta_max)) {
                psi_vals_mat[offset2d(i_max, j, i_max+1)] = - Dx_Deta_Deta_max / Dx_Deta_max;
            }
        } else {
            psi_vals_mat[offset2d(i_min, j, i_max+1)] = psi_valuse;
            psi_vals_mat[offset2d(i_max, j, i_max+1)] = psi_valuse;
        }
    }

    /* eq 5 */
    for (i = 0; i < i_max+1; i++) {
        if (phi_valuse == -1) {
            Dx_Dxai_min = first_deriv(x_vals_mat, 'i', i_min, i_max, j_min, j_max, i, j_min);
            Dy_Dxai_min = first_deriv(y_vals_mat, 'i', i_min, i_max, j_min, j_max, i, j_min);
            Dx_Dxai_max = first_deriv(x_vals_mat, 'i', i_min, i_max, j_min, j_max, i, j_max);
            Dy_Dxai_max = first_deriv(y_vals_mat, 'i', i_min, i_max, j_min, j_max, i, j_max);
            Dx_Dxai_Dxai_min = second_deriv(x_vals_mat, 'i', i_min, i_max, j_min, j_max, i, j_min);
            Dy_Dxai_Dxai_min = second_deriv(y_vals_mat, 'i', i_min, i_max, j_min, j_max, i, j_min);
            Dx_Dxai_Dxai_max = second_deriv(x_vals_mat, 'i', i_min, i_max, j_min, j_max, i, j_max);
            Dy_Dxai_Dxai_max = second_deriv(y_vals_mat, 'i', i_min, i_max, j_min, j_max, i, j_max);

            if (fabs(Dx_Dxai_min) > fabs(Dy_Dxai_min)) {
                phi_vals_mat[offset2d(i, j_min, i_max+1)] = - Dx_Dxai_Dxai_min / Dx_Dxai_min;
            }
            if (fabs(Dx_Dxai_min) < fabs(Dy_Dxai_min)) {
                phi_vals_mat[offset2d(i, j_min, i_max+1)] = - Dy_Dxai_Dxai_min / Dy_Dxai_min;
            }
            if (fabs(Dx_Dxai_max) > fabs(Dy_Dxai_max)) {
                phi_vals_mat[offset2d(i, j_max, i_max+1)] = - Dx_Dxai_Dxai_max / Dx_Dxai_max;
            }
            if (fabs(Dx_Dxai_max) < fabs(Dy_Dxai_max)) {
                phi_vals_mat[offset2d(i, j_max, i_max+1)] = - Dy_Dxai_Dxai_max / Dy_Dxai_max;
            }
        } else {
            phi_vals_mat[offset2d(i_min, j, i_max+1)] = phi_valuse;
            phi_vals_mat[offset2d(i_max, j, i_max+1)] = phi_valuse;
        }
    }

    interpulat_mat(psi_vals_mat, 'i', i_max, j_max, i_min, j_min);
    interpulat_mat(phi_vals_mat, 'j', i_max, j_max, i_min, j_min);
}

/* copys matrix src to matrix src 
argument list:
dst - 1D array for the destination valus 
src - 1D array of the source valus */
void copy_mat(double *dst, double *src, int i_max, int j_max)
{
    int i, j;
    
    for (i = 0; i < i_max+1; i++) {
        for (j = 0; j < j_max+1; j++) {
            dst[offset2d(i, j, i_max+1)] = src[offset2d(i, j, i_max+1)];
        }
    }
}

/* copys the row 'src' into the j row in the destination matrix */
void copy_row_to_mat(double *dst, double *src, int row_num, int i_max)
{
    int i;
    
    for (i = 0; i < i_max+1; i++) {
        dst[offset2d(i, row_num, i_max+1)] = src[i];
    }
}

void copy_col_to_mat(double *dst, double *src, int col_num, int i_max, int j_max)
{
    int j;
    
    for (j = 0; j < j_max+1; j++) {
        dst[offset2d(col_num, j, i_max+1)] = src[j];
    }
}

/* returns the valus of L_x according to equation 10
argumetn list:
x_vals_mat - 1D array of the x valus
alpha_vals_mat - 1D array of the alpha valus 
phi_vals_mat - 1D array of the phi valus
beta_vals_mat - 1D array of the beta valus
gama_vals_mat - 1D array of the gama valus
psi_vals_mat - 1D array of the psi valus 
i, j = the point coordinate */
double L_x(double *x_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, int i_max, int j_max, int i, int j)
{
    double first_element, second_element, third_element, x_i_plus1_j, x_i_j,
    x_i_minus1_j, x_i_plus1_j_plus1, x_i_plus1_j_minus1, x_i_minus1_j_plus1,
    x_i_minus1_j_minus1, x_i_j_plus1, x_i_j_minus1;

    if (i == 0 || i == i_max || j == 0 || j == j_max) {
        return 0;
    }

    x_i_plus1_j = x_vals_mat[offset2d(i+1, j, i_max+1)];
    x_i_j = x_vals_mat[offset2d(i, j, i_max+1)];
    x_i_minus1_j = x_vals_mat[offset2d(i-1, j, i_max+1)];
    x_i_plus1_j_plus1 = x_vals_mat[offset2d(i+1,j+1, i_max+1)];
    x_i_plus1_j_minus1 = x_vals_mat[offset2d(i+1, j-1, i_max+1)];
    x_i_minus1_j_plus1 = x_vals_mat[offset2d(i-1, j+1, i_max+1)];
    x_i_minus1_j_minus1 = x_vals_mat[offset2d(i-1, j-1, i_max+1)];
    x_i_j_plus1 = x_vals_mat[offset2d(i, j+1, i_max+1)];
    x_i_j_minus1 = x_vals_mat[offset2d(i, j-1, i_max+1)];

    first_element = alpha_vals_mat[offset2d(i, j, i_max+1)] *
                    ((x_i_plus1_j - 2 * x_i_j + x_i_minus1_j) +
                    0.5 * phi_vals_mat[offset2d(i, j, i_max+1)] *
                    (x_i_plus1_j - x_i_minus1_j));
    second_element = 0.5 * beta_vals_mat[offset2d(i, j, i_max+1)] *
                     (x_i_plus1_j_plus1 - x_i_plus1_j_minus1 -
                     x_i_minus1_j_plus1 + x_i_minus1_j_minus1);
    third_element = gama_vals_mat[offset2d(i, j, i_max+1)] *
                    ((x_i_j_plus1 - 2 * x_i_j + x_i_j_minus1) +
                    0.5 * psi_vals_mat[offset2d(i, j, i_max+1)] *
                    (x_i_j_plus1 - x_i_j_minus1));

    return first_element - second_element + third_element;
}

/* returns the valus of L_y according to equation 11
argumetn list:
y_vals_mat - 1D array of the y valus
alpha_vals_mat - 1D array of the alpha valus 
phi_vals_mat - 1D array of the phi valus
beta_vals_mat - 1D array of the beta valus
gama_vals_mat - 1D array of the gama valus
phi_vals_mat - 1D array of the phi valus
psi_vals_mat - 1D array of the psi valus 
i, j = the point coordinate */
double L_y(double *y_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, int i_max, int j_max, int i, int j)
{
    double first_element, second_element, third_element, y_i_plus1_j, y_i_j,
    y_i_minus1_j, y_i_plus1_j_plus1, y_i_plus1_j_minus1, y_i_minus1_j_plus1,
    y_i_minus1_j_minus1, y_i_j_plus1, y_i_j_minus1;

    if (i == 0 || i == i_max || j == 0 || j == j_max) {
        return 0;
    }

    y_i_plus1_j = y_vals_mat[offset2d(i+1, j, i_max+1)];
    y_i_j = y_vals_mat[offset2d(i, j, i_max+1)];
    y_i_minus1_j = y_vals_mat[offset2d(i-1, j, i_max+1)];
    y_i_plus1_j_plus1 = y_vals_mat[offset2d(i+1,j+1, i_max+1)];
    y_i_plus1_j_minus1 = y_vals_mat[offset2d(i+1, j-1, i_max+1)];
    y_i_minus1_j_plus1 = y_vals_mat[offset2d(i-1, j+1, i_max+1)];
    y_i_minus1_j_minus1 = y_vals_mat[offset2d(i-1, j-1, i_max+1)];
    y_i_j_plus1 = y_vals_mat[offset2d(i, j+1, i_max+1)];
    y_i_j_minus1 = y_vals_mat[offset2d(i, j-1, i_max+1)];

    first_element = alpha_vals_mat[offset2d(i, j, i_max+1)] *
                    ((y_i_plus1_j - 2 * y_i_j + y_i_minus1_j) +
                    0.5 * phi_vals_mat[offset2d(i, j, i_max+1)] *
                    (y_i_plus1_j - y_i_minus1_j));
    second_element = 0.5 * beta_vals_mat[offset2d(i, j, i_max+1)] *
                     (y_i_plus1_j_plus1 - y_i_plus1_j_minus1 -
                     y_i_minus1_j_plus1 + y_i_minus1_j_minus1);
    third_element = gama_vals_mat[offset2d(i, j, i_max+1)] *
                    ((y_i_j_plus1 - 2 * y_i_j + y_i_j_minus1) +
                    0.5 * psi_vals_mat[offset2d(i, j, i_max+1)] *
                    (y_i_j_plus1 - y_i_j_minus1));

    return first_element - second_element + third_element;
}

/* doing the first sweep according to equation 12;
returns 0 on success
argument list:
fx_vals_mat - 1D array for the fx valus 
fy_vals_mat - 1D array for the fy valus
x_vals_mat - 1D array of the x valus
y_vals_mat - 1D array of the y valus
alpha_vals_mat - 1D array of the alpha valus 
phi_vals_mat - 1D array of the phi valus
beta_vals_mat - 1D array of the beta valus
gama_vals_mat - 1D array of the gama valus
psi_vals_mat - 1D array of the psi valus */
int sweep1(double *fx_vals_mat, double *fy_vals_mat, double *x_vals_mat, double *y_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, double *A, double *B, double *C,  double *D, double *temp_row, int i_max, int j_max, double r, double omega)
{
    int j_index, success = 0;

    /* solving for each j */
    for (j_index = 0; j_index < j_max+1; j_index++) {
        LHS_sweep1(A, B, C, alpha_vals_mat, j_index, i_max, r);
        RHS_sweep1_x(D, x_vals_mat, alpha_vals_mat, phi_vals_mat, beta_vals_mat, gama_vals_mat, psi_vals_mat, j_index, i_max, j_max, r, omega); 
        BC_sweep1(A, B, C, D, i_max);
        success = tridiag(A, B, C, D, temp_row, 0, i_max);
        if (success == 1) {
            break;
        }
        copy_row_to_mat(fx_vals_mat, temp_row, j_index, i_max);

        LHS_sweep1(A, B, C, alpha_vals_mat, j_index, i_max, r);
        RHS_sweep1_y(D, y_vals_mat, alpha_vals_mat, phi_vals_mat, beta_vals_mat, gama_vals_mat, psi_vals_mat, j_index, i_max, j_max, r, omega);
        BC_sweep1(A, B, C, D, i_max);
        success = tridiag(A, B, C, D, temp_row, 0, i_max);
        if (success == 1) {
            break;
        }
        copy_row_to_mat(fy_vals_mat, temp_row, j_index, i_max);
    }
    if (success == 1) {
        return 1;
    }

    return 0;
}

/* doing the first sweep according to equation 13;
returns 0 on success
argument list:
Cx_vals_mat - 1D array for the Cx valus 
Cy_vals_mat - 1D array for the Cy valus
fx_vals_mat - 1D array of the fx valus
fy_vals_mat - 1D array of the fy valus
gama_vals_mat - 1D array of the gama valus */
int sweep2(double *Cx_vals_mat, double *Cy_vals_mat, double *fx_vals_mat, double *fy_vals_mat, double *gama_vals_mat, double *A, double *B, double *C, double *D, double *temp_row, int i_max, int j_max, double r)
{
    int i_index, success = 0;

    /* solving for each i */
    for (i_index = 0; i_index < i_max+1; i_index++) {
        LHS_sweep2(A, B, C, gama_vals_mat, i_index, i_max, j_max, r);
        RHS_sweep2_x(D, fx_vals_mat, i_index, i_max, j_max);
        BC_sweep2(A, B, C, D, j_max);
        /*test*/
        // for (int index = 0; index < j_max+1; index++) {
        //     printf("%g\n", B[index]);
        // }
        /*test*/
        success = tridiag(A, B, C, D, temp_row, 0, j_max);
        if (success == 1) {
            printf("1\n");
            break;
        }
        copy_col_to_mat(Cx_vals_mat, temp_row, i_index, i_max, j_max);

        LHS_sweep2(A, B, C, gama_vals_mat, i_index, i_max, j_max, r);
        RHS_sweep2_y(D, fy_vals_mat, i_index, i_max, j_max);
        BC_sweep2(A, B, C, D, j_max);
        /*test*/
        // for (int index = 0; index < j_max+1; index++) {
        //     printf("%g\n", Bx2[index]);
        // }
        /*test*/
        success = tridiag(A, B, C, D, temp_row, 0, j_max);
        if (success == 1) {
            printf("2\n");
            break;
        }
    copy_col_to_mat(Cy_vals_mat, temp_row, i_index, i_max, j_max);
    }
    if (success == 1) {
        return 1;
    }

    return 0;
}

/* populates the A, B, C vectors accoding to eq 15 
argument list:
A - 1D array for the A valus 
B - 1D array for the B valus 
C - 1D array for the C valus
alpha_vals_mat - 1D array of the alpha valus
j - the row number */
void LHS_sweep1(double *A, double *B, double *C, double *alpha_vals_mat, int j, int i_max, double r)
{
    int i;

    for (i = 0; i < i_max+1; i++) {
        A[i] = -alpha_vals_mat[offset2d(i, j, i_max+1)];
        B[i] = r + 2 * alpha_vals_mat[offset2d(i, j, i_max+1)];
        C[i] = -alpha_vals_mat[offset2d(i, j, i_max+1)];
    }
}

/* populates the A, B, C vectors accoding to eq 15 
argument list:
A - 1D array for the A valus 
B - 1D array for the B valus 
C - 1D array for the C valus
alpha_vals_mat - 1D array of the alpha valus
j - the row number */
void LHS_sweep2(double *A, double *B, double *C, double *gama_vals_mat, int i, int i_max, int j_max, double r)
{
    int j;

    for (j = 0; j < j_max+1; j++) {
        A[j] = -gama_vals_mat[offset2d(i, j, i_max+1)];
        B[j] = r + 2 * gama_vals_mat[offset2d(i, j, i_max+1)];
        C[j] = -gama_vals_mat[offset2d(i, j, i_max+1)];
    }
}

/* populates the D vectors accoding to eq 15 
argument list:
D - 1D array for the D valus 
x_vals_mat - 1D array of the x valus
alpha_vals_mat - 1D array of the alpha valus
phi_vals_mat - 1D array of the phi valus
beta_vals_mat - 1D array of the beta valus
gama_vals_mat - 1D array of the gama valus
psi_vals_mat - 1D array of the psi valus 
j - the row number */
void RHS_sweep1_x(double *D, double *x_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, int j, int i_max, int j_max, double r, double omega)
{
    int i;
    
    for (i = 0; i < i_max+1; i++) {
        D[i] = r * omega * L_x(x_vals_mat, alpha_vals_mat, phi_vals_mat, beta_vals_mat, gama_vals_mat, psi_vals_mat, i_max, j_max, i, j);
    }
}

/* populates the D vectors accoding to eq 15 
argument list:
D - 1D array for the D valus 
fx_vals_mat - 1D array of the fx valus
i - the row number */
void RHS_sweep2_x(double *D,double *fx_vals_mat, int i, int i_max, int j_max)
{
    int j;
    
    for (j = 0; j < j_max+1; j++) {
        D[j] = fx_vals_mat[offset2d(i, j, i_max+1)];
    }
}

/* populates the D vectors accoding to eq 15 
argument list:
D - 1D array for the D valus 
y_vals_mat - 1D array of the y valus
alpha_vals_mat - 1D array of the alpha valus
phi_vals_mat - 1D array of the phi valus
beta_vals_mat - 1D array of the beta valus
gama_vals_mat - 1D array of the gama valus
psi_vals_mat - 1D array of the psi valus 
j - the row number */
void RHS_sweep1_y(double *D, double *y_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, int j, int i_max, int j_max, double r, double omega)
{
    int i;
    
    for (i = 0; i < i_max+1; i++) {
        D[i] = r * omega * L_y(y_vals_mat, alpha_vals_mat, phi_vals_mat, beta_vals_mat, gama_vals_mat, psi_vals_mat, i_max, j_max, i, j);
    }
}

/* populates the D vectors accoding to eq 15 
argument list:
D - 1D array for the D valus 
fy_vals_mat - 1D array of the fy valus
i - the row number */
void RHS_sweep2_y(double *D,double *fy_vals_mat, int i, int i_max, int j_max)
{
    int j;
    
    for (j = 0; j < j_max+1; j++) {
        D[j] = fy_vals_mat[offset2d(i, j, i_max+1)];
    }
}

/* applaying the boundary conditions on the A, B, C and D vectors according to eq 16
argument list:
A, B, C are the tridaig diaganols
D is the RHS vector */
void BC_sweep1(double *A, double *B, double *C, double *D, int i_max)
{
    A[0] = 0;
    B[0] = 1;
    C[0] = 0;
    D[0] = 0;
    A[i_max] = 0;
    B[i_max] = 1;
    C[i_max] = 0;
    D[i_max] = 0;
}

/* applaying the boundary conditions on the A, B, C and D vectors according to eq 16
argument list:
A, B, C are the tridaig diaganols
D is the RHS vector */
void BC_sweep2(double *A, double *B, double *C, double *D, int j_max)
{
    A[0] = 0;
    B[0] = 1;
    C[0] = 0;
    D[0] = 0;
    A[j_max] = 0;
    B[j_max] = 1;
    C[j_max] = 0;
    D[j_max] = 0;
}

/* a, b, c, are the vectors of the diagonal and the two
off-diagonals. The vector d is the RHS vector, the vector
u is the solution vector, "is" is the starting point, and
ie is the last point */
int tridiag(double *a, double *b, double *c, double *d, double *u, int is, int ie)
{
    int i;
    double beta;

    for (i = is + 1; i <= ie; i++)
    {
        if (b[i - 1] == 0.)
            return (1);
        beta = a[i] / b[i - 1];
        b[i] = b[i] - c[i - 1] * beta;
        d[i] = d[i] - d[i - 1] * beta;
    }

    u[ie] = d[ie] / b[ie];
    for (i = ie - 1; i >= is; i--)
    {
        u[i] = (d[i] - c[i] * u[i + 1]) / b[i];
    }
    return (0);
} 

/* doing the first sweep according to eq 9-17;
returns Vec2 with Lx_max, Ly_max (the residual)
argument list:
Cx_vals_mat - 1D array of the Cx valus 
Cy_vals_mat - 1D array of the Cy valus
fx_vals_mat - 1D array of the fx valus 
fy_vals_mat - 1D array of the fy valus
x_vals_mat_current - 1D array of the current x valus
x_vals_mat_next - 1D array of the next x valus
y_vals_mat_current - 1D array of the current y valus
y_vals_mat_next - 1D array of the next y valus
alpha_vals_mat - 1D array of the alpha valus 
phi_vals_mat - 1D array of the phi valus
beta_vals_mat - 1D array of the beta valus
gama_vals_mat - 1D array of the gama valus
psi_vals_mat - 1D array of the psi valus
A-D are temperary vectors for inverting the tri-diag matrices */
Vec2 step(double *Cx_vals_mat, double *Cy_vals_mat, double *fx_vals_mat, double *fy_vals_mat, double *x_vals_mat_current, double *x_vals_mat_next, double *y_vals_mat_current, double *y_vals_mat_next, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, double *A, double *B, double *C, double *D, double *temp_row, int i_max, int i_min, int j_max, int j_min, double r, double omega)
{
    int success, i, j, index; 
    Vec2 ans = {.x = 0, .y = 0};

    alpha_beta_gama(alpha_vals_mat, beta_vals_mat, gama_vals_mat, x_vals_mat_current, y_vals_mat_current, i_max, i_min, j_max, j_min);

    success = sweep1(fx_vals_mat, fy_vals_mat, x_vals_mat_current, y_vals_mat_current, alpha_vals_mat, phi_vals_mat, beta_vals_mat, gama_vals_mat, psi_vals_mat, A, B, C, D, temp_row, i_max, j_max, r, omega);
    if (success != 0) {
        ans.x = 1;
        return ans;
    }
    success = sweep2(Cx_vals_mat, Cy_vals_mat, fx_vals_mat, fy_vals_mat, gama_vals_mat, A, B, C, D, temp_row, i_max, j_max, r);
    if (success != 0) {
        ans.x = 2;
        return ans;
    }

    for (i = 0; i < i_max+1; i++) {
        for (j = 0; j < j_max+1; j++) {
            index = offset2d(i, j, i_max+1);
            x_vals_mat_next[index] = Cx_vals_mat[index] + x_vals_mat_current[index];
            y_vals_mat_next[index] = Cy_vals_mat[index] + y_vals_mat_current[index];
        }
    }

    ans.x = calculate_max_L_x(x_vals_mat_current, alpha_vals_mat, phi_vals_mat, beta_vals_mat, gama_vals_mat, psi_vals_mat, i_max, j_max);
    ans.y = calculate_max_L_y(y_vals_mat_current, alpha_vals_mat, phi_vals_mat, beta_vals_mat, gama_vals_mat, psi_vals_mat, i_max, j_max);

    return ans;
}

/* printing a 1D array to a file
argument list:
fp - file pointer
data - 1D array */
void mat_print_to_file(FILE *fp, double *data, int i_max, int j_max)
{
    int j_index, i_index;
    for (j_index = 0; j_index < j_max+1; j_index++) {
        for (i_index = 0; i_index < i_max+1; i_index++) {
            fprintf(fp, "%g ", data[offset2d(i_index, j_index, i_max+1)]);
        }
        fprintf(fp, "\n");
    }
}

/* printing a 1D array to the commend line
argument list:
data - 1D array */
void mat_print(double *data, int i_max, int j_max)
{
    int j_index, i_index;
    for (j_index = 0; j_index < j_max+1; j_index++) {
        for (i_index = 0; i_index < i_max+1; i_index++) {
            printf("%g ", data[offset2d(i_index, j_index, i_max+1)]);
        }
        printf("\n");
    }
}

/* returning the absolut maximum valus at the Lx matrix
argument list:
x_vals_mat - 1D array of the x valus
alpha_vals_mat - 1D array of the alpha valus 
phi_vals_mat - 1D array of the phi valus
beta_vals_mat - 1D array of the beta valus
gama_vals_mat - 1D array of the gama valus
psi_vals_mat - 1D array of the psi valus */
double calculate_max_L_x(double *x_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, int i_max, int j_max)
{
    int j_index, i_index;
    double max_L_x = 0, current_L_x;

    for (j_index = 0; j_index < j_max+1; j_index++) {
        for (i_index = 0; i_index < i_max+1; i_index++) {
            current_L_x = L_x(x_vals_mat, alpha_vals_mat, phi_vals_mat, beta_vals_mat, gama_vals_mat, psi_vals_mat, i_max, j_max, i_index, j_index);
            if (fabs(current_L_x) > max_L_x) {
                max_L_x = fabs(current_L_x);
            }
        }
    }
    return max_L_x;
}

/* returning the absolut maximum valus at the Ly matrix
argument list:
y_vals_mat - 1D array of the y valus
alpha_vals_mat - 1D array of the alpha valus 
phi_vals_mat - 1D array of the phi valus
beta_vals_mat - 1D array of the beta valus
gama_vals_mat - 1D array of the gama valus
psi_vals_mat - 1D array of the psi valus */
double calculate_max_L_y(double *y_vals_mat, double *alpha_vals_mat, double *phi_vals_mat, double *beta_vals_mat, double *gama_vals_mat, double *psi_vals_mat, int i_max, int j_max)
{
    int j_index, i_index;
    double max_L_y = 0, current_L_y;

    for (j_index = 0; j_index < j_max+1; j_index++) {
        for (i_index = 0; i_index < i_max+1; i_index++) {
            current_L_y = L_y(y_vals_mat, alpha_vals_mat, phi_vals_mat, beta_vals_mat, gama_vals_mat, psi_vals_mat, i_max, j_max, i_index, j_index);
            if (fabs(current_L_y) > max_L_y) {
                max_L_y = fabs(current_L_y);
            }
        }
    }
    return max_L_y;
}
