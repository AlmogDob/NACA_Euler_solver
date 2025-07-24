#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sqlite3.h>
#include "mesher.h"
#include "solver.h"
#define MATRIX2D_IMPLEMENTATION
#include "Matrix2D.h"

sqlite3 *setup_DB(char *db_name);
int get_data_from_DB(double **x_2Dmat, double **y_2Dmat, double **rho_2Dmat, double **u_2Dmat, double **v_2Dmat, double **e_2Dmat, int *ni, int *nj, int *num_points_on_airfoil, double *Gamma, double *Mach_inf, double *rho_inf, double *p_inf, int NACA);

int main()
{
    double *x_2Dmat, *y_2Dmat, *rho_2Dmat, *u_2Dmat, *v_2Dmat, *e_2Dmat, Gamma, rho_inf, p_inf, Mach_inf;
    int ni, nj, num_points_on_airfoil;

    get_data_from_DB(&x_2Dmat, &y_2Dmat, &rho_2Dmat, &u_2Dmat, &v_2Dmat, &e_2Dmat, &ni, &nj, &num_points_on_airfoil, &Gamma, &Mach_inf, &rho_inf, &p_inf, 12);
    dprintINT(ni);
    dprintINT(nj);
    dprintINT(num_points_on_airfoil);

    double speed_of_sound_inf = sqrt(Gamma * p_inf / rho_inf);

    double velocity_inf = Mach_inf * speed_of_sound_inf;

    const int i_LE  = (ni-1) / 2;
    const int i_TEL = i_LE - num_points_on_airfoil / 2;
    const int i_TEU = i_LE + num_points_on_airfoil / 2;

    Mat2D airfoil_points = mat2D_alloc(num_points_on_airfoil, 2);
    Mat2D parallel_to_airfoil_vecs = mat2D_alloc(num_points_on_airfoil-1, 2);
    Mat2D normal_to_airfoil_vecs = mat2D_alloc(num_points_on_airfoil-1, 2);
    Mat2D force_vecs = mat2D_alloc(num_points_on_airfoil-1, 2);
    Vec2 tot_force = {0};

    for (int i = i_TEL; i <= i_TEU; i++) {
        MAT2D_AT(airfoil_points, i-i_TEL, 0) = x_2Dmat[offset2d_solver(i, 0, ni, nj)];
        MAT2D_AT(airfoil_points, i-i_TEL, 1) = y_2Dmat[offset2d_solver(i, 0, ni, nj)];
    }
    // MAT2D_PRINT(airfoil_points);

    /* calculating the parallel vector */
    for (int i = i_TEL; i < i_LE; i++) {
        MAT2D_AT(parallel_to_airfoil_vecs, i-i_TEL, 0) = MAT2D_AT(airfoil_points, i-i_TEL+1, 0) - MAT2D_AT(airfoil_points, i-i_TEL, 0);
        MAT2D_AT(parallel_to_airfoil_vecs, i-i_TEL, 1) = MAT2D_AT(airfoil_points, i-i_TEL+1, 1) - MAT2D_AT(airfoil_points, i-i_TEL, 1);
    }
    for (int i = i_TEU; i > i_LE; i--) {
        MAT2D_AT(parallel_to_airfoil_vecs, i-i_TEL-1, 0) = MAT2D_AT(airfoil_points, i-i_TEL, 0) - MAT2D_AT(airfoil_points, i-i_TEL-1, 0);
        MAT2D_AT(parallel_to_airfoil_vecs, i-i_TEL-1, 1) = MAT2D_AT(airfoil_points, i-i_TEL, 1) - MAT2D_AT(airfoil_points, i-i_TEL-1, 1);
    }
    for (int i = 0; i < num_points_on_airfoil-1; i++) {
        /* normalizing the vector */
        double size = sqrt(MAT2D_AT(parallel_to_airfoil_vecs, i, 0) * MAT2D_AT(parallel_to_airfoil_vecs, i, 0) + MAT2D_AT(parallel_to_airfoil_vecs, i, 1) * MAT2D_AT(parallel_to_airfoil_vecs, i, 1));
        MAT2D_AT(parallel_to_airfoil_vecs, i, 0) /= size;
        MAT2D_AT(parallel_to_airfoil_vecs, i, 1) /= size;

        /* calculating the normal vector */
        MAT2D_AT(normal_to_airfoil_vecs, i, 0) = - MAT2D_AT(parallel_to_airfoil_vecs, i, 1);
        MAT2D_AT(normal_to_airfoil_vecs, i, 1) = MAT2D_AT(parallel_to_airfoil_vecs, i, 0);

    }
    /* calculating the force */
    for (int i = i_TEL; i < i_LE; i++) {
        double size = sqrt(MAT2D_AT(parallel_to_airfoil_vecs, i-i_TEL, 0) * MAT2D_AT(parallel_to_airfoil_vecs, i-i_TEL, 0) + MAT2D_AT(parallel_to_airfoil_vecs, i-i_TEL, 1) * MAT2D_AT(parallel_to_airfoil_vecs, i-i_TEL, 1));

        double pressure = calculate_p(e_2Dmat[offset2d_solver(i, 0, ni, nj)], rho_2Dmat[offset2d_solver(i, 0, ni, nj)], u_2Dmat[offset2d_solver(i, 0, ni, nj)], v_2Dmat[offset2d_solver(i, 0, ni, nj)], Gamma);
        double force_size = pressure * size;
        // dprintD(force_size);
        MAT2D_AT(force_vecs, i-i_TEL, 0) = MAT2D_AT(normal_to_airfoil_vecs, i-i_TEL, 0) * force_size;
        MAT2D_AT(force_vecs, i-i_TEL, 1) = MAT2D_AT(normal_to_airfoil_vecs, i-i_TEL, 1) * force_size;
    }
    for (int i = i_TEU; i > i_LE; i--) {
        double size = sqrt(MAT2D_AT(parallel_to_airfoil_vecs, i-i_TEL-1, 0) * MAT2D_AT(parallel_to_airfoil_vecs, i-i_TEL-1, 0) + MAT2D_AT(parallel_to_airfoil_vecs, i-i_TEL-1, 1) * MAT2D_AT(parallel_to_airfoil_vecs, i-i_TEL-1, 1));

        double pressure = calculate_p(e_2Dmat[offset2d_solver(i, 0, ni, nj)], rho_2Dmat[offset2d_solver(i, 0, ni, nj)], u_2Dmat[offset2d_solver(i, 0, ni, nj)], v_2Dmat[offset2d_solver(i, 0, ni, nj)], Gamma);
        double force_size = pressure * size;
        // dprintD(force_size);
        MAT2D_AT(force_vecs, i-i_TEL-1, 0) = MAT2D_AT(normal_to_airfoil_vecs, i-i_TEL-1, 0) * force_size;
        MAT2D_AT(force_vecs, i-i_TEL-1, 1) = MAT2D_AT(normal_to_airfoil_vecs, i-i_TEL-1, 1) * force_size;
    }
    // MAT2D_PRINT(parallel_to_airfoil_vecs);
    // MAT2D_PRINT(normal_to_airfoil_vecs);
    // MAT2D_PRINT(force_vecs);

    for (int i = 0; i < num_points_on_airfoil-1; i++) {
        tot_force.x += MAT2D_AT(force_vecs, i, 0);
        tot_force.y += MAT2D_AT(force_vecs, i, 1);
    }
    printf("tot force: (%f, %f)\n", tot_force.x, tot_force.y);
    
    Vec2 tot_coeff = {.x = tot_force.x / (0.5 * rho_inf * velocity_inf * velocity_inf), .y = tot_force.y / (0.5 * rho_inf * velocity_inf * velocity_inf)}; /* L = 0.5*rho*v^2*c*CL */
    printf("tot coeff: (%f, %f)\n", tot_coeff.x, tot_coeff.y);

    return 0;
}

sqlite3 *setup_DB(char *db_name)
{
    sqlite3 *db;
    int rc = sqlite3_open(db_name, &db);

    if (rc != SQLITE_OK) {
        fprintf(stderr, "%s:%d: [ERROR] cannot open database %s\n", __FILE__, __LINE__, sqlite3_errmsg(db));
        return NULL;
    }

    return db;
}

/* returns zero on success */
int get_data_from_DB(double **x_2Dmat, double **y_2Dmat, double **rho_2Dmat, double **u_2Dmat, double **v_2Dmat, double **e_2Dmat, int *ni, int *nj, int *num_points_on_airfoil, double *Gamma, double *Mach_inf, double *rho_inf, double *p_inf, int NACA)
{
    sqlite3 *db = setup_DB("NACA.db");
    if (!db) {
        return 1;
    }
    char temp_sql[MAXWORD];
    sprintf(temp_sql, "select ni, nj, num_points_on_airfoil, Gamma, Mach_inf, density, environment_pressure, x_2Dmat, y_2Dmat, rho_2Dmat, u_2Dmat, v_2Dmat, e_2Dmat from NACA_data where NACA = %d;", NACA);
    sqlite3_stmt *statement_pointer;
    int rc = sqlite3_prepare_v2(db, temp_sql, -1, &statement_pointer, 0);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "%s:%d: [ERROR] failed to prepare statement %s\n", __FILE__, __LINE__, temp_sql);
        return SQLITE_ERROR;
    }
    rc = sqlite3_step(statement_pointer);
    if (rc != SQLITE_ROW) {
        fprintf(stderr, "%s:%d: [ERROR] failed to step the statement\n", __FILE__, __LINE__);
        return SQLITE_ERROR;
    }
    *ni = sqlite3_column_int(statement_pointer, 0);
    *nj = sqlite3_column_int(statement_pointer, 1);
    *num_points_on_airfoil = sqlite3_column_int(statement_pointer, 2); 
    *Gamma = sqlite3_column_double(statement_pointer, 3);
    *Mach_inf = sqlite3_column_double(statement_pointer, 4);
    *rho_inf = sqlite3_column_double(statement_pointer, 5);
    *p_inf = sqlite3_column_double(statement_pointer, 6);

    *x_2Dmat = (double *)malloc(sizeof(double) * (*ni) * (*nj));
    for (int i = 0; i < *ni; i++) {   /* filling the matrix with zeros */
        for (int j = 0; j < *nj; j++) {
            (*x_2Dmat)[offset2d_mesher(i, j, *ni)] = 0;
        }
    }
    *y_2Dmat = (double *)malloc(sizeof(double) * (*ni) * (*nj));
    for (int i = 0; i < *ni; i++) {   /* filling the matrix with zeros */
        for (int j = 0; j < *nj; j++) {
            (*y_2Dmat)[offset2d_mesher(i, j, *ni)] = 0;
        }
    }
    *rho_2Dmat = (double *)malloc(sizeof(double) * (*ni) * (*nj));
    for (int i = 0; i < *ni; i++) {   /* filling the matrix with zeros */
        for (int j = 0; j < *nj; j++) {
            (*rho_2Dmat)[offset2d_mesher(i, j, *ni)] = 0;
        }
    }
    *u_2Dmat = (double *)malloc(sizeof(double) * (*ni) * (*nj));
    for (int i = 0; i < *ni; i++) {   /* filling the matrix with zeros */
        for (int j = 0; j < *nj; j++) {
            (*u_2Dmat)[offset2d_mesher(i, j, *ni)] = 0;
        }
    }
    *v_2Dmat = (double *)malloc(sizeof(double) * (*ni) * (*nj));
    for (int i = 0; i < *ni; i++) {   /* filling the matrix with zeros */
        for (int j = 0; j < *nj; j++) {
            (*v_2Dmat)[offset2d_mesher(i, j, *ni)] = 0;
        }
    }
    *e_2Dmat = (double *)malloc(sizeof(double) * (*ni) * (*nj));
    for (int i = 0; i < *ni; i++) {   /* filling the matrix with zeros */
        for (int j = 0; j < *nj; j++) {
            (*e_2Dmat)[offset2d_mesher(i, j, *ni)] = 0;
        }
    }

    const void *blob_data_x = sqlite3_column_blob(statement_pointer, 7);
    if (blob_data_x == NULL) {
        fprintf(stderr, "Error: x_2Dmat is NULL\n");
        return 1;
    }
    memcpy(*x_2Dmat, sqlite3_column_blob(statement_pointer,7), *ni * *nj * sizeof(double));

    const void *blob_data_y = sqlite3_column_blob(statement_pointer, 8);
    if (blob_data_y == NULL) {
        fprintf(stderr, "Error: y_2Dmat is NULL\n");
        return 1;
    }
    memcpy(*y_2Dmat, sqlite3_column_blob(statement_pointer,8), *ni * *nj * sizeof(double));

    const void *blob_data_rho = sqlite3_column_blob(statement_pointer, 9);
    if (blob_data_rho == NULL) {
        fprintf(stderr, "Error: rho_2Dmat is NULL\n");
        return 1;
    }
    memcpy(*rho_2Dmat, sqlite3_column_blob(statement_pointer,9), *ni * *nj * sizeof(double));

    const void *blob_data_u = sqlite3_column_blob(statement_pointer, 10);
    if (blob_data_u == NULL) {
        fprintf(stderr, "Error: u_2Dmat is NULL\n");
        return 1;
    }
    memcpy(*u_2Dmat, sqlite3_column_blob(statement_pointer,10), *ni * *nj * sizeof(double));

    const void *blob_data_v = sqlite3_column_blob(statement_pointer, 11);
    if (blob_data_v == NULL) {
        fprintf(stderr, "Error: v_2Dmat is NULL\n");
        return 1;
    }
    memcpy(*v_2Dmat, sqlite3_column_blob(statement_pointer,11), *ni * *nj * sizeof(double));

    const void *blob_data_e = sqlite3_column_blob(statement_pointer, 12);
    if (blob_data_e == NULL) {
        fprintf(stderr, "Error: e_2Dmat is NULL\n");
        return 1;
    }
    memcpy(*e_2Dmat, sqlite3_column_blob(statement_pointer,12), *ni * *nj * sizeof(double));

    // mat_print(*x_2Dmat, *ni-1, *nj-1);

    return 0;
}
