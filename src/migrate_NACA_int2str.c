#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sqlite3.h>
#include "mesher.h"
#include "Almog_Dynamic_Array.h"

typedef struct {
    char NACA[256];
    int ni;
    int nj;
    int num_points_on_airfoil;
    double delta_y;
    double XSF;
    double YSF;
    double r;
    double omega;
    double Mach_inf;
    double angle_of_attack_deg;
    double density;
    double environment_pressure;
    double delta_t;
    double Gamma;
    double epse;
    double CL;
    double CD;
    double *x_2Dmat;
    double *y_2Dmat;
    double *rho_2Dmat;
    double *u_2Dmat;
    double *v_2Dmat;
    double *e_2Dmat;
} Params;

typedef struct {
    size_t length;
    size_t capacity;
    int* elements;
} ada_int_array;

sqlite3 *setup_DB_int(char *db_name);
sqlite3 *setup_DB_str(char *db_name);
int get_data_from_DB_int(sqlite3 *db, int ID, Params *params);
int save_to_DB(sqlite3 *db, Params params);
int get_ID_list_from_db(sqlite3 *db, ada_int_array *IDs);

int main()
{
    Params params;
    ada_int_array IDs;
    ada_init_array(int, IDs);

    /* setup DB */
    printf("[INFO] setting up DB str\n");
    sqlite3 *db_str = setup_DB_str("NACA.db");
    if (!db_str) {
        return 1;
    }
    printf("[INFO] setting up DB int\n");
    sqlite3 *db_int = setup_DB_int("NACA_int.db");
    if (!db_int) {
        return 1;
    }

    printf("[INFO] getting data from int DB and saving to str DB\n");
    if (get_ID_list_from_db(db_int, &IDs) != SQLITE_OK) {
        return 1;
    }

    for (size_t i = 0; i < IDs.length; i++) {
        printf("ID: %d\n", IDs.elements[i]);

        if (get_data_from_DB_int(db_int, IDs.elements[i], &params) != SQLITE_OK) {
            return 1;
        }
        if (save_to_DB(db_str, params) != SQLITE_OK) {
            return 1;
        }
    }



    sqlite3_close(db_str);
    free(params.x_2Dmat);
    free(params.y_2Dmat);
    free(params.rho_2Dmat);
    free(params.u_2Dmat);
    free(params.v_2Dmat);
    free(params.e_2Dmat);

    return 0;
}

sqlite3 *setup_DB_int(char *db_name)
{
    sqlite3 *db;
    char *err_msg = 0;
    int rc = sqlite3_open(db_name, &db);

    if (rc != SQLITE_OK) {
        fprintf(stderr, "%s:%d: [ERROR] cannot open database %s\n", __FILE__, __LINE__, sqlite3_errmsg(db));
        return NULL;
    }

    sqlite3_busy_timeout(db, 5e3);

    char *sql = "CREATE TABLE IF NOT EXISTS NACA_data(ID INTEGER PRIMARY KEY, NACA INTEGER NOT NULL, ni INTEGER NOT NULL, nj INTEGER NOT NULL, num_points_on_airfoil INTEGER NOT NULL, delta_y REAL NOT NULL, XSF REAL NOT NULL, YSF REAL NOT NULL, r REAL NOT NULL, omega REAL NOT NULL, Mach_inf REAL NOT NULL, angle_of_attack_deg REAL NOT NULL, density REAL NOT NULL, environment_pressure REAL NOT NULL, delta_t REAL NOT NULL, Gamma REAL NOT NULL, epse REAL NOT NULL, CL REAL NOT NULL, CD REAL NOT NULL, x_2Dmat BLOB, y_2Dmat BLOB, rho_2Dmat BLOB, u_2Dmat BLOB, v_2Dmat BLOB, e_2Dmat BLOB)";
    rc = sqlite3_exec(db, sql, 0, 0, &err_msg);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "%s:%d: [ERROR] cannot exec command: %s\n", __FILE__, __LINE__, sqlite3_errmsg(db));
        sqlite3_free(err_msg);
        sqlite3_close(db);
        return NULL;
    }

    return db;
}

sqlite3 *setup_DB_str(char *db_name)
{
    sqlite3 *db;
    char *err_msg = 0;
    int rc = sqlite3_open(db_name, &db);

    if (rc != SQLITE_OK) {
        fprintf(stderr, "%s:%d: [ERROR] cannot open database %s\n", __FILE__, __LINE__, sqlite3_errmsg(db));
        return NULL;
    }

    sqlite3_busy_timeout(db, 5e3);

    char *sql = "CREATE TABLE IF NOT EXISTS NACA_data(ID INTEGER PRIMARY KEY, NACA TEXT NOT NULL, ni INTEGER NOT NULL, nj INTEGER NOT NULL, num_points_on_airfoil INTEGER NOT NULL, delta_y REAL NOT NULL, XSF REAL NOT NULL, YSF REAL NOT NULL, r REAL NOT NULL, omega REAL NOT NULL, Mach_inf REAL NOT NULL, angle_of_attack_deg REAL NOT NULL, density REAL NOT NULL, environment_pressure REAL NOT NULL, delta_t REAL NOT NULL, Gamma REAL NOT NULL, epse REAL NOT NULL, CL REAL NOT NULL, CD REAL NOT NULL, x_2Dmat BLOB, y_2Dmat BLOB, rho_2Dmat BLOB, u_2Dmat BLOB, v_2Dmat BLOB, e_2Dmat BLOB)";
    rc = sqlite3_exec(db, sql, 0, 0, &err_msg);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "%s:%d: [ERROR] cannot exec command: %s\n", __FILE__, __LINE__, sqlite3_errmsg(db));
        sqlite3_free(err_msg);
        sqlite3_close(db);
        return NULL;
    }

    return db;
}

/* returns zero on success */
int get_data_from_DB_int(sqlite3 *db, int ID, Params *params)
{
    char temp_sql[2048];
    sprintf(temp_sql, "select NACA, ni, nj, num_points_on_airfoil, delta_y, XSF, YSF, r, omega, Mach_inf, angle_of_attack_deg, density, environment_pressure, delta_t, Gamma, epse, CL, CD, x_2Dmat, y_2Dmat, rho_2Dmat, u_2Dmat, v_2Dmat, e_2Dmat from NACA_data where ID = %d;", ID);
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

    int param_index = 0;
    int NACA_int = sqlite3_column_int(statement_pointer, param_index++);
    if (NACA_int < 100) {
        sprintf(params->NACA, "00%d", NACA_int);
    } else {
        sprintf(params->NACA, "%d", NACA_int);
    }
    params->ni = sqlite3_column_int(statement_pointer, param_index++);
    params->nj = sqlite3_column_int(statement_pointer, param_index++);
    params->num_points_on_airfoil = sqlite3_column_int(statement_pointer, param_index++); 
    params->delta_y = sqlite3_column_double(statement_pointer, param_index++);
    params->XSF = sqlite3_column_double(statement_pointer, param_index++);
    params->YSF = sqlite3_column_double(statement_pointer, param_index++);
    params->r = sqlite3_column_double(statement_pointer, param_index++);
    params->omega = sqlite3_column_double(statement_pointer, param_index++);
    params->Mach_inf = sqlite3_column_double(statement_pointer, param_index++);
    params->angle_of_attack_deg = sqlite3_column_double(statement_pointer, param_index++);
    params->density = sqlite3_column_double(statement_pointer, param_index++);
    params->environment_pressure = sqlite3_column_double(statement_pointer, param_index++);
    params->delta_t = sqlite3_column_double(statement_pointer, param_index++);
    params->Gamma = sqlite3_column_double(statement_pointer, param_index++);
    params->epse = sqlite3_column_double(statement_pointer, param_index++);
    params->CL = sqlite3_column_double(statement_pointer, param_index++);
    params->CD = sqlite3_column_double(statement_pointer, param_index++);

    params->x_2Dmat = (double *)malloc(sizeof(double) * (params->ni) * (params->nj));
    for (int i = 0; i < params->ni; i++) {   /* filling the matrix with zeros */
        for (int j = 0; j < params->nj; j++) {
            (params->x_2Dmat)[offset2d_mesher(i, j, params->ni)] = 0;
        }
    }
    params->y_2Dmat = (double *)malloc(sizeof(double) * (params->ni) * (params->nj));
    for (int i = 0; i < params->ni; i++) {   /* filling the matrix with zeros */
        for (int j = 0; j < params->nj; j++) {
            (params->y_2Dmat)[offset2d_mesher(i, j, params->ni)] = 0;
        }
    }
    params->rho_2Dmat = (double *)malloc(sizeof(double) * (params->ni) * (params->nj));
    for (int i = 0; i < params->ni; i++) {   /* filling the matrix with zeros */
        for (int j = 0; j < params->nj; j++) {
            (params->rho_2Dmat)[offset2d_mesher(i, j, params->ni)] = 0;
        }
    }
    params->u_2Dmat = (double *)malloc(sizeof(double) * (params->ni) * (params->nj));
    for (int i = 0; i < params->ni; i++) {   /* filling the matrix with zeros */
        for (int j = 0; j < params->nj; j++) {
            (params->u_2Dmat)[offset2d_mesher(i, j, params->ni)] = 0;
        }
    }
    params->v_2Dmat = (double *)malloc(sizeof(double) * (params->ni) * (params->nj));
    for (int i = 0; i < params->ni; i++) {   /* filling the matrix with zeros */
        for (int j = 0; j < params->nj; j++) {
            (params->v_2Dmat)[offset2d_mesher(i, j, params->ni)] = 0;
        }
    }
    params->e_2Dmat = (double *)malloc(sizeof(double) * (params->ni) * (params->nj));
    for (int i = 0; i < params->ni; i++) {   /* filling the matrix with zeros */
        for (int j = 0; j < params->nj; j++) {
            (params->e_2Dmat)[offset2d_mesher(i, j, params->ni)] = 0;
        }
    }

    const void *blob_data_x = sqlite3_column_blob(statement_pointer, param_index);
    if (blob_data_x == NULL) {
        fprintf(stderr, "Error: x_2Dmat is NULL\n");
        return 1;
    }
    memcpy(params->x_2Dmat, sqlite3_column_blob(statement_pointer,param_index++), params->ni * params->nj * sizeof(double));

    const void *blob_data_y = sqlite3_column_blob(statement_pointer, param_index);
    if (blob_data_y == NULL) {
        fprintf(stderr, "Error: y_2Dmat is NULL\n");
        return 1;
    }
    memcpy(params->y_2Dmat, sqlite3_column_blob(statement_pointer,param_index++), params->ni * params->nj * sizeof(double));

    const void *blob_data_rho = sqlite3_column_blob(statement_pointer, param_index);
    if (blob_data_rho == NULL) {
        fprintf(stderr, "Error: rho_2Dmat is NULL\n");
        return 1;
    }
    memcpy(params->rho_2Dmat, sqlite3_column_blob(statement_pointer,param_index++), params->ni * params->nj * sizeof(double));

    const void *blob_data_u = sqlite3_column_blob(statement_pointer, param_index);
    if (blob_data_u == NULL) {
        fprintf(stderr, "Error: u_2Dmat is NULL\n");
        return 1;
    }
    memcpy(params->u_2Dmat, sqlite3_column_blob(statement_pointer,param_index++), params->ni * params->nj * sizeof(double));

    const void *blob_data_v = sqlite3_column_blob(statement_pointer, param_index);
    if (blob_data_v == NULL) {
        fprintf(stderr, "Error: v_2Dmat is NULL\n");
        return 1;
    }
    memcpy(params->v_2Dmat, sqlite3_column_blob(statement_pointer,param_index++), params->ni * params->nj * sizeof(double));

    const void *blob_data_e = sqlite3_column_blob(statement_pointer, param_index);
    if (blob_data_e == NULL) {
        fprintf(stderr, "Error: e_2Dmat is NULL\n");
        return 1;
    }
    memcpy(params->e_2Dmat, sqlite3_column_blob(statement_pointer,param_index++), params->ni * params->nj * sizeof(double));


    return 0;
}

/* saving input param mesh and solution to the DB and deleting if there are duplicates;
returning the error return code. */
int save_to_DB(sqlite3 *db, Params params)
{
    /* saving to DB */
    char temp_sql[2048];
    char *err_msg = 0;

    strcpy(temp_sql, "");
    sprintf(temp_sql, "insert into NACA_data(NACA, ni, nj, num_points_on_airfoil, delta_y, XSF, YSF, r, omega, Mach_inf, angle_of_attack_deg, density, environment_pressure, delta_t, Gamma, epse, CL, CD, x_2Dmat, y_2Dmat, rho_2Dmat, u_2Dmat, v_2Dmat, e_2Dmat) values(%s, %d, %d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, ?, ?, ?, ?, ?, ?);", params.NACA, params.ni, params.nj, params.num_points_on_airfoil, params.delta_y, params.XSF, params.YSF, params.r, params.omega, params.Mach_inf, params.angle_of_attack_deg, params.density, params.environment_pressure, params.delta_t, params.Gamma, params.epse, params.CL, params.CD);
    sqlite3_stmt *statement_pointer;
    int rc = sqlite3_prepare_v2(db, temp_sql, -1, &statement_pointer, 0);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "%s:%d: [ERROR] failed to prepare statement %s\n", __FILE__, __LINE__, temp_sql);
        return SQLITE_ERROR;
    }
    
    /* binding mesh */
    rc = sqlite3_bind_blob(statement_pointer, 1, params.x_2Dmat, params.ni * params.nj * sizeof(double), SQLITE_STATIC);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "%s:%d: [ERROR] failed to bind x values\n", __FILE__, __LINE__);
        return SQLITE_ERROR;
    }
    rc = sqlite3_bind_blob(statement_pointer, 2, params.y_2Dmat, params.ni * params.nj * sizeof(double), SQLITE_STATIC);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "%s:%d: [ERROR] failed to bind y values\n", __FILE__, __LINE__);
        return SQLITE_ERROR;
    }
    rc = sqlite3_bind_blob(statement_pointer, 3, params.rho_2Dmat, params.ni * params.nj * sizeof(double), SQLITE_STATIC);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "%s:%d: [ERROR] failed to bind rho values\n", __FILE__, __LINE__);
        return SQLITE_ERROR;
    }
    rc = sqlite3_bind_blob(statement_pointer, 4, params.u_2Dmat, params.ni * params.nj * sizeof(double), SQLITE_STATIC);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "%s:%d: [ERROR] failed to bind u values\n", __FILE__, __LINE__);
        return SQLITE_ERROR;
    }
    rc = sqlite3_bind_blob(statement_pointer, 5, params.v_2Dmat, params.ni * params.nj * sizeof(double), SQLITE_STATIC);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "%s:%d: [ERROR] failed to bind v values\n", __FILE__, __LINE__);
        return SQLITE_ERROR;
    }
    rc = sqlite3_bind_blob(statement_pointer, 6, params.e_2Dmat, params.ni * params.nj * sizeof(double), SQLITE_STATIC);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "%s:%d: [ERROR] failed to bind e values\n", __FILE__, __LINE__);
        return SQLITE_ERROR;
    }

    /* executing statement */
    rc = sqlite3_step(statement_pointer);
    if (rc != SQLITE_DONE) {
        fprintf(stderr, "%s:%d: [ERROR] failed to step the statement of x and y and Q mat\n", __FILE__, __LINE__);
        return SQLITE_ERROR;
    }
    rc = sqlite3_finalize(statement_pointer);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "%s:%d: [ERROR] failed to finalize the statement of x and y and Q mat\n", __FILE__, __LINE__);
        return SQLITE_ERROR;
    }

    /* deleting duplicates */
    strcpy(temp_sql, "");
    sprintf(temp_sql, "DELETE FROM NACA_data WHERE ID NOT IN (SELECT MIN(ID) FROM NACA_data GROUP BY NACA, ni, nj, num_points_on_airfoil, delta_y, XSF, YSF, r, omega, Mach_inf, angle_of_attack_deg, density, environment_pressure, Gamma, epse);");
    err_msg = 0;
    rc = sqlite3_exec(db ,temp_sql, 0, 0, &err_msg);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "%s:%d: [ERROR] cannot write to DB %s\n", __FILE__, __LINE__, err_msg);
        return SQLITE_ERROR;
    }

    return rc;
}

int get_ID_list_from_db(sqlite3 *db, ada_int_array *IDs)
{
    sqlite3_stmt *stmt;
    const char *sql = "select ID from NACA_data;";
    int rc = sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "%s:%d: [ERROR] cannot read DB: %s\n", __FILE__, __LINE__, sqlite3_errmsg(db));
        return SQLITE_ERROR;
    }

    ada_int_array temp_IDs = *IDs;
    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) 
    {
        int ID = sqlite3_column_int(stmt, 0);
        ada_appand(int, temp_IDs, ID);
    }

    *IDs = temp_IDs;
    
    return SQLITE_OK;
}
