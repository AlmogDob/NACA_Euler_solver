#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sqlite3.h>
#include "mesher.h"

sqlite3 *setup_DB(char *db_name);

int main()
{
    sqlite3 *db = setup_DB("NACA.db");
    if (!db) {
        return 1;
    }
    char *temp_sql = "select ni, nj, x_mat, y_mat from NACA_data where NACA = 2250;";
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
    int ni = sqlite3_column_int(statement_pointer, 0);
    int nj = sqlite3_column_int(statement_pointer, 1);
    dprintINT(ni);
    dprintINT(nj);

    double *x_mat = (double *)malloc(sizeof(double) * (ni) * (nj));
    for (int i = 0; i < ni; i++) {   /* filling the matrix with zeros */
        for (int j = 0; j < nj; j++) {
            x_mat[offset2d(i, j, ni)] = 0;
        }
    }
    double *y_mat = (double *)malloc(sizeof(double) * (ni) * (nj));
    for (int i = 0; i < ni; i++) {   /* filling the matrix with zeros */
        for (int j = 0; j < nj; j++) {
            y_mat[offset2d(i, j, ni)] = 0;
        }
    }

    const void *blob_data = sqlite3_column_blob(statement_pointer, 2);
    int blob_size = sqlite3_column_bytes(statement_pointer, 2);

    if (blob_data == NULL) {
        fprintf(stderr, "Error: blob data is NULL\n");
        return 1;
    }

    memcpy(x_mat, sqlite3_column_blob(statement_pointer,2), ni * nj * sizeof(double));
    memcpy(y_mat, sqlite3_column_blob(statement_pointer,3), ni * nj * sizeof(double));

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
