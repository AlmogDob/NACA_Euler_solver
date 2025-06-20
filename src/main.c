
#include "mesher.h"
#include <sqlite3.h>

typedef struct {
    int NACA;
    int ni;
    int nj;
    int num_points_on_airfoil;
    double delta_y;
    double XSF;
    double YSF;
    double r;
    double omega;
} Input_param;

#define ON_LINUX 0
#ifdef __linux__
    #include <sys/stat.h>
    #ifndef __USE_MISC
        #define __USE_MISC
    #endif
    #include <dirent.h>
    #undef ON_LINUX
    #define ON_LINUX 1
    int create_empty_dir(char *parent_directory);
    int create_output_dir(char *output_dir, Input_param input_param);
#endif

void read_input(char *input_file, Input_param *input_param);
void output_metadata(char *output_dir, Input_param input_param);
void mat_output_to_file(FILE *fp, double *data, Input_param input_param);
void output_mesh(char *output_dir, double *x_vals_mat, double *y_vals_mat, Input_param input_param);
sqlite3 *setup_DB(char * db_name);
int output_input_param_to_DB(sqlite3 *db, Input_param input_param);

int main(int argc, char const *argv[])
{
    /* Geting the input file and output directory */
    char input_file[MAXDIR], output_dir[MAXDIR]; 

    if (--argc != 2) {
        fprintf(stderr, "%s:%d: [ERROR] not right usage\nUsage: main 'input file'\n", __FILE__, __LINE__);
        return 1;
    }

    strncpy(input_file, (*(++argv)), MAXDIR);
    if (input_file[MAXDIR-1] != '\0') {
        fprintf(stderr, "%s:%d: [ERROR] input too long\n", __FILE__, __LINE__);
        return 1;
    }
    strncpy(output_dir, (*(++argv)), MAXDIR);
    if (output_dir[MAXDIR-1] != '\0') {
        fprintf(stderr, "%s:%d: [ERROR] input too long\n", __FILE__, __LINE__);
        return 1;
    }

    /* Input variables */
    Input_param input_param;

    read_input(input_file, &input_param);

    /* Checking that I got the right input */
    dprintINT(input_param.NACA);
    dprintINT(input_param.ni);
    dprintINT(input_param.nj);
    dprintINT(input_param.num_points_on_airfoil);
    dprintD(input_param.delta_y);
    dprintD(input_param.XSF);
    dprintD(input_param.YSF);
    dprintD(input_param.r);
    dprintD(input_param.omega);
    printf("--------------------\n");

    /* creating output directory */
    printf("[INFO] creating output directory\n");
    if (ON_LINUX) {
        if (create_output_dir(output_dir, input_param)) {
            return 1;
        }
    }

    /* creating mesh */
    printf("[INFO] meshing\n");
    double *x_mat, *y_mat;

    int mesh_rc = create_mesh(&x_mat, &y_mat, input_param.NACA, input_param.ni, input_param.nj, input_param.num_points_on_airfoil, input_param.delta_y, input_param.XSF, input_param.YSF, input_param.r, input_param.omega, output_dir);
    if (mesh_rc != 0) {
        fprintf(stderr, "%s:%d: [ERROR] creating mesh\n", __FILE__, __LINE__);
        return 1;
    }

    /* saving mesh */
    printf("[INFO] saving mesh\n");
    output_mesh(output_dir, x_mat, y_mat, input_param);
    output_metadata(output_dir, input_param);
    char *err_msg = 0;
    sqlite3 *db = setup_DB("NACA.db");
    if (!db) {
        return 1;
    }
    char *sql = "CREATE TABLE IF NOT EXISTS NACA_data(ID INTEGER PRIMARY KEY, NACA INTEGER, ni INTEGER, nj INTEGER, num_points_on_airfoil INTEGER, delta_y REAL, XSF REAL, YSF REAL, r REAL, omega REAL)";
    int rc = sqlite3_exec(db, sql, 0, 0, &err_msg);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "%s:%d: [ERROR] cannot exec command: %s\n", __FILE__, __LINE__, sqlite3_errmsg(db));
        sqlite3_free(err_msg);
        sqlite3_close(db);
        return 1;
    }

    if (output_input_param_to_DB(db, input_param) != SQLITE_OK) {
        return 1;
    }
    sqlite3_close(db);

    return 0;
}

#if ON_LINUX
/* create empty dir at 'parent directory'.
// if allready exists, delete all the files inside
returns 0 on success
this function handles the errors so on fail just quit
argument list:
parent_directory - char pointer to the directory name */
int create_empty_dir(char *parent_directory)
{
    char path_to_remove[BUFSIZ];

    if (mkdir(parent_directory, 0777) == -1) {
        if (errno == EEXIST) {
            DIR *dir = opendir(parent_directory);
            if (dir == NULL) {
                fprintf(stderr, "%s:%d: [ERROR] problem opening '%s': %s\n", __FILE__, __LINE__, parent_directory, strerror(errno));
                return 1;
            }
            struct dirent* entity;
            entity = readdir(dir);
            while (entity != NULL) {   /* directory is not empty */
                strncpy(path_to_remove, parent_directory, BUFSIZ);
                strncat(path_to_remove, "/", BUFSIZ/2);
                strncat(path_to_remove, entity->d_name, BUFSIZ);
                if (entity->d_type == DT_REG) {
                    if (remove(path_to_remove) != 0) {
                        fprintf(stderr, "%s:%d: [ERROR] problem removing '%s': %s\n", __FILE__, __LINE__, path_to_remove, strerror(errno));
                        return 1;
                    }
                }
                entity = readdir(dir);
            }
            closedir(dir);
            return 0;
        }

        fprintf(stderr, "%s:%d: [ERROR] problem making '%s': %s\n", __FILE__, __LINE__, parent_directory, strerror(errno));
        return 1;
    }
    return 0;
}

/* return non zero value on error */
int create_output_dir(char *output_dir, Input_param input_param)
{
    char temp_word[MAXDIR];

    if (create_empty_dir(output_dir) != 0) {
        fprintf(stderr, "%s:%d: [ERROR] creating ouput directory\n", __FILE__, __LINE__);
        return 1;
    }
    sprintf(temp_word, "/NACA%d", input_param.NACA);
    strcat(output_dir, temp_word);
    if (create_empty_dir(output_dir) != 0) {
        fprintf(stderr, "%s:%d: [ERROR] creating ouput directory\n", __FILE__, __LINE__);
        return 1;
    }
    sprintf(temp_word, "/ni%d", input_param.ni);
    strcat(output_dir, temp_word);
    if (create_empty_dir(output_dir) != 0) {
        fprintf(stderr, "%s:%d: [ERROR] creating ouput directory\n", __FILE__, __LINE__);
        return 1;
    }
    sprintf(temp_word, "/nj%d", input_param.nj);
    strcat(output_dir, temp_word);
    if (create_empty_dir(output_dir) != 0) {
        fprintf(stderr, "%s:%d: [ERROR] creating ouput directory\n", __FILE__, __LINE__);
        return 1;
    }
    sprintf(temp_word, "/num_points_on_airfoil%d", input_param.num_points_on_airfoil);
    strcat(output_dir, temp_word);
    if (create_empty_dir(output_dir) != 0) {
        fprintf(stderr, "%s:%d: [ERROR] creating ouput directory\n", __FILE__, __LINE__);
        return 1;
    }
    sprintf(temp_word, "/delta_y%g", input_param.delta_y);
    strcat(output_dir, temp_word);
    if (create_empty_dir(output_dir) != 0) {
        fprintf(stderr, "%s:%d: [ERROR] creating ouput directory\n", __FILE__, __LINE__);
        return 1;
    }
    sprintf(temp_word, "/XSF%g", input_param.XSF);
    strcat(output_dir, temp_word);
    if (create_empty_dir(output_dir) != 0) {
        fprintf(stderr, "%s:%d: [ERROR] creating ouput directory\n", __FILE__, __LINE__);
        return 1;
    }
    sprintf(temp_word, "/YSF%g", input_param.YSF);
    strcat(output_dir, temp_word);
    if (create_empty_dir(output_dir) != 0) {
        fprintf(stderr, "%s:%d: [ERROR] creating ouput directory\n", __FILE__, __LINE__);
        return 1;
    }
    sprintf(temp_word, "/r%g", input_param.r);
    strcat(output_dir, temp_word);
    if (create_empty_dir(output_dir) != 0) {
        fprintf(stderr, "%s:%d: [ERROR] creating ouput directory\n", __FILE__, __LINE__);
        return 1;
    }
    sprintf(temp_word, "/omega%g", input_param.omega);
    strcat(output_dir, temp_word);
    if (create_empty_dir(output_dir) != 0) {
        fprintf(stderr, "%s:%d: [ERROR] creating ouput directory\n", __FILE__, __LINE__);
        return 1;
    }
    return 0;
}
#endif

/* sets 'flags' and variables according to the input file
argument list:
dir - the directory of the input file */
void read_input(char *input_file, Input_param *input_param)
{
    FILE *fp = fopen(input_file, "rt");
    char current_word[MAXWORD];
    float temp;
    if (!fp) {
        fprintf(stderr, "%s:%d: [ERROR] opening file: %s\n", __FILE__, __LINE__, strerror(errno));
        exit(1);
    }

    /* Seting the input variables according to the input file */
    while(fscanf(fp, "%s", current_word) != EOF) {  
        if (!strcmp(current_word, "NACA")) {
            fscanf(fp, "%d", &(input_param->NACA));
        } else if (!strcmp(current_word, "ni")) {
            fscanf(fp, "%d", &(input_param->ni));
        } else if (!strcmp(current_word, "nj")) {
            fscanf(fp, "%d ", &(input_param->nj));
        } else if (!strcmp(current_word, "num_points_on_airfoil")) {
            fscanf(fp, "%d ", &(input_param->num_points_on_airfoil));
        } else if (!strcmp(current_word, "delta_y")) {
            fscanf(fp, "%g", &temp); /* fscanf returns a float, so I cast it to double */
            input_param->delta_y = (double)temp;
        } else if (!strcmp(current_word, "XSF")) {
            fscanf(fp, "%g", &temp); /* fscanf returns a float, so I cast it to double */
            input_param->XSF = (double)temp;
        } else if (!strcmp(current_word, "YSF")) {
            fscanf(fp, "%g", &temp); /* fscanf returns a float, so I cast it to double */
            input_param->YSF = (double)temp;
        } else if (!strcmp(current_word, "r")) {
            fscanf(fp, "%g", &temp);
            input_param->r = (double)temp;
        } else if (!strcmp(current_word, "omega")) {
            fscanf(fp, "%g", &temp);
            input_param->omega = (double)temp;
        }
    }

    if (!(input_param->ni % 2)) {
        fprintf(stderr, "%s:%d: [ERROR] i_max must be even\n", __FILE__, __LINE__);
        exit(1);
    }
    if (!(input_param->num_points_on_airfoil % 2)) {
        fprintf(stderr, "%s:%d: [ERROR] num_points_on_airfoil must be odd\n", __FILE__, __LINE__);
        exit(1);
    }

    fclose(fp);
}

void output_metadata(char *output_dir, Input_param input_param)
{
    char temp_word[MAXWORD];
    FILE *metadata_file;

    strcpy(temp_word, output_dir);
    strcat(temp_word, "/metadata.txt");
    metadata_file = fopen(temp_word, "wt");
    if (!metadata_file) {
        fprintf(stderr, "%s:%d: [ERROR] opening file: %s\n", __FILE__, __LINE__, strerror(errno));
        exit(1);
    }

    fprintf(metadata_file, "%s, %s, %s, %s, %s, %s, %s, %s, %s\n", "NACA", "ni", "nj", "num_points_on_airfoil", "delta_y", "XSF", "YSF", "r", "omega");
    fprintf(metadata_file, "%d, %d, %d, %d, %f, %f, %f, %f, %f\n", input_param.NACA, input_param.ni, input_param.nj, input_param.num_points_on_airfoil, input_param.delta_y, input_param.XSF, input_param.YSF, input_param.r, input_param.omega);

    fclose(metadata_file);
}

void mat_output_to_file(FILE *fp, double *data, Input_param input_param)
{
    int i, j;
    
    for (j = 0; j < input_param.nj; j++) {
        for (i = 0; i < input_param.ni; i++) {
            fprintf(fp, "%g ", data[offset2d(i, j, input_param.ni)]);
        }
        fprintf(fp, "\n");
    }
}

void output_mesh(char *output_dir, double *x_vals_mat, double *y_vals_mat, Input_param input_param)
{
    char temp_word[MAXDIR];
    FILE *fp;
    strcpy(temp_word, output_dir);
    strcat(temp_word, "/mesh.txt");
    fp = fopen(temp_word, "wt");
    if (!fp) {
        fprintf(stderr, "%s:%d: [ERROR] opening file: %s\n", __FILE__, __LINE__, strerror(errno));
        exit(1);
    }

    fprintf(fp, "ni\n%d\n\nnj\n%d\n\n", input_param.ni, input_param.nj);
    fprintf(fp, "x_vals\n");
    mat_output_to_file(fp, x_vals_mat, input_param);
    fprintf(fp, "\n");
    fprintf(fp, "y_vals\n");
    mat_output_to_file(fp, y_vals_mat, input_param);

    fclose(fp);
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

/* saving input param to the DB and deleting if there are duplicates;
returning the return code form the 'sqlite3_exec' function. */
int output_input_param_to_DB(sqlite3 *db, Input_param input_param)
{
    char temp_sql[MAXWORD];
    sprintf(temp_sql, "insert into NACA_data(NACA, ni, nj, num_points_on_airfoil, delta_y, XSF, YSF, r, omega) values(%d, %d, %d, %d, %f, %f, %f, %f, %f);", input_param.NACA, input_param.ni, input_param.nj, input_param.num_points_on_airfoil, input_param.delta_y, input_param.XSF, input_param.YSF, input_param.r, input_param.omega);

    char *err_msg = 0;
    int rc = sqlite3_exec(db ,temp_sql, 0, 0, &err_msg);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "%s:%d: [ERROR] cannot write to DB %s\n", __FILE__, __LINE__, err_msg);
        return SQLITE_ERROR;
    }

    strcpy(temp_sql, "");
    sprintf(temp_sql, "DELETE FROM NACA_data WHERE ID NOT IN (SELECT MIN(ID) FROM NACA_data GROUP BY NACA, ni, nj, num_points_on_airfoil, delta_y, XSF, YSF, r, omega);");
    err_msg = 0;
    rc = sqlite3_exec(db ,temp_sql, 0, 0, &err_msg);
    if (rc != SQLITE_OK) {
        fprintf(stderr, "%s:%d: [ERROR] cannot write to DB %s\n", __FILE__, __LINE__, err_msg);
        return SQLITE_ERROR;
    }

    return rc;
}
