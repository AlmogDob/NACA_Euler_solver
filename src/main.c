#include <stdio.h>
#include "mesher.h"

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
#endif

void read_input(char *input_file, int *NACA, int *ni, int *nj, int *num_points_on_airfoil, double *delta_y, double *XSF, double *YSF, double *r, double *omega);
void output_metadata(char *output_dir, int NACA, int ni, int nj, int num_points_on_airfoil, double delta_y, double XSF, double YSF, double r, double omega);
void mat_output_to_file(FILE *fp, double *data, int ni, int nj);
void output_mesh(char *output_dir, double *x_vals_mat, double *y_vals_mat, int ni, int nj);

int main(int argc, char const *argv[])
{
    /* Geting the input file and output directory */
    char input_file[MAXDIR], output_dir[MAXDIR], temp_word[MAXDIR]; 

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
    double delta_y, XSF, YSF, r, omega;
    int NACA, ni, nj, num_points_on_airfoil;

    read_input(input_file, &NACA, &ni, &nj, &num_points_on_airfoil, &delta_y, &XSF, &YSF, &r, &omega);

    /* Checking that I got the right input */
    dprintINT(NACA);
    dprintINT(ni);
    dprintINT(nj);
    dprintINT(num_points_on_airfoil);
    dprintD(delta_y);
    dprintD(XSF);
    dprintD(YSF);
    dprintD(r);
    dprintD(omega);
    printf("--------------------\n");

    /* creating output directory */
    printf("[INFO] creating output directory\n");
    if (ON_LINUX) {
        if (create_empty_dir(output_dir) != 0) {
            fprintf(stderr, "%s:%d: [ERROR] creating ouput directory\n", __FILE__, __LINE__);
            return 1;
        }
        sprintf(temp_word, "/NACA%d", NACA);
        strcat(output_dir, temp_word);
        if (create_empty_dir(output_dir) != 0) {
            fprintf(stderr, "%s:%d: [ERROR] creating ouput directory\n", __FILE__, __LINE__);
            return 1;
        }
        sprintf(temp_word, "/ni%d", ni);
        strcat(output_dir, temp_word);
        if (create_empty_dir(output_dir) != 0) {
            fprintf(stderr, "%s:%d: [ERROR] creating ouput directory\n", __FILE__, __LINE__);
            return 1;
        }
        sprintf(temp_word, "/nj%d", nj);
        strcat(output_dir, temp_word);
        if (create_empty_dir(output_dir) != 0) {
            fprintf(stderr, "%s:%d: [ERROR] creating ouput directory\n", __FILE__, __LINE__);
            return 1;
        }
        sprintf(temp_word, "/num_points_on_airfoil%d", num_points_on_airfoil);
        strcat(output_dir, temp_word);
        if (create_empty_dir(output_dir) != 0) {
            fprintf(stderr, "%s:%d: [ERROR] creating ouput directory\n", __FILE__, __LINE__);
            return 1;
        }
        sprintf(temp_word, "/delta_y%g", delta_y);
        strcat(output_dir, temp_word);
        if (create_empty_dir(output_dir) != 0) {
            fprintf(stderr, "%s:%d: [ERROR] creating ouput directory\n", __FILE__, __LINE__);
            return 1;
        }
        sprintf(temp_word, "/XSF%g", XSF);
        strcat(output_dir, temp_word);
        if (create_empty_dir(output_dir) != 0) {
            fprintf(stderr, "%s:%d: [ERROR] creating ouput directory\n", __FILE__, __LINE__);
            return 1;
        }
        sprintf(temp_word, "/YSF%g", YSF);
        strcat(output_dir, temp_word);
        if (create_empty_dir(output_dir) != 0) {
            fprintf(stderr, "%s:%d: [ERROR] creating ouput directory\n", __FILE__, __LINE__);
            return 1;
        }
        sprintf(temp_word, "/r%g", r);
        strcat(output_dir, temp_word);
        if (create_empty_dir(output_dir) != 0) {
            fprintf(stderr, "%s:%d: [ERROR] creating ouput directory\n", __FILE__, __LINE__);
            return 1;
        }
        sprintf(temp_word, "/omega%g", omega);
        strcat(output_dir, temp_word);
        if (create_empty_dir(output_dir) != 0) {
            fprintf(stderr, "%s:%d: [ERROR] creating ouput directory\n", __FILE__, __LINE__);
            return 1;
        }
    }

    /* creating mesh */
    printf("[INFO] Meshing\n");
    double *x_mat, *y_mat;

    int mesh_rc = create_mesh(&x_mat, &y_mat, NACA, ni, nj, num_points_on_airfoil, delta_y, XSF, YSF, r, omega, output_dir);
    if (mesh_rc != 0) {
        fprintf(stderr, "%s:%d: [ERROR] creating mesh\n", __FILE__, __LINE__);
        return 1;
    }
    output_mesh(output_dir, x_mat, y_mat, ni, nj);

    output_metadata(output_dir, NACA, ni, nj, num_points_on_airfoil, delta_y, XSF, YSF, r, omega);

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
#endif

/* sets 'flags' and variables according to the input file
argument list:
dir - the directory of the input file */
void read_input(char *input_file, int *NACA, int *ni, int *nj, int *num_points_on_airfoil, double *delta_y, double *XSF, double *YSF, double *r, double *omega)
{
    FILE *fp = fopen(input_file, "rt");
    char current_word[MAXWORD];
    float temp;
    if (!fp) {
        fprintf(stderr, "%s:%d: [ERROR] opening file: %s\n", __FILE__, __LINE__, strerror(errno));
        exit(1);
    }

    /* Seting the input varibles according to the input file */
    while(fscanf(fp, "%s", current_word) != EOF) {  
        if (!strcmp(current_word, "NACA")) {
            fscanf(fp, "%d", NACA);
        } else if (!strcmp(current_word, "ni")) {
            fscanf(fp, "%d", ni);
        } else if (!strcmp(current_word, "nj")) {
            fscanf(fp, "%d ", nj);
        } else if (!strcmp(current_word, "num_points_on_airfoil")) {
            fscanf(fp, "%d ", num_points_on_airfoil);
        } else if (!strcmp(current_word, "delta_y")) {
            fscanf(fp, "%g", &temp); /* fscanf returns a float, so I cast it to double */
            *delta_y = (double)temp;
        } else if (!strcmp(current_word, "XSF")) {
            fscanf(fp, "%g", &temp); /* fscanf returns a float, so I cast it to double */
            *XSF = (double)temp;
        } else if (!strcmp(current_word, "YSF")) {
            fscanf(fp, "%g", &temp); /* fscanf returns a float, so I cast it to double */
            *YSF = (double)temp;
        } else if (!strcmp(current_word, "r")) {
            fscanf(fp, "%g", &temp);
            *r = (double)temp;
        } else if (!strcmp(current_word, "omega")) {
            fscanf(fp, "%g", &temp);
            *omega = (double)temp;
        }
    }

    if (!(*ni % 2)) {
        fprintf(stderr, "%s:%d: [ERROR] i_max must be even\n", __FILE__, __LINE__);
        exit(1);
    }
    if (!(*num_points_on_airfoil % 2)) {
        fprintf(stderr, "%s:%d: [ERROR] num_points_on_airfoil must be odd\n", __FILE__, __LINE__);
        exit(1);
    }

    fclose(fp);
}

void output_metadata(char *output_dir, int NACA, int ni, int nj, int num_points_on_airfoil, double delta_y, double XSF, double YSF, double r, double omega)
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
    fprintf(metadata_file, "%d, %d, %d, %d, %f, %f, %f, %f, %f\n", NACA, ni, nj, num_points_on_airfoil, delta_y, XSF, YSF, r, omega);

    fclose(metadata_file);
}

void mat_output_to_file(FILE *fp, double *data, int ni, int nj)
{
    int i, j;
    
    for (j = 0; j < nj; j++) {
        for (i = 0; i < ni; i++) {
            fprintf(fp, "%g ", data[offset2d(i, j, ni)]);
        }
        fprintf(fp, "\n");
    }
}

void output_mesh(char *output_dir, double *x_vals_mat, double *y_vals_mat, int ni, int nj)
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

    fprintf(fp, "ni\n%d\n\nnj\n%d\n\n", ni, nj);
    fprintf(fp, "x_vals\n");
    mat_output_to_file(fp, x_vals_mat, ni, nj);
    fprintf(fp, "\n");
    fprintf(fp, "y_vals\n");
    mat_output_to_file(fp, y_vals_mat, ni, nj);

    fclose(fp);
}
