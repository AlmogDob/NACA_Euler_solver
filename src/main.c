#include "mesher.h"

void read_input(char *input_file, int *NACA, int *ni, int *nj, int *num_points_on_airfoil, double *delta_y, double *XSF, double *YSF, double *r, double *omega);

int main(int argc, char const *argv[])
{
    char input_file[MAXDIR]; 
    /* Geting the input and output directories */
    if (--argc != 1) {
        fprintf(stderr, "ERROR: not right usage\nUsage: main 'input file'\n");
        return -1;
    }

    strncpy(input_file, (*(++argv)), MAXDIR);

    if (input_file[MAXDIR-1] != '\0') {
        fprintf(stderr, "Error: input too long\n");
        return -1;
    }

    /* Input variables */
    double delta_y, XSF, YSF, r, omega;
    int NACA, ni, nj, num_points_on_airfoil;

    read_input(input_file, &NACA, &ni, &nj, &num_points_on_airfoil, &delta_y, &XSF, &YSF, &r, &omega);

    double *x_mat, *y_mat;
    
    int mesh_rc = create_mesh(&x_mat, &y_mat, NACA, ni, nj, num_points_on_airfoil, delta_y, XSF, YSF, r, omega);

    mat_print(y_mat);

    return mesh_rc;
}

/* sets 'flags' and variables according to the input file
argument list:
dir - the directory of the input file */
void read_input(char *input_file, int *NACA, int *ni, int *nj, int *num_points_on_airfoil, double *delta_y, double *XSF, double *YSF, double *r, double *omega)
{
    FILE *fp = fopen(input_file, "rt");
    char current_word[MAXWORD];
    float temp;

    if (!fp) {
        fprintf(stderr, "Error opening file: %s\n", strerror(errno));
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


