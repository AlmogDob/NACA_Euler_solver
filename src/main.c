#include "mesher.h"

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

    double *x_mat, *y_mat;
    
    int mesh_rc = create_mesh(&x_mat, &y_mat, input_file);

    mat_print(y_mat);

    return mesh_rc;
}
