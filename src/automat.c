#include <stdio.h>
#include <sys/stat.h>
#include <errno.h>
#define __USE_MISC 
#include <dirent.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <stdint.h>
#include <unistd.h>
#include <stdarg.h>

#define dfprintINT(fp, expr) do{fprintf(fp, #expr "\n%d\n\n", expr);} while(0)     /* macro for easy debuging*/
#define dfprintD(fp, expr) do{fprintf(fp, #expr "\n%g\n\n", expr);} while(0)     /* macro for easy debuging*/
#define dfprintS(fp, expr) do{fprintf(fp, #expr "\n%s\n\n", expr);} while(0)     /* macro for easy debuging*/

int create_empty_dir(char *parent_directory);
void print_command_to_file(FILE *fp, char *program, ...);
void create_input_file(char *file_name, char *NACA, int ni, int nj, int num_points_on_airfoil, double delta_y, double XSF, double YSF, double r, double omega, double Mach_inf, double angle_of_attack_deg, double density, double environment_pressure, double delta_t, double Gamma, double epse, int max_iteration);

int main()
{
    char temp_dir[BUFSIZ], temp1[BUFSIZ], temp_input[BUFSIZ];
    int output_counter = 0;

    /* creating root output directory */
    char parent_dir[] = "./auto";
    if (create_empty_dir(parent_dir) != 0) {
        return 1;
    }

    /* opening file to print command to */
    strncpy(temp_dir, parent_dir, BUFSIZ);
    strncat(temp_dir, "/command_to_run.txt", BUFSIZ/2);
    FILE *fp = fopen(temp_dir, "wt");

    /* Makefile command to build the solver */
    fprintf(fp, "make build_and_link_main\n");

    /* creating input files and commands */
    // int NACAs[] = {5410, 12, 20};
    char *NACAs[] = {"5512", "4512", "3512", "2512", "1512", "55012", "45012", "35012", "25012", "15012", "55112", "45112", "35112", "25112", "15112"};
    // double Mach_infs[] = {0.6, 0.7, 0.8, 0.9, 1};
    // double Mach_infs[] = {1.1, 1.2 ,1.3, 1.4, 1.5};
    double Mach_infs[] = {0.8, 0.9, 1, 1.1, 1.2 ,1.3, 1.4, 1.5};


    for (int Mach_index = 4; Mach_index < 5; Mach_index++) {
        for (int NACA_index = 5; NACA_index < 7; NACA_index++) {
            for (int alphas = -8; alphas <= 8; alphas++) {
                // printf("%d\n", alphas);
                char *NACA = NACAs[NACA_index];
                int ni = 135;
                int nj = 88;
                int num_points_on_airfoil = 71;
                double delta_y = 0.003;
                double XSF = 1.05;
                double YSF = 1.05;
                double r = 0.003;
                double omega = 1.5;
                double Mach_inf = Mach_infs[Mach_index];
                double angle_of_attack_deg = alphas;
                double density = 1.225;
                double environment_pressure = 101325;
                double delta_t = 1e-5;
                double Gamma = 1.4;
                double epse = 0.06;
                int max_iteration = 3e4;

                strncpy(temp_input, parent_dir, BUFSIZ);
                strncat(temp_input, "/input", BUFSIZ/2);
                sprintf(temp1, "%d.txt", output_counter++);
                strncat(temp_input, temp1, BUFSIZ/2);
                create_input_file(temp_input, NACA, ni, nj, num_points_on_airfoil, delta_y, XSF, YSF, r, omega, Mach_inf, angle_of_attack_deg, density, environment_pressure, delta_t, Gamma, epse, max_iteration);

                strncpy(temp_dir, parent_dir, BUFSIZ);
                strncat(temp_dir, "/results", BUFSIZ/2);
                
                print_command_to_file(fp, "./build/main", temp_input, temp_dir, NULL);
            }
        }
    }
    




    fprintf(fp, "make clean_main\n");

    return 0;
}

/* if allready exists, delete all the files inside 
returns 0 on success
this function handles the errors so on fail just quit */
int create_empty_dir(char *parent_directory)
{
    char path_to_remove[BUFSIZ];

    if (mkdir(parent_directory, 0777) == -1) {
        if (errno == EEXIST) {
            DIR *dir = opendir(parent_directory);
            if (dir == NULL) {
                fprintf(stderr, "%s:%d: [Error] problem opening '%s': %s\n", __FILE__, __LINE__, parent_directory, strerror(errno));
                return 1;
            }
            struct dirent* entity;
            entity = readdir(dir);
            printf("\n");
            while (entity != NULL) {   /* directory is not empty */
                strncpy(path_to_remove, parent_directory, BUFSIZ);
                strncat(path_to_remove, "/", BUFSIZ/2);
                strncat(path_to_remove, entity->d_name, BUFSIZ);
                printf("%hhd: %s\n", entity->d_type, path_to_remove);
                if (entity->d_type == DT_REG) {
                    if (remove(path_to_remove) != 0) {
                        fprintf(stderr, "%s:%d: [Error] problem removing '%s': %s\n", __FILE__, __LINE__, path_to_remove, strerror(errno));
                        return 1;
                    }
                    printf("remove %s\n", path_to_remove);
                }
                entity = readdir(dir);
            }


            printf("\ndirectory already exist\n\n");

            closedir(dir);

            return 0;
        }

        fprintf(stderr, "%s:%d: [Error] problem making '%s': %s\n", __FILE__, __LINE__, parent_directory, strerror(errno));
        return 1;
    }
    return 0;
}

void print_command_to_file(FILE *fp, char *program, ...)
{
    fprintf(fp, "%s ", program);
    va_list args;
    va_start(args, program);
    char *temp_string;
    while ((temp_string = va_arg(args, char *)) != NULL) {
        fprintf(fp, "%s ", temp_string);
    }
    va_end(args);
    fprintf(fp, "\n");
}

void create_input_file(char *file_name, char *NACA, int ni, int nj, int num_points_on_airfoil, double delta_y, double XSF, double YSF, double r, double omega, double Mach_inf, double angle_of_attack_deg, double density, double environment_pressure, double delta_t, double Gamma, double epse, int max_iteration)
{
    FILE *input_fp = fopen(file_name, "wt");

    dfprintS(input_fp, NACA);
    dfprintINT(input_fp, ni);
    dfprintINT(input_fp, nj);
    dfprintINT(input_fp, num_points_on_airfoil);
    dfprintD(input_fp, delta_y);
    dfprintD(input_fp, XSF);
    dfprintD(input_fp, YSF);
    dfprintD(input_fp, r);
    dfprintD(input_fp, omega);
    dfprintD(input_fp, Mach_inf);
    dfprintD(input_fp, angle_of_attack_deg);
    dfprintD(input_fp, density);
    dfprintD(input_fp, environment_pressure);
    dfprintD(input_fp, delta_t);
    dfprintD(input_fp, Gamma);
    dfprintD(input_fp, epse);
    dfprintINT(input_fp, max_iteration);

    fclose(input_fp);
}
