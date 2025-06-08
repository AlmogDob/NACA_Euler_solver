/* This one-file library is heavily inspired by Tsoding's nn.h implementation of matrix
creation and operation. you can find the source code in:
https://github.com/tsoding/nn.h .
featured in this video of his:
https://youtu.be/L1TbWe8bVOc?list=PLpM-Dvs8t0VZPZKggcql-MmjaBdZKeDMw .*/


#ifndef MATRIX_H_
#define MATRIX_H_

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef MATRIX_MALLOC
#define MATRIX_MALLOC malloc
#endif //MATRIX_MALLOC

#ifndef MATRIX_ASSERT
#include <assert.h>
#define MATRIX_ASSERT assert
#endif //MATRIX_ASSERT

typedef struct {
    size_t rows;
    size_t cols;
    size_t stride;
    double *elements;
    
} Mat;

#define MAT_AT(m, i, j) (m).elements[(i)*(m).cols + (j)]
#define MAT_PRINT(m) mat_print(m, #m, 0)

double rand_double(void);
Mat mat_alloc(size_t rows, size_t cols);
void mat_fill(Mat m, double x);
void mat_rand(Mat m, double low, double high);
void mat_dot(Mat dst, Mat a, Mat b);
void mat_sum(Mat dst, Mat a);
void mat_mult(Mat m, size_t factor);
void mat_print(Mat m, const char *name, size_t padding);
void mat_mat_print_to_file(FILE *fp, Mat m, const char *name);
void mat_identity_mat(Mat m);

#endif // MATRIX_H_

#ifdef MATRIX_IMPLEMENTATION

double rand_double(void)
{
    return (double) rand() / (double) RAND_MAX;
}

Mat mat_alloc(size_t rows, size_t cols)
{
    Mat m;
    m.rows = rows;
    m.cols = cols;
    m.stride = cols;
    m.elements = (double*)MATRIX_MALLOC(sizeof(*m.elements)*rows*cols);
    MATRIX_ASSERT(m.elements != NULL);
    return m;    
}

void mat_fill(Mat m, double x)
{
    for (size_t i = 0; i < m.rows; ++i) {
        for (size_t j = 0; j < m.cols; ++j) {
            MAT_AT(m, i, j) = x;
        }
    }
}

void mat_rand(Mat m, double low, double high)
{
    for (size_t i = 0; i < m.rows; ++i) {
        for (size_t j = 0; j < m.cols; ++j) {
            MAT_AT(m, i, j) = rand_double()*(high - low) + low;
        }
    }
}

void mat_dot(Mat dst, Mat a, Mat b)
{
    MATRIX_ASSERT(a.cols == b.rows);
    size_t n = a.cols;
    MATRIX_ASSERT(a.rows == dst.rows);
    MATRIX_ASSERT(b.cols == dst.cols);

    for (size_t i = 0; i < dst.rows; i++) {
        for (size_t j = 0; j < dst.cols; j++) {
            for (size_t k = 0; k < n; k++) {
                MAT_AT(dst, i, j) += MAT_AT(a, i, k)*MAT_AT(b, k, j);
            }
        }
    }

}

void mat_sum(Mat dst, Mat a)
{
    MATRIX_ASSERT(dst.rows == a.rows);
    MATRIX_ASSERT(dst.cols == a.cols);
    for (size_t i = 0; i < dst.rows; ++i) {
        for (size_t j = 0; j < dst.cols; ++j) {
            MAT_AT(dst, i, j) += MAT_AT(a, i, j);
        }
    }
}

void mat_mult(Mat m, size_t factor)
{
    for (size_t i = 0; i < m.rows; ++i) {
        for (size_t j = 0; j < m.cols; ++j) {
            MAT_AT(m, i, j) = MAT_AT(m, i, j) * factor;
        }
    }
}

void mat_print(Mat m, const char *name, size_t padding)
{
    printf("%*s%s = [\n", (int) padding, "", name);
    for (size_t i = 0; i < m.rows; ++i) {
        printf("%*s    ", (int) padding, "");
        for (size_t j = 0; j < m.cols; ++j) {
            printf("%g ", MAT_AT(m, m.rows -1 -i, j));
        }
        printf("\n");
    }
    printf("%*s]\n", (int) padding, "");
}

void mat_mat_print_to_file(FILE *fp, Mat m, const char *name)
{
    for (size_t i = 0; i < m.rows; ++i) {
        for (size_t j = 0; j < m.cols; ++j) {
            fprintf(fp, "%f ", MAT_AT(m, i, j));
        }
        fprintf(fp, "\n");
    }
    name = (void *)name;
}

void mat_identity_mat(Mat m)
{
    MATRIX_ASSERT(m.cols == m.rows);
    for (size_t i = 0; i < m.rows; ++i) {
        for (size_t j = 0; j < m.cols; ++j) {
            if (i == j) {
                MAT_AT(m, i, j) = 1;
            }
            else {
                MAT_AT(m, i, j) = 0;
            }
        }
    }
}

#endif // MATRIX_IMPLEMENTATION