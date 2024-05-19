#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    int i, j;    
    m.rows = row;
    m.cols = col;
    for (i = 0; i < row; i++){
        for (j = 0; j< col; j++){
            m.data[i][j] = 0;
        }
    }
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    Matrix c;
    int i, j;

    if (a.rows!=b.rows || a.cols!=b.cols){
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    c = create_matrix(a.rows, a.cols);
    for (i=0; i<c.rows; i++){
        for (j=0; j<c.rows; j++){
            c.data[i][j] = a.data[i][j] + b.data[i][j];
        }
    }
    return c;
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    Matrix c;
    int i, j;

    if (a.rows!=b.rows || a.cols!=b.cols){
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    c = create_matrix(a.rows, a.cols);
    for (i=0; i<c.rows; i++){
        for (j=0; j<c.rows; j++){
            c.data[i][j] = a.data[i][j] - b.data[i][j];
        }
    }
    return c;
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    Matrix c;
    int i, j ,k;

    if (a.cols != b.rows){
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);
    }
    c = create_matrix(a.rows, b.cols);
    for (i=0; i<c.rows; i++){
        for (j=0; j<a.cols; j++){
            for (k=0; k<a.cols; k++){
                c.data[i][j] += a.data[i][k]*b.data[k][j];
            }
        }
    }
    return c;
}

Matrix scale_matrix(Matrix a, double k)
{
    int i, j;
    for (i=0; i<a.rows; i++){
        for (j=0; j<a.cols; j++){
            a.data[i][j] *= k;
        }
    }
    return a;
}

Matrix transpose_matrix(Matrix a)
{
    int i, j, temp;
    Matrix b;
    b = create_matrix(a.cols, a.rows);
    for (i = 0; i < b.rows; i++){
        for (j = 0; j < b.cols; j++){
            b.data[i][j] = a.data[j][i];
        }
    }
    return b;
}

double det_matrix(Matrix a)
{
    double result = 1;
    int i;

    if (a.rows != a.cols){
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    if (a.rows == 1){
        return a.data[0][0];
    }
    for (i=0; i<a.rows; i++){
        result += Pow(-1,i)*a.data[0][i]*det_matrix(cal_matrix(a, 0, i));   // Leplace定理
    }
    return result;
}

Matrix inv_matrix(Matrix a)
{
    double det;   // 矩阵a的秩
    Matrix matrix;
    int i, j;

    if (a.rows != a.cols){   // 不是方阵
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }
    det = det_matrix(a);
    if (det == 0){   // 矩阵的秩为0，逆矩阵不存在
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }
    /* 计算伴随矩阵 */
    matrix = create_matrix(a.rows, a.cols);
    for (i=0; i<matrix.rows; i++){
        for (j=0; j<matrix.cols; j++){
            matrix.data[i][j] = Pow(-1, i+j)*det_matrix(cal_matrix(a, i, j));
        }
    }
    /* 计算逆矩阵并返回 */
    for (i=0; i<matrix.rows; i++){
        for (j=0; j<matrix.cols; j++){
            matrix.data[i][j] /= det; 
        }
    }
    return matrix;
}

int rank_matrix(Matrix a)
{
    // ToDo
    return 0;
}

double trace_matrix(Matrix a)
{
    double result = 1;
    int i;

    if (a.rows != a.cols){
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    for (i=0; i<a.rows; i++){
        result *= a.data[i][i];
    }
    return result;
}

void print_matrix(Matrix a)
{
    int i, j;
    for (i = 0; i < a.rows; i++)
    {
        for (j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}

double Pow(int a, int n)
{
    int i;
    double result = 1;
    for (i=0; i<n; i++){
        result *= a;
    }
    return result;
}

Matrix cal_matrix(Matrix a, int row, int col)
{
    Matrix b;
    int i, j;

    b = create_matrix(a.rows-1, a.cols-1);
    for (i=0; i<b.rows; i++){
        for (j=0; j<b.cols; j++){
            if (i>=row && j>=col){
                b.data[i][j] = a.data[i+1][j+1];
            }else if (i>=row && j<col){
                b.data[i][j] = a.data[i+1][j];
            }else if (i<row && j>=col){
                b.data[i][j] = a.data[i][j+1];
            }else{
                b.data[i][j] = a.data[i][j];
            }
        }
    }
    return b;
}