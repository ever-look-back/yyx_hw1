#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{   /* 初始化矩阵元素为0 */
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
{   /* 矩阵加法 */
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
{   /* 矩阵减法 */
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
{   /* 矩阵乘法 */
    Matrix c;
    int i, j, k;

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
{   /* 矩阵的数乘 */
    int i, j;
    for (i=0; i<a.rows; i++){
        for (j=0; j<a.cols; j++){
            a.data[i][j] *= k;
        }
    }
    return a;
}

Matrix transpose_matrix(Matrix a)
{   /* 矩阵的转置 */
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
{   /* 矩阵的行列式 */
    double result = 0;
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
{   /* 逆矩阵 */
    double det;   
    Matrix matrix;
    int i, j;

    if (a.rows != a.cols){   // 不是方阵
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }
    det = det_matrix(a);   // 矩阵a的秩
    if (det == 0){   // 矩阵的秩为0，逆矩阵不存在
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }
    /* 计算伴随矩阵 */
    matrix = create_matrix(a.rows, a.cols);
    for (i=0; i<matrix.rows; i++){
        for (j=0; j<matrix.cols; j++){
            matrix.data[i][j] = Pow(-1, i+j)*det_matrix(cal_matrix(a, j, i));
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
{   /* 矩阵的秩 */
    int rank = 0;
    int i, j, k;
    double factor;
    Matrix r = create_matrix(a.rows, a.cols);   // 创建矩阵副本

    /* 复制原矩阵 */
    for (i = 0; i < a.rows; i++) {
        for (j = 0; j < a.cols; j++) {
            r.data[i][j] = a.data[i][j];
        }
    }
    /* 进行行操作以达到简化行阶梯形状 */
    for (k = 0; k < a.cols; k++) {
        int maxi = k;
        /* 寻找绝对值最大的行 */
        for (i = k + 1; i < a.rows; i++) {
            if (fabs(r.data[i][k]) > fabs(r.data[maxi][k])) {
                maxi = i;
            }
        }
        /* 如果找到非零元素，说明当前列有贡献，秩+1 */
        if (r.data[maxi][k] != 0) {
            /* 交换当前行和最大值行 */
            if (k != maxi) {
                for (j = k; j < a.cols; j++) {
                    double temp = r.data[k][j];
                    r.data[k][j] = r.data[maxi][j];
                    r.data[maxi][j] = temp;
                }
            }
            rank++;
            /* 使下方元素变为0 */
            for (i = k + 1; i < a.rows; i++) {
                factor = r.data[i][k] / r.data[k][k];
                for (j = k; j < a.cols; j++) {
                    r.data[i][j] -= r.data[k][j] * factor;
                }
            }
        }
    }
    return rank;
}

double trace_matrix(Matrix a)
{   /* 矩阵的迹：对角线元素之和（方针）*/
    double result = 0;
    int i;

    if (a.rows != a.cols){
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    for (i=0; i<a.rows; i++){
        result += a.data[i][i];
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
{   /* 返回将矩阵a的第row行和第col列删除后得到的新矩阵（row，cal从0开始计算） */
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