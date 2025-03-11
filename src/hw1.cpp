#include "../include/hw1.h"
#include <random>
#include <iomanip>
#include <iostream>

namespace algebra{

std::default_random_engine re(std::random_device{}());

Matrix zeros(std::size_t n, std::size_t m){
    Matrix data(n, std::vector<double>(m, 0));
    return data;
}

Matrix ones(std::size_t n, std::size_t m){
    Matrix data(n, std::vector<double>(m, 1));
    return data;
}

double getRandomNumber(double min, double max){
   std::uniform_real_distribution<double> unif(min,max);
    double a_random_double = unif(re);
   return a_random_double;
}

Matrix random(std::size_t n, std::size_t m, double min, double max){
    if (min > max){
        throw std::logic_error("The min is greater than max");
    }
    Matrix data(n, std::vector<double>(m));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            data[i][j] = getRandomNumber(min, max);
        }
    }
    return data;
}

void show(const Matrix& matrix){
    int n = matrix.size();
    int m = matrix[0].size();
    for (int i = 0; i < n; i++){
        std::cout << "| ";
        for (int j = 0; j < m; j++){
            std::cout << std::setw(8) << std::setprecision(3) <<  matrix[i][j];
        }
        std::cout << " |" << std::endl;
    }
}

Matrix multiply(const Matrix& matrix, double c){
    int n = matrix.size();
    int m = matrix[0].size();
    Matrix mul_scalar_matrix = matrix;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            mul_scalar_matrix[i][j] = matrix[i][j] * c;
        }
    }
    return mul_scalar_matrix;
}

Matrix multiply(const Matrix& matrix1, const Matrix& matrix2){
    if (matrix1.empty() || matrix2.empty()){
        return Matrix{};
    }
    int n1 = matrix1.size();
    int m1 = matrix1[0].size();
    int n2 = matrix2.size();
    int m2 = matrix2[0].size();
    if (m1 != n2){
        throw std::logic_error("The dimensions are incorrect for these two matrices");
    }
    //multiple
    Matrix matrix3(n1, std::vector<double>(m2));
    for (int i = 0; i < n1; i++){
        for (int j = 0; j < m2; j++){
            double temp = 0;
            for (int k = 0; k < m1; k++){
                temp += matrix1[i][k] * matrix2[k][j];
            }
            matrix3[i][j] = temp;
        }
    }
    return matrix3;

}

Matrix sum(const Matrix& matrix1, double c){
    if (matrix1.empty()){
        return Matrix{};
    } else{
        Matrix sum_matrix = matrix1;
        for (int i = 0; i < matrix1.size(); i++){
            for (int j = 0; j < matrix1[0].size(); j++){
                sum_matrix[i][j] += c;
            }
        }
        return sum_matrix;
    }
}


Matrix sum(const Matrix& matrix1, const Matrix& matrix2){
    if (matrix1.empty() && matrix2.empty()){
        return Matrix{};
    }

    if (matrix1.empty() || matrix2.empty()){
        throw std::logic_error("one of matrix is empty.");
    }

    int n1 = matrix1.size();
    int m1 = matrix1[0].size();
    int n2 = matrix2.size();
    int m2 = matrix2[0].size();
    if (n1 != n2 || m1 != m2){
        throw std::logic_error("The dimension of these two matrices are incorrect.");
    }
    //after that we check empty again
    
    Matrix sum_matrix (n1, std::vector<double>(m1));
    for (int i = 0; i < n1; i++){
        for (int j = 0; j < m1; j++){
            sum_matrix[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
    return sum_matrix;
}

Matrix transpose(const Matrix& matrix){
    if (matrix.empty()){
        return Matrix{};
    }
    int n = matrix.size();
    int m = matrix[0].size();
    Matrix transpose_matrix (m, std::vector<double>(n));
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            transpose_matrix[i][j] = matrix[j][i];
        }
    }
    return transpose_matrix;
}

Matrix minor(const Matrix& matrix, std::size_t n, std::size_t m){
    if (matrix.empty()){
        return Matrix{};
    }
    int row = matrix.size();
    int column = matrix[0].size();
    if (n >= row || n < 0 || m >= column || m < 0){
        throw std::logic_error("The nth row or mth column is not existed.");
    }
    Matrix minor_matrix (row - 1, std::vector<double>(column - 1));
    for (int i = 0; i < row; i++){
        for (int j = 0; j < column; j++){
            if (i < n && j < m){
                minor_matrix[i][j] = matrix[i][j];
            } else if (i < n && j > m){
                minor_matrix[i][j-1] = matrix[i][j];
            } else if (i > n && j < m){
                minor_matrix[i-1][j] = matrix[i][j];
            } else if (i > n && j > m){
                minor_matrix[i-1][j-1] = matrix[i][j];
            }
        }
    }
    return minor_matrix;

}

double determinant(const Matrix& matrix){
    if (matrix.empty()){
        return 1;
    }
    if (matrix.size() == 1 && matrix[0].size() == 1){
        return matrix[0][0];
    }
    int n = matrix.size();
    int m = matrix[0].size();
    if (n != m){
        throw std::logic_error("The number of rows is not equal to the number columns.");
    }
    if (n == 2 && m == 2){
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    }
    // 3 * 3 matrix or more
    // fix the row is 0
    double det = 0;
    for (int j = 0; j < m; j++){
        Matrix minor_matrix = minor(matrix, 0, j);
        det += (j % 2 == 0 ? 1 : -1) * matrix[0][j] * determinant(minor_matrix);
    }
    return det;
}

Matrix inverse(const Matrix& matrix){
    if (matrix.empty() || matrix[0].empty()){
        return Matrix{};
    }
    int n = matrix.size();
    int m = matrix[0].size();
    if (n != m){
        throw std::logic_error("It is not a squared matrix.");
    }
    if (determinant(matrix) == 0){
        throw std::logic_error("This is singular matrix.");
    }
    // Matrix of minors
    Matrix matrix_minors (n, std::vector<double>(m));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            matrix_minors[i][j] = determinant(minor(matrix, i, j));
        }
    }
    // matrix of cofactors
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            if ((i + j) % 2 == 1){
                matrix_minors[i][j] = -matrix_minors[i][j];
            }
        }
    }
    //adjugate
    Matrix adjugate_matrix = transpose(matrix_minors);
    Matrix inverse_matrix = multiply(adjugate_matrix, 1/determinant(matrix));
    return inverse_matrix;
}

Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis){
    if (matrix1.empty() && matrix2.empty()){
        return Matrix{};
    }
    if (matrix1.empty()){
        return matrix2;
    }
    if (matrix2.empty()){
        return matrix1;
    }
    int n1 = matrix1.size();
    int m1 = matrix1[0].size();
    int n2 = matrix2.size();
    int m2 = matrix2[0].size();
    if (axis == 0){
        if (m1 != m2){
            throw std::logic_error("The dimension of columns are different");
        }
        Matrix matrix3 (n1+n2, std::vector<double>(m1));
        for (int i = 0; i < n1+n2; i++){
            for (int j = 0; j < m1; j++){
                if (i >= n1){
                    matrix3[i][j] = matrix2[i-n1][j];
                } else{
                    matrix3[i][j] = matrix1[i][j];
                }
            }
        }
        return matrix3;
    }else if (axis == 1){
        if (n1 != n2){
            throw std::logic_error("The dimension of rows are different");
        }
        Matrix matrix3 (n1, std::vector<double>(m1+m2));
        for (int i = 0; i < n1; i++){
            for (int j = 0; j < m1+m2; j++){
                if (j >= m1){
                    matrix3[i][j] = matrix2[i][j-m1];
                } else{
                    matrix3[i][j] = matrix1[i][j];
                }
            }
        }
        return matrix3;
    }else{
        throw std::logic_error("axis is not 0 or 1.");
    }
}

Matrix ero_swap(const Matrix& matrix, std::size_t r1, std::size_t r2){
    if (matrix.empty()){
        return Matrix{};
    }
    int n = matrix.size();
    int m = matrix[0].size();
    if (r1 < 0 || r1 >= n || r2 < 0 || r2 >= n){
        throw std::logic_error("the row is out of range.");
    }
    if (r1 == r2){
        return matrix;
    }
    Matrix swap_matrix = matrix;
    std::vector<double> temp(m);
    for (int j = 0; j < m; j++){
        temp[j] = matrix[r1][j];
    }
    for (int j = 0; j < m; j++){
        swap_matrix[r1][j] = matrix[r2][j];
        swap_matrix[r2][j] = temp[j];
    }
    return swap_matrix;
}

Matrix ero_multiply(const Matrix& matrix, std::size_t r, double c){
    if (matrix.empty()){
        return Matrix{};
    }
    int n = matrix.size();
    int m = matrix[0].size();
    if (r < 0 || r >= n){
        throw std::logic_error("the row is out of range.");
    }
    Matrix matrix2 = matrix;
    for (int j = 0; j < m; j++){
        matrix2[r][j] = matrix[r][j] * c;
    }
    return matrix2;
}

Matrix ero_sum(const Matrix& matrix, std::size_t r1, double c, std::size_t r2){
    if (matrix.empty()){
        return Matrix{};
    }
    int n = matrix.size();
    int m = matrix[0].size();
    if (r1 < 0 || r1 >= n || r2 < 0 || r2 >= n){
        throw std::logic_error("the row is out of range.");
    }
    std::vector<double> temp(m);
    for (int j = 0; j < m; j++){
        temp[j] = matrix[r1][j] * c;
    }
    Matrix matrix2 = matrix;
    for (int j = 0; j < m; j++){
        matrix2[r2][j] += temp[j];
    }
    return matrix2;
}

Matrix upper_triangular(const Matrix& matrix){
    if (matrix.empty()){
        return Matrix{};
    }
    int n = matrix.size();
    int m = matrix[0].size();
    if (n != m){
        throw std::logic_error("It is non-square matrix.");
    }
    Matrix matrix2 = matrix;

    //if the matrix2[i][i] is 0, then we can swap with the below rows
    for (int i = 0; i < n-1; i++){
        if (matrix2[i][i] == 0){
            for (int j = i+1; j < n; j++){
                matrix2 = ero_swap(matrix2, i, j);
                if (matrix2[i][i] != 0){
                    break;
                }
            }   
        }
    }

    for (int i = 0; i < n-1; i++){
        
        for (int j = i+1; j < m; j++){
            matrix2 = ero_sum(matrix2, i, -matrix2[j][i]/matrix2[i][i],j);
        }
    }
    return matrix2;
}


}
