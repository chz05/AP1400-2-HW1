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


}
