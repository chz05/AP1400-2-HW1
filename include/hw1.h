#ifndef AP_HW1_H
#define AP_HW1_H

#include <vector>

using Matrix = std::vector<std::vector<double>>;

//the namespace is good for using similar function names or other similar things in different namespace
//we can use them in algebra::function_name.
namespace algebra{
    Matrix zeros(std::size_t n, std::size_t m);
    Matrix ones(std::size_t n, std::size_t m);
    Matrix random(std::size_t n, std::size_t m, double min, double max);
    void show(const Matrix& matrix);
    Matrix multiply(const Matrix& matrix, double c);
    Matrix multiply(const Matrix& matrix1, const Matrix& matrix2);
}

#endif //AP_HW1_H
