#ifndef UTILS
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <random>


using matrix = std::vector<std::vector<double>>;
using sorted_matrix = std::vector<std::vector<std::pair<double, size_t>>>;
using int_matrix = std::vector<std::vector<size_t>>;
using vec_int_float = std::vector<std::pair<size_t, double>>;

namespace utils
{
    double sqr(double);

    matrix fill_matrix(const double* const, const double* const, const size_t);

    std::pair < matrix, sorted_matrix> fill_matrix_and_sort(const double* const, const double* const, const size_t);

    double length_rout(const std::vector<size_t>&, const matrix&);

    double length_rout(const vec_int_float&, const matrix&);

    double length_rout(const std::vector<size_t>&, const matrix&, const size_t, const size_t);

    double length_rout_before(const std::vector<size_t>&, const matrix&, const size_t);

    double length_rout_after(const std::vector<size_t>&, const matrix&, const size_t);

    double length_rout_fict(const vec_int_float&, const matrix&, const size_t);

    double length_rout_fict(const std::vector<size_t>&, const matrix&, const size_t);

    double length_routs(const int_matrix&, const matrix&);

}
#endif
#define UTILS