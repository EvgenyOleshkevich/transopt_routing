#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <random>


using matrix = std::vector<std::vector<double>>;
using sorted_matrix = std::vector<std::vector<std::pair<double, size_t>>>;
using int_matrix = std::vector<std::vector<size_t>>;

namespace utils
{
    double sqr(double);

    matrix fill_matrix(const double const*, const double const*, const size_t);

    std::pair < matrix, sorted_matrix> fill_matrix_and_sort(const double const*, const double const*, const size_t);

    double length_rout(const std::vector<size_t>&, const matrix&);

    double length_routs(const int_matrix&, const matrix&);

}