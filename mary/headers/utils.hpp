#ifndef UTILS
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <random>
#include <cmath>
#include <stdexcept>
#include <chrono> // для std::random_device и std::mt19937


using std::vector;
using std::pair;
using matrix = vector<vector<double>>;
using sorted_matrix = vector<vector<pair<double, size_t>>>;
using int_matrix = vector<vector<size_t>>;
using vec_int_float = vector<pair<size_t, double>>;
//template <typename type> using matrix = vector < vector <type>>;
constexpr double EPS = 0.000001;


namespace utils
{
    double sqr(double);

    matrix fill_matrix(const double* const, const double* const, const size_t);

    matrix fill_matrix(const vector<double>&, const vector<double>&, const size_t);

    pair < matrix, sorted_matrix> fill_matrix_and_sort(const vector<double>&, const vector<double>&);

    sorted_matrix fill_matrix_and_sort(const vector<double>&, const vector<double>&, const matrix&);

    matrix fill_matrix_with_end_point(const vector<double>&, const vector<double>&);

    double length_rout(const vector<size_t>&, const matrix&);

    double length_rout_0(const vector<size_t>&, const matrix&);

    double length_rout(const vec_int_float&, const matrix&);

    double length_rout(const vector<size_t>&, const matrix&, const size_t, const size_t);

    double length_rout_before(const vector<size_t>&, const matrix&, const size_t);

    double length_rout_after(const vector<size_t>&, const matrix&, const size_t);

    double length_rout_fict(const vec_int_float&, const matrix&, const size_t);

    double length_rout_fict(const vector<size_t>&, const matrix&, const size_t);

    double length_routs(const int_matrix&, const matrix&);

}
#endif
#define UTILS