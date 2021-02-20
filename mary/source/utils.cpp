#include "../headers/balancedVRP.hpp"
#include "../headers/utils.hpp"
#include <math.h>

namespace utils
{
    double sqr(double x)
    {
        return x * x;
    }

    matrix fill_matrix(const double const* x, const double const* y,
        const size_t count_point)
    {
        matrix dist_mat = std::vector<std::vector<double>>(count_point, std::vector<double>(count_point));
        for (size_t i = 0; i < count_point; ++i)
            for (size_t j = 0; j < count_point; ++j)
                dist_mat[i][j] = std::sqrt(sqr(x[i] - x[j]) + sqr(y[i] - y[j]));
        return dist_mat;
    }

    double length_rout(const std::vector<size_t>& rout, const matrix& dist_mat)
    {
        double lenght = dist_mat[0][rout[0]] + dist_mat[rout.back()][0];
        for (size_t i = 1; i < rout.size(); ++i)
            lenght += dist_mat[rout[i - 1]][rout[i]];
        return lenght;
    }

    double length_routs(const int_matrix& routs, const matrix& dist_mat)
    {
        double lenght = 0;
        for (const std::vector<size_t>& rout : routs)
            lenght += length_rout(rout, dist_mat);
        return lenght;
    }


}