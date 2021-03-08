#include "../headers/balancedVRP.hpp"
#include "../headers/utils.hpp"
#include <random>
#include <math.h>

namespace utils
{
    double sqr(double x)
    {
        return x * x;
    }

    matrix fill_matrix(const double* const x, const double* const y,
        const size_t count_point)
    {
        matrix dist_mat = std::vector<std::vector<double>>(count_point, std::vector<double>(count_point));
        for (size_t i = 0; i < count_point; ++i)
            for (size_t j = 0; j < count_point; ++j)
                dist_mat[i][j] = std::sqrt(sqr(x[i] - x[j]) + sqr(y[i] - y[j]));
        return dist_mat;
    }

    std::pair<matrix, sorted_matrix> fill_matrix_and_sort(
        const double* const x, const double* const y,
        const size_t count_point)
    {
        sorted_matrix sorted_edges(count_point, std::vector<std::pair<double, size_t>>(count_point));
        matrix dist_mat = std::vector<std::vector<double>>(count_point, std::vector<double>(count_point));
        for (size_t i = 0; i < count_point; ++i)
        {
            for (size_t j = 0; j < count_point; ++j)
            {
                dist_mat[i][j] = std::sqrt(sqr(x[i] - x[j]) + sqr(y[i] - y[j]));
                sorted_edges[i][j] = { dist_mat[i][j] , j };
            }
            std::sort(sorted_edges[i].begin(), sorted_edges[i].end()/*,
                [](const std::pair<double, size_t>& a, const std::pair<double, size_t>& b)
                {
                    return a.first > b.first;
                }*/);
        }
        return { dist_mat, sorted_edges };
    }

    double length_rout(const std::vector<size_t>& rout, const matrix& dist_mat)
    {
        double lenght = dist_mat[0][rout[0]] + dist_mat[rout.back()][0];
        for (size_t i = 1; i < rout.size(); ++i)
            lenght += dist_mat[rout[i - 1]][rout[i]];
        return lenght;
    }

    double length_rout(const vec_int_float& rout, const matrix& dist_mat)
    {
        double lenght = dist_mat[0][rout[0].first] + dist_mat[rout.back().first][0];
        for (size_t i = 1; i < rout.size(); ++i)
            lenght += dist_mat[rout[i - 1].first][rout[i].first];
        return lenght;
    }

    double length_rout(const std::vector<size_t>& rout, const matrix& dist_mat,
        const size_t beg, const size_t end)
    {
        double lenght = 0;
        for (size_t i = beg; i < end; ++i)
            lenght += dist_mat[rout[i]][rout[i + 1]];
        return lenght;
    }

    double length_rout_before(const std::vector<size_t>& rout,
        const matrix& dist_mat, const size_t end)
    {
        double lenght = dist_mat[0][rout[0]];
        for (size_t i = 0; i < end; ++i)
            lenght += dist_mat[rout[i]][rout[i + 1]];
        return lenght;
    }

    double length_rout_after(const std::vector<size_t>& rout,
        const matrix& dist_mat, const size_t beg)
    {
        double lenght = dist_mat[rout.back()][0];
        for (size_t i = beg + 1; i < rout.size(); ++i)
            lenght += dist_mat[rout[i - 1]][rout[i]];
        return lenght;
    }

    double length_rout_fict(const vec_int_float& rout_fict,
        const matrix& dist_mat, const size_t fict)
    {
        std::vector<size_t> rout;
        for (auto v : rout_fict)
            if (v.first != fict)
                rout.push_back(v.first);
        if (rout.size() == 0)
            return 0;
        return length_rout(rout, dist_mat);
    }

    double length_rout_fict(const std::vector<size_t>& rout_fict,
        const matrix& dist_mat, const size_t fict)
    {
        std::vector<size_t> rout;
        for (auto v : rout_fict)
            if (v != fict)
                rout.push_back(v);
        if (rout.size() == 0)
            return 0;
        return length_rout(rout, dist_mat);
    }

    double length_routs(const int_matrix& routs, const matrix& dist_mat)
    {
        double lenght = 0;
        for (const std::vector<size_t>& rout : routs)
            lenght += length_rout(rout, dist_mat);
        return lenght;
    }
}