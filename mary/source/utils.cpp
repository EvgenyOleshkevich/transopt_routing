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
        matrix dist_mat = vector<vector<double>>(count_point, vector<double>(count_point));
        for (size_t i = 0; i < count_point; ++i)
            for (size_t j = 0; j < count_point; ++j)
                dist_mat[i][j] = std::sqrt(sqr(x[i] - x[j]) + sqr(y[i] - y[j]));
        return dist_mat;
    }

    matrix fill_matrix(const vector<double>& x, const vector<double>& y)
    {
        matrix dist_mat = vector<vector<double>>(x.size(), vector<double>(x.size()));
        for (size_t i = 0; i < x.size(); ++i)
            for (size_t j = 0; j < x.size(); ++j)
                dist_mat[i][j] = std::sqrt(sqr(x[i] - x[j]) + sqr(y[i] - y[j]));
        return dist_mat;
    }

    double mean_edge(const matrix& dist_mat)
    {
        double mean_edge = 0;
        size_t count_edge = 0;
        for (size_t i = 1; i < dist_mat.size(); ++i)
            for (size_t j = 1; j < dist_mat.size(); ++j)
            {
                ++count_edge;
                mean_edge += dist_mat[i][j];
            }
        return mean_edge /= count_edge;
    }

    std::pair<matrix, sorted_matrix> fill_matrix_and_sort(
        const vector<double>& x, const vector<double>& y)
    {
        if (x.size() != y.size())
            return {matrix(), sorted_matrix() };

        const size_t count_point = x.size();
        sorted_matrix sorted_edges(count_point, vector<std::pair<double, size_t>>(count_point));
        matrix dist_mat = vector<vector<double>>(count_point, vector<double>(count_point));
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

    sorted_matrix fill_sort_matrix( const matrix& dist_mat)
    {
        const size_t count_point = dist_mat.size();
        sorted_matrix sorted_edges(count_point, vector<std::pair<double, size_t>>(count_point));
        for (size_t i = 0; i < count_point; ++i)
        {
            for (size_t j = 0; j < count_point; ++j)
                sorted_edges[i][j] = { dist_mat[i][j] , j };
            std::sort(sorted_edges[i].begin(), sorted_edges[i].end()/*,
                [](const std::pair<double, size_t>& a, const std::pair<double, size_t>& b)
                {
                    return a.first > b.first;
                }*/);
        }
        return sorted_edges;
    }

    matrix fill_matrix_with_end_point(const vector<double>& x, const vector<double>& y)
    {
        if (x.size() != y.size())
            return matrix();

        matrix dist_mat = vector<vector<double>>(x.size(), vector<double>(x.size()));
        for (size_t i = 0; i < x.size(); ++i)
            for (size_t j = 0; j < x.size(); ++j)
                dist_mat[i][j] = std::sqrt(sqr(x[i] - x[j]) + sqr(y[i] - y[j]));

        for (size_t i = 1; i < x.size(); ++i)
        {
            double extra_dist = 10;
            double d1 = std::sqrt(sqr(x[i] - (x[0] + extra_dist)) + sqr(y[i] - (y[0] + extra_dist)));
            double d2 = std::sqrt(sqr(x[i] - (x[0] - extra_dist)) + sqr(y[i] - (y[0] + extra_dist)));
            double d3 = std::sqrt(sqr(x[i] - (x[0] + extra_dist)) + sqr(y[i] - (y[0] - extra_dist)));
            double d4 = std::sqrt(sqr(x[i] - (x[0] - extra_dist)) + sqr(y[i] - (y[0] - extra_dist)));

            dist_mat[i][0] = std::min(d1, std::min(d2, std::min(d3, d4)));
        }

        return dist_mat;
    }

    double length_rout(const vector<size_t>& rout, const matrix& dist_mat)
    {
        if (rout.empty())
            return 0;
        double length = dist_mat[0][rout[0]] + dist_mat[rout.back()][0];
        for (size_t i = 1; i < rout.size(); ++i)
            length += dist_mat[rout[i - 1]][rout[i]];
        return length;
    }

    double length_rout_0(const vector<size_t>& rout, const matrix& dist_mat)
    {
        double length = 0;
        for (size_t i = 1; i < rout.size(); ++i)
            length += dist_mat[rout[i - 1]][rout[i]];
        return length;
    }

    double length_rout(const vec_int_float& rout, const matrix& dist_mat)
    {
        double length = dist_mat[0][rout[0].first] + dist_mat[rout.back().first][0];
        for (size_t i = 1; i < rout.size(); ++i)
            length += dist_mat[rout[i - 1].first][rout[i].first];
        return length;
    }

    double length_rout(const vector<size_t>& rout, const matrix& dist_mat,
        const size_t beg, const size_t end)
    {
        double length = 0;
        for (size_t i = beg; i < end; ++i)
            length += dist_mat[rout[i]][rout[i + 1]];
        return length;
    }

    double length_rout_before(const vector<size_t>& rout,
        const matrix& dist_mat, const size_t end)
    {
        double length = dist_mat[0][rout[0]];
        for (size_t i = 0; i < end; ++i)
            length += dist_mat[rout[i]][rout[i + 1]];
        return length;
    }

    double length_rout_after(const vector<size_t>& rout,
        const matrix& dist_mat, const size_t beg)
    {
        double length = dist_mat[rout.back()][0];
        for (size_t i = beg + 1; i < rout.size(); ++i)
            length += dist_mat[rout[i - 1]][rout[i]];
        return length;
    }

    double length_rout_fict(const vec_int_float& rout_fict,
        const matrix& dist_mat, const size_t fict)
    {
        vector<size_t> rout;
        for (auto v : rout_fict)
            if (v.first != fict)
                rout.push_back(v.first);
        if (rout.size() == 0)
            return 0;
        return length_rout(rout, dist_mat);
    }

    double length_rout_fict(const vector<size_t>& rout_fict,
        const matrix& dist_mat, const size_t fict)
    {
        vector<size_t> rout;
        for (auto v : rout_fict)
            if (v != fict)
                rout.push_back(v);
        if (rout.size() == 0)
            return 0;
        return length_rout(rout, dist_mat);
    }

    double length_routs(const int_matrix& routs, const matrix& dist_mat)
    {
        double length = 0;
        for (const vector<size_t>& rout : routs)
            length += length_rout(rout, dist_mat);
        return length;
    }
}