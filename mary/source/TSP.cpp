#include "../headers/balancedVRP.hpp"
#include "../headers/TSP.hpp"
#include "../headers/utils.hpp"
#include <math.h>
#include <stack>

namespace TSP
{
    namespace Lin_Kernighan
    {
        std::vector<size_t> Lin_Kernighan(const double const* x, const double const* y, const size_t size)
        {
            auto p = utils::fill_matrix_and_sort(x, y, size);
            auto dist_mat = p.first;
            auto sorted_edges = p.second;

            auto rout = balancedVRP::clustering::radian_sort(x, y, size);
            std::cout << "[" << 0 << ", ";
            for (auto vertex : rout)
                std::cout << vertex << ", ";
            std::cout << 0 << "]," << std::endl;
            std::cout << "sweaping: lenght= " <<
                utils::length_rout(rout, dist_mat) << std::endl;

            auto rout2 = rout;
            local_opt::TSP_2_opt(rout2, dist_mat);
            local_opt::TSP_3_opt(rout2, dist_mat);
            local_opt::TSP_2_opt(rout2, dist_mat);
            local_opt::TSP_3_opt(rout2, dist_mat);

            std::cout << "[" << 0 << ", ";
            for (auto vertex : rout2)
                std::cout << vertex << ", ";
            std::cout << 0 << "]," << std::endl;
            std::cout << "local_opt: lenght= " <<
                utils::length_rout(rout2, dist_mat) << std::endl;

            Lin_Kernighan_by_rout(rout, dist_mat, sorted_edges);
            return rout;
        }

        void Lin_Kernighan_by_rout(std::vector<size_t>& rout,
            const matrix& dist_mat, const sorted_matrix& sorted_edges)
        {
            auto best_len = utils::length_rout(rout, dist_mat);
            std::vector<size_t> best_rout;
            std::vector<size_t> used(rout.size() + 1, 0);
            used[0] = 1;
            size_t t1 = 0;
            size_t t2 = rout[0];

            std::stack<size_t> stack_used_k;
            size_t k = 1;
            stack_used_k.push(k);
            size_t i = 0;
            do
            {
                k = stack_used_k.top();
                stack_used_k.pop();
                bool can_choose = false;
                for (; k < sorted_edges[t1].size() &&
                    sorted_edges[t1][k].first <= dist_mat[t1][t2]; ++k)
                {
                    if (1 == used[sorted_edges[t1][k].second])
                        continue;
                    size_t ty2 = sorted_edges[t1][k].second;
                    best_rout.push_back(ty2);
                    used[ty2] = 1;
                    t1 = ty2;
                    stack_used_k.push(k + 1);
                    can_choose = true;
                    break;
                }

                if (best_rout.size() == rout.size())
                {
                    double len = utils::length_rout(best_rout, dist_mat);
                    if (best_len > len)
                    {
                        best_len = len;
                        rout = best_rout;
                    }
                    used[best_rout.back()] = 0;
                    best_rout.pop_back();
                    used[best_rout.back()] = 0;
                    best_rout.pop_back();
                    stack_used_k.pop();
                    t1 = best_rout.back();
                }
                else if (can_choose)
                {
                    stack_used_k.push(1);
                }
                else {
                    if (best_rout.size() == 0)
                        break;
                    used[best_rout.back()] = 0;
                    best_rout.pop_back();
                    if (best_rout.size() > 0)
                        t1 = best_rout.back();
                    else
                        t1 = 0;
                }

                if (0 == used[rout[best_rout.size()]])
                    t2 = rout[best_rout.size()];
                else
                    for (size_t j = 1; j < used.size(); ++j)
                        if (0 == used[j])
                        {
                            t2 = j;
                            break;
                        }
                ++i;
            } while (!stack_used_k.empty() && i < 1000000);
        }
    }

    namespace local_opt
    {
        void TSP_2_opt(std::vector<size_t>& rout, const matrix& dist_mat)
        {
            auto best_len = utils::length_rout(rout, dist_mat);
            auto best_rout = rout;
            for (int i = 0; i < rout.size() - 2; ++i)
                for (int j = i + 1; j < rout.size() - 1; ++j)
                {
                    rout = std::vector<size_t>();
                    for (int k = 0; k < i; ++k)
                        rout.push_back(best_rout[k]);
                    for (int k = j; k >= i; --k)
                        rout.push_back(best_rout[k]);
                    for (int k = j + 1; k < best_rout.size(); ++k)
                        rout.push_back(best_rout[k]);

                    auto len = utils::length_rout(rout, dist_mat);
                    if (best_len > len)
                    {
                        best_len = len;
                        best_rout = rout;
                    }
                }
            rout = best_rout;
        }

        void TSP_3_opt(std::vector<size_t>& rout, const matrix& dist_mat)
        {
            auto best_len = utils::length_rout(rout, dist_mat);
            auto best_rout = rout;
            // 3 ребра смежны
            for (size_t i = 0; i < rout.size() - 3; ++i)
            {
                std::swap(rout[i + 1], rout[i + 2]);

                auto len = utils::length_rout(rout, dist_mat);
                if (best_len > len)
                {
                    best_len = len;
                    best_rout = rout;
                }
                else
                    std::swap(rout[i + 1], rout[i + 2]);
            }

            // 2 ребра смежны, пусть вторые 2

            for (size_t i = 0; i < rout.size() - 4; ++i)
                for (size_t j = i + 2; j < rout.size() - 1; ++j)
                {
                    std::vector<size_t> rout1;
                    std::vector<size_t> rout2;

                    for (size_t k = 0; k <= i; ++k)
                    {
                        rout1.push_back(best_rout[k]);
                        rout2.push_back(best_rout[k]);
                    }

                    rout1.push_back(best_rout[j + 1]);
                    for (size_t k = i + 1; k <= j; ++k)
                        rout1.push_back(best_rout[k]);
                    for (size_t k = j + 2; k < best_rout.size(); ++k)
                        rout1.push_back(best_rout[k]);

                    for (size_t k = j; k > i; --k)
                        rout2.push_back(best_rout[k]);
                    for (size_t k = j + 1; k < best_rout.size(); ++k)
                        rout2.push_back(best_rout[k]);

                    auto len = utils::length_rout(rout1, dist_mat);
                    if (best_len > len)
                    {
                        best_len = len;
                        best_rout = rout1;
                    }

                    len = utils::length_rout(rout2, dist_mat);
                    if (best_len > len)
                    {
                        best_len = len;
                        best_rout = rout2;
                    }
                }

            // 2 ребра смежны, пусть первые 2

            for (size_t j = 0; j < rout.size() - 1; ++j)
                for (size_t i = j + 3; i < rout.size() - 4; ++i)
                {
                    std::vector<size_t> rout1;
                    std::vector<size_t> rout2;

                    for (size_t k = 0; k <= j; ++k)
                    {
                        rout1.push_back(best_rout[k]);
                        rout2.push_back(best_rout[k]);
                    }


                    for (size_t k = j + 2; k <= i; ++k)
                        rout1.push_back(best_rout[k]);
                    rout1.push_back(best_rout[j + 1]);
                    for (size_t k = i + 1; k < best_rout.size(); ++k)
                        rout1.push_back(best_rout[k]);

                    for (size_t k = i; k > j; --k)
                        rout2.push_back(best_rout[k]);
                    for (size_t k = i + 1; k < best_rout.size(); ++k)
                        rout2.push_back(best_rout[k]);

                    auto len = utils::length_rout(rout1, dist_mat);
                    if (best_len > len)
                    {
                        best_len = len;
                        best_rout = rout1;
                    }

                    len = utils::length_rout(rout2, dist_mat);
                    if (best_len > len)
                    {
                        best_len = len;
                        best_rout = rout2;
                    }
                }


            // 3 ребра не смежны

            for (size_t i = 0; i < rout.size() - 5; ++i)
                for (size_t j = i + 2; j < rout.size() - 3; ++j)
                    for (size_t k = j + 2; i < rout.size() - 1; ++i)
                    {
                        std::vector<size_t> rout1;
                        std::vector<size_t> rout2;

                        for (size_t l = 0; l <= i; ++l)
                        {
                            rout1.push_back(best_rout[l]);
                            rout2.push_back(best_rout[l]);
                        }

                        for (size_t l = k; l >= j + 1; --l)
                            rout1.push_back(best_rout[l]);
                        for (size_t l = i + 1; l <= j; ++l)
                            rout1.push_back(best_rout[l]);
                        for (size_t l = k + 1; l < best_rout.size(); ++l)
                            rout1.push_back(best_rout[l]);


                        for (size_t l = j; l >= i + 1; --l)
                            rout2.push_back(best_rout[l]);
                        for (size_t l = k; l >= j + 1; --l)
                            rout2.push_back(best_rout[l]);
                        for (size_t l = k + 1; l < best_rout.size(); ++l)
                            rout2.push_back(best_rout[l]);

                        auto len = utils::length_rout(rout1, dist_mat);
                        if (best_len > len)
                        {
                            best_len = len;
                            best_rout = rout1;
                        }

                        len = utils::length_rout(rout2, dist_mat);
                        if (best_len > len)
                        {
                            best_len = len;
                            best_rout = rout2;
                        }
                    }
            rout = best_rout;
        }
    }
}