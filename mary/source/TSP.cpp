#include "../headers/balancedVRP.hpp"
#include "../headers/TSP.hpp"
#include "../headers/utils.hpp"
#include <math.h>
#include <stack>

namespace TSP
{
    namespace Lin_Kernighan
    {
        std::vector<size_t> Lin_Kernighan(const std::vector<double>& x, const std::vector<double>& y, const size_t size)
        {
            auto p = utils::fill_matrix_and_sort(x, y, size);
            auto dist_mat = p.first;
            auto sorted_edges = p.second;

            auto rout = balancedVRP::clustering::radian_sort(x, y, size);

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
                /*
                t2 = 0;
                for (size_t j = 0; j < rout.size() - 1; ++j)
                    if (t1 == rout[j])
                    {
                        if (0 == used[rout[j + 1]])
                            t2 = rout[j + 1];
                        break;
                    }
                if (t2 == 0)
                    for (size_t j = 1; j < used.size(); ++j)
                        if (0 == used[j])
                        {
                            t2 = j;
                            break;
                        }
                */
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
        void opt_2(std::vector<size_t>& rout, const matrix& dist_mat)
        {
            if (rout.size() < 3)
                return;
            auto best_len = utils::length_rout(rout, dist_mat);
            auto best_rout = rout;
            for (int i = 1; i < (int)rout.size() - 2; ++i)
                for (int j = i + 1; j < (int)rout.size() - 1; ++j)
                {
                    rout = std::vector<size_t>();
                    for (int k = 0; k < i; ++k)
                        rout.push_back(best_rout[k]);
                    for (int k = j; k >= i; --k)
                        rout.push_back(best_rout[k]);
                    for (size_t k = j + 1; k < best_rout.size(); ++k)
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

        void opt_2_fast(std::vector<size_t>& rout, const matrix& dist_mat)
        {
            if (rout.size() < 3)
                return;
            auto best_len = utils::length_rout(rout, dist_mat);
            auto best_rout = rout;
            for (int i = 1; i < (int)rout.size() - 2; ++i)
                for (int j = i + 1; j < (int)rout.size() - 1; ++j)
                {
                    size_t index = 0;

                    for (int k = 0; k < i; ++k)
                    {
                        rout[index] = best_rout[k];
                        ++index;
                    }
                    for (int k = j; k >= i; --k)
                    {
                        rout[index] = best_rout[k];
                        ++index;
                    }
                    for (size_t k = j + 1; k < best_rout.size(); ++k)
                    {
                        rout[index] =  best_rout[k];
                        ++index;
                    }

                    double lenght = dist_mat[0][rout[0]] + dist_mat[rout.back()][0];
                    for (size_t k = 1; k < rout.size(); ++k)
                        lenght += dist_mat[rout[k - 1]][rout[k]];
                    if (best_len > lenght)
                    {
                        best_len = lenght;
                        best_rout = rout;
                    }
                }
            rout = best_rout;
        }

        void opt_2_symmetrical(std::vector<size_t>& rout, const matrix& dist_mat)
        {
            size_t segment1 = 0;
            size_t segment2 = 0;
            double best_len = utils::length_rout(rout, dist_mat);
            double base_len = best_len;
            auto best_rout = rout;
            // intial and final position are fixed (initial/final node remains 0)
            for (int a = 1; a < rout.size() - 2; a++)
            {
                size_t A1 = rout[a - 1];
                size_t A2 = rout[a];
                for (int b = a + 1; b < rout.size() - 1; b++)
                {
                    size_t B1 = rout[b];
                    size_t B2 = rout[b + 1];
                    double len = base_len + dist_mat[A1][B1] + dist_mat[A2][B2]
                         - dist_mat[A1][A2] - dist_mat[B1][B2];
                    if (best_len > len)
                    {
                        best_len = len;
                        segment1 = a;
                        segment2 = b;
                    }
                }
            }

            if (segment2 == 0)
                return;

            size_t index = 0;

            for (int k = 0; k < segment1; ++k)
            {
                rout[index] = best_rout[k];
                ++index;
            }
            for (int k = segment2; k >= segment1; --k)
            {
                rout[index] = best_rout[k];
                ++index;
            }
            for (size_t k = segment2 + 1; k < best_rout.size(); ++k)
            {
                rout[index] = best_rout[k];
                ++index;
            }
        }

        void opt_2_fast2(std::vector<size_t>& rout, const matrix& dist_mat)
        {
            size_t segment1 = 0;
            size_t segment2 = 0;
            double best_len = utils::length_rout(rout, dist_mat);
            double base_len = best_len;
            double A_B_reverse_len = 0;
            auto best_rout = rout;
            // intial and final position are fixed (initial/final node remains 0)
            for (int i = 1; i < rout.size() - 2; i++)
            {
                size_t A1 = rout[i - 1];
                size_t A2 = rout[i];
                A_B_reverse_len = 0;
                for (int j = i + 1; j < rout.size() - 1; j++)
                {
                    
                    size_t B1 = rout[j];
                    size_t B2 = rout[j + 1];
                    A_B_reverse_len += dist_mat[B1][B1 - 1] - dist_mat[B1 - 1][B1];
                    double len = base_len + A_B_reverse_len
                        + dist_mat[A1][B1] + dist_mat[A2][B2]
                        - dist_mat[A1][A2] - dist_mat[B1][B2];
                    if (best_len > len)
                    {
                        best_len = len;
                        segment1 = i;
                        segment2 = j;
                    }
                }
            }

            if (segment2 == 0)
                return;

            size_t index = 0;

            for (int k = 0; k < segment1; ++k)
            {
                rout[index] = best_rout[k];
                ++index;
            }
            for (int k = segment2; k >= segment1; --k)
            {
                rout[index] = best_rout[k];
                ++index;
            }
            for (size_t k = segment2 + 1; k < best_rout.size(); ++k)
            {
                rout[index] = best_rout[k];
                ++index;
            }
        }

        void opt_3(std::vector<size_t>& rout, const matrix& dist_mat)
        {
            if (rout.size() < 4)
                return;
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

            for (size_t i = 0; false && i < rout.size() - 5; ++i)
                for (size_t j = i + 2; j < rout.size() - 3; ++j)
                    for (size_t k = j + 2; k < rout.size() - 1; ++k)
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

        void opt_3_fast(std::vector<size_t>& rout, const matrix& dist_mat)
        {
            if (rout.size() < 4)
                return;
            auto best_len = utils::length_rout(rout, dist_mat);
            double len = 0;
            
            // 3 ребра смежны
            for (size_t i = 0; i < rout.size() - 4; ++i)
            {
                // swap(rout[i + 1], rout[i + 2]);
                
                len = best_len
                    - dist_mat[rout[i]][rout[i + 1]]
                    - dist_mat[rout[i + 1]][rout[i + 2]]
                    - dist_mat[rout[i + 2]][rout[i + 3]]
                    + dist_mat[rout[i]][rout[i + 2]]
                    + dist_mat[rout[i + 2]][rout[i + 1]]
                    + dist_mat[rout[i + 1]][rout[i + 3]];

                if (best_len > len)
                {
                    best_len = len;
                    std::swap(rout[i + 1], rout[i + 2]);
                }
            }

            // on the ends:
            {
                len = best_len
                    - dist_mat[0][rout[0]]
                    - dist_mat[rout[0]][rout[1]]
                    - dist_mat[rout[1]][rout[2]]
                    + dist_mat[0][rout[1]]
                    + dist_mat[rout[1]][rout[0]]
                    + dist_mat[rout[0]][rout[2]];

                if (best_len > len)
                {
                    best_len = len;
                    std::swap(rout[0], rout[1]);
                }

                size_t end_index = rout.size() - 3;
                len = best_len
                    - dist_mat[rout[end_index]][rout[end_index + 1]]
                    - dist_mat[rout[end_index + 1]][rout[end_index + 2]]
                    - dist_mat[rout[end_index + 2]][0]
                    + dist_mat[rout[end_index]][rout[end_index + 2]]
                    + dist_mat[rout[end_index + 2]][rout[end_index + 1]]
                    + dist_mat[rout[end_index + 1]][0];

                if (best_len > len)
                {
                    best_len = len;
                    std::swap(rout[end_index + 1], rout[end_index + 2]);
                }
            }

            // 2 ребра смежны, пусть вторые 2
            std::vector<size_t> rout1(rout.size());
            std::vector<size_t> rout2(rout.size());
            size_t index1 = 0;
            size_t index2 = 0;
            for (size_t i = 0; i < rout.size() - 4; ++i)
                for (size_t j = i + 2; j < rout.size() - 1; ++j)
                {
                    index1 = 0;
                    index2 = 0;

                    for (size_t k = 0; k <= i; ++k)
                    {
                        rout1[index1] = rout[k];
                        rout2[index2] = rout[k];
                        ++index1;
                        ++index2;
                    }

                    rout1[index1] = rout[j + 1];
                    ++index1;
                    for (size_t k = i + 1; k <= j; ++k)
                    {
                        rout1[index1] = rout[k];
                        ++index1;
                    }
                    for (size_t k = j + 2; k < rout.size(); ++k)
                    {
                        rout1[index1] = rout[k];
                        ++index1;
                    }

                    for (size_t k = j; k > i; --k)
                    {
                        rout2[index2] = rout[k];
                        ++index2;
                    }
                    for (size_t k = j + 1; k < rout.size(); ++k)
                    {
                        rout2[index2] = rout[k];
                        ++index2;
                    }

                    len = dist_mat[0][rout1[0]] + dist_mat[rout1.back()][0];
                    for (size_t k = 1; k < rout.size(); ++k)
                        len += dist_mat[rout1[k - 1]][rout1[k]];
                    if (best_len > len)
                    {
                        best_len = len;
                        rout = rout1;
                    }

                    len = dist_mat[0][rout2[0]] + dist_mat[rout2.back()][0];
                    for (size_t k = 1; k < rout.size(); ++k)
                        len += dist_mat[rout2[k - 1]][rout2[k]];
                    if (best_len > len)
                    {
                        best_len = len;
                        rout = rout2;
                    }
                }

            // 2 ребра смежны, пусть первые 2

            for (size_t j = 0; j < rout.size() - 1; ++j)
                for (size_t i = j + 3; i < rout.size() - 4; ++i)
                {
                    index1 = 0;
                    index2 = 0;

                    for (size_t k = 0; k <= j; ++k)
                    {
                        rout1[index1] = rout[k];
                        rout2[index2] = rout[k];
                        ++index1;
                        ++index2;
                    }

                    for (size_t k = j + 2; k <= i; ++k)
                    {
                        rout1[index1] = rout[k];
                        ++index1;
                    }
                    rout1[index1] = rout[j + 1];
                    ++index1;
                    for (size_t k = i + 1; k < rout.size(); ++k)
                    {
                        rout1[index1] = rout[k];
                        ++index1;
                    }

                    for (size_t k = i; k > j; --k)
                    {
                        rout2[index2] = rout[k];
                        ++index2;
                    }
                    for (size_t k = i + 1; k < rout.size(); ++k)
                    {
                        rout2[index2] = rout[k];
                        ++index2;
                    }

                    len = dist_mat[0][rout1[0]] + dist_mat[rout1.back()][0];
                    for (size_t k = 1; k < rout.size(); ++k)
                        len += dist_mat[rout1[k - 1]][rout1[k]];
                    if (best_len > len)
                    {
                        best_len = len;
                        rout = rout1;
                    }

                    len = dist_mat[0][rout2[0]] + dist_mat[rout2.back()][0];
                    for (size_t k = 1; k < rout.size(); ++k)
                        len += dist_mat[rout2[k - 1]][rout2[k]];
                    if (best_len > len)
                    {
                        best_len = len;
                        rout = rout2;
                    }
                }


            // 3 ребра не смежны

            for (size_t i = 0; false && i < rout.size() - 5; ++i)
                for (size_t j = i + 2; j < rout.size() - 3; ++j)
                    for (size_t k = j + 2; k < rout.size() - 1; ++k)
                    {
                        index1 = 0;
                        index2 = 0;

                        for (size_t l = 0; l <= i; ++l)
                        {
                            rout1[index1] = rout[l];
                            rout2[index2] = rout[l];
                            ++index1;
                            ++index2;
                        }

                        for (size_t l = k; l >= j + 1; --l)
                        {
                            rout1[index1] = rout[l];
                            ++index1;
                        }
                        for (size_t l = i + 1; l <= j; ++l)
                        {
                            rout1[index1] = rout[l];
                            ++index1;
                        }
                        for (size_t l = k + 1; l < rout.size(); ++l)
                        {
                            rout1[index1] = rout[l];
                            ++index1;
                        }


                        for (size_t l = j; l >= i + 1; --l)
                        {
                            rout2[index2] = rout[l];
                            ++index2;
                        }
                        for (size_t l = k; l >= j + 1; --l)
                        {
                            rout2[index2] = rout[l];
                            ++index2;
                        }
                        for (size_t l = k + 1; l < rout.size(); ++l)
                        {
                            rout2[index2] = rout[l];
                            ++index2;
                        }

                        len = dist_mat[0][rout1[0]] + dist_mat[rout1.back()][0];
                        for (size_t l = 1; l < rout.size(); ++l)
                            len += dist_mat[rout1[l - 1]][rout1[l]];
                        if (best_len > len)
                        {
                            best_len = len;
                            rout = rout1;
                        }

                        len = dist_mat[0][rout2[0]] + dist_mat[rout2.back()][0];
                        for (size_t l = 1; l < rout.size(); ++l)
                            len += dist_mat[rout2[l - 1]][rout2[l]];
                        if (best_len > len)
                        {
                            best_len = len;
                            rout = rout2;
                        }
                    }
        }

        void opt_3_fast2(std::vector<size_t>& rout, const matrix& dist_mat)
        {
            if (rout.size() < 4)
                return;
            auto best_len = utils::length_rout(rout, dist_mat);
            double len = 0;

            // 3 ребра смежны
            for (size_t i = 0; i < rout.size() - 4; ++i)
            {
                // swap(rout[i + 1], rout[i + 2]);

                len = best_len
                    - dist_mat[rout[i]][rout[i + 1]]
                    - dist_mat[rout[i + 1]][rout[i + 2]]
                    - dist_mat[rout[i + 2]][rout[i + 3]]
                    + dist_mat[rout[i]][rout[i + 2]]
                    + dist_mat[rout[i + 2]][rout[i + 1]]
                    + dist_mat[rout[i + 1]][rout[i + 3]];

                if (best_len > len)
                {
                    best_len = len;
                    std::swap(rout[i + 1], rout[i + 2]);
                }
            }

            // on the ends:
            {
                len = best_len
                    - dist_mat[0][rout[0]]
                    - dist_mat[rout[0]][rout[1]]
                    - dist_mat[rout[1]][rout[2]]
                    + dist_mat[0][rout[1]]
                    + dist_mat[rout[1]][rout[0]]
                    + dist_mat[rout[0]][rout[2]];

                if (best_len > len)
                {
                    best_len = len;
                    std::swap(rout[0], rout[1]);
                }

                size_t end_index = rout.size() - 3;
                len = best_len
                    - dist_mat[rout[end_index]][rout[end_index + 1]]
                    - dist_mat[rout[end_index + 1]][rout[end_index + 2]]
                    - dist_mat[rout[end_index + 2]][0]
                    + dist_mat[rout[end_index]][rout[end_index + 2]]
                    + dist_mat[rout[end_index + 2]][rout[end_index + 1]]
                    + dist_mat[rout[end_index + 1]][0];

                if (best_len > len)
                {
                    best_len = len;
                    std::swap(rout[end_index + 1], rout[end_index + 2]);
                }
            }

            // 2 ребра смежны, пусть вторые 2
            std::vector<size_t> rout1(rout.size());
            std::vector<size_t> rout2(rout.size());
            size_t index1 = 0;
            size_t index2 = 0;
            for (size_t i = 0; i < rout.size() - 4; ++i)
                for (size_t j = i + 2; j < rout.size() - 1; ++j)
                {
                    index1 = 0;
                    index2 = 0;

                    for (size_t k = 0; k <= i; ++k)
                    {
                        rout1[index1] = rout[k];
                        rout2[index2] = rout[k];
                        ++index1;
                        ++index2;
                    }

                    rout1[index1] = rout[j + 1];
                    ++index1;
                    for (size_t k = i + 1; k <= j; ++k)
                    {
                        rout1[index1] = rout[k];
                        ++index1;
                    }
                    for (size_t k = j + 2; k < rout.size(); ++k)
                    {
                        rout1[index1] = rout[k];
                        ++index1;
                    }

                    for (size_t k = j; k > i; --k)
                    {
                        rout2[index2] = rout[k];
                        ++index2;
                    }
                    for (size_t k = j + 1; k < rout.size(); ++k)
                    {
                        rout2[index2] = rout[k];
                        ++index2;
                    }

                    len = dist_mat[0][rout1[0]] + dist_mat[rout1.back()][0];
                    for (size_t k = 1; k < rout.size(); ++k)
                        len += dist_mat[rout1[k - 1]][rout1[k]];
                    if (best_len > len)
                    {
                        best_len = len;
                        rout = rout1;
                    }

                    len = dist_mat[0][rout2[0]] + dist_mat[rout2.back()][0];
                    for (size_t k = 1; k < rout.size(); ++k)
                        len += dist_mat[rout2[k - 1]][rout2[k]];
                    if (best_len > len)
                    {
                        best_len = len;
                        rout = rout2;
                    }
                }

            // 2 ребра смежны, пусть первые 2

            for (size_t j = 0; j < rout.size() - 1; ++j)
                for (size_t i = j + 3; i < rout.size() - 4; ++i)
                {
                    index1 = 0;
                    index2 = 0;

                    for (size_t k = 0; k <= j; ++k)
                    {
                        rout1[index1] = rout[k];
                        rout2[index2] = rout[k];
                        ++index1;
                        ++index2;
                    }

                    for (size_t k = j + 2; k <= i; ++k)
                    {
                        rout1[index1] = rout[k];
                        ++index1;
                    }
                    rout1[index1] = rout[j + 1];
                    ++index1;
                    for (size_t k = i + 1; k < rout.size(); ++k)
                    {
                        rout1[index1] = rout[k];
                        ++index1;
                    }

                    for (size_t k = i; k > j; --k)
                    {
                        rout2[index2] = rout[k];
                        ++index2;
                    }
                    for (size_t k = i + 1; k < rout.size(); ++k)
                    {
                        rout2[index2] = rout[k];
                        ++index2;
                    }

                    len = dist_mat[0][rout1[0]] + dist_mat[rout1.back()][0];
                    for (size_t k = 1; k < rout.size(); ++k)
                        len += dist_mat[rout1[k - 1]][rout1[k]];
                    if (best_len > len)
                    {
                        best_len = len;
                        rout = rout1;
                    }

                    len = dist_mat[0][rout2[0]] + dist_mat[rout2.back()][0];
                    for (size_t k = 1; k < rout.size(); ++k)
                        len += dist_mat[rout2[k - 1]][rout2[k]];
                    if (best_len > len)
                    {
                        best_len = len;
                        rout = rout2;
                    }
                }


            // 3 ребра не смежны

            double len_i = dist_mat[0][rout1[0]];
            double len1_k_j = 0;
            double len_j2 = 0;
            double len_k1 = 0;
            double len_k2 = 0;
            for (size_t i = 0; i < rout.size() - 5; ++i)
            {
                rout1[i] = rout[i];
                rout2[i] = rout[i];
                if (i > 0)
                    len_i += dist_mat[rout[i - 1]][rout[i]];
                for (size_t j = i + 2; j < rout.size() - 3; ++j)
                {
                    for (size_t k = j + 2; k < rout.size() - 1; ++k)
                    {
                        index1 = i + 1;
                        index2 = i + 1;
                        index1 = 0;
                        index2 = 0;


                        for (size_t l = 0; l <= i; ++l)
                        {
                            rout1[index1] = rout[l];
                            rout2[index2] = rout[l];
                            ++index1;
                            ++index2;
                        }


                        if (k - 2 > j)
                            len1_k_j += dist_mat[rout[k - 1]][k]
                            - dist_mat[rout[j]][k];
                        else
                            len1_k_j = 2;
                        for (size_t l = k; l >= j + 1; --l)
                        {
                            rout1[index1] = rout[l];
                            ++index1;
                        }
                        for (size_t l = i + 1; l <= j; ++l)
                        {
                            rout1[index1] = rout[l];
                            ++index1;
                        }
                        for (size_t l = k + 1; l < rout.size(); ++l)
                        {
                            rout1[index1] = rout[l];
                            ++index1;
                        }


                        for (size_t l = j; l >= i + 1; --l)
                        {
                            rout2[index2] = rout[l];
                            ++index2;
                        }
                        for (size_t l = k; l >= j + 1; --l)
                        {
                            rout2[index2] = rout[l];
                            ++index2;
                        }
                        for (size_t l = k + 1; l < rout.size(); ++l)
                        {
                            rout2[index2] = rout[l];
                            ++index2;
                        }

                        len = dist_mat[0][rout1[0]] + dist_mat[rout1.back()][0];
                        for (size_t l = 1; l < rout.size(); ++l)
                            len += dist_mat[rout1[l - 1]][rout1[l]];
                        if (best_len > len)
                        {
                            best_len = len;
                            rout = rout1;
                        }

                        len = dist_mat[0][rout2[0]] + dist_mat[rout2.back()][0];
                        for (size_t l = 1; l < rout.size(); ++l)
                            len += dist_mat[rout2[l - 1]][rout2[l]];
                        if (best_len > len)
                        {
                            best_len = len;
                            rout = rout2;
                        }
                    }
                }
            }
        }

        void swap(std::vector<size_t>& rout, const matrix& dist_mat)
        {
            if (rout.size() < 3)
                return;
            auto best_len = utils::length_rout(rout, dist_mat);
            auto best_rout = rout;

            for (int i = 0; i < (int)rout.size() - 1; ++i)
                for (int j = i + 1; j < (int)rout.size(); ++j)
                {
                    std::swap(rout[i], rout[j]);
                    auto len = utils::length_rout(rout, dist_mat);
                    if (best_len > len)
                    {
                        best_len = len;
                        best_rout = rout;
                    }
                    std::swap(rout[i], rout[j]);
                }
            rout = best_rout;
        }

        void shift(std::vector<size_t>& rout, const matrix& dist_mat)
        {
            if (rout.size() < 3)
                return;
            auto best_len = utils::length_rout(rout, dist_mat);
            auto best_rout = rout;

            for (int i = 0; i < (int)rout.size() - 1; ++i)
            {
                rout = best_rout;
                for (int j = i; j < (int)rout.size() - 1; ++j)
                {
                    std::swap(rout[j], rout[j + 1]);
                    auto len = utils::length_rout(rout, dist_mat);
                    if (best_len > len)
                    {
                        best_len = len;
                        best_rout = rout;
                    }
                }
            }
            rout = best_rout;
        }

        void exchange(std::vector<size_t>& rout, const matrix& dist_mat)
        {
            size_t shift = rout.size() / 5;
            if (rout.size() < 3)
                return;
            auto best_len = utils::length_rout(rout, dist_mat);
            auto best_rout = rout;

            for (int i = 0; i < (int)rout.size() - 1; ++i)
                for (int j = i + shift; j < (int)rout.size(); ++j)
                {
                    


                    auto len = utils::length_rout(rout, dist_mat);
                    if (best_len > len)
                    {
                        best_len = len;
                        best_rout = rout;
                    }
                    std::swap(rout[i], rout[j]);
                }
            rout = best_rout;
        }
    }

    namespace local_opt_for_VND_STS
    {
        void TSP_2_opt(vec_int_float& rout, const matrix& dist_mat)
        {
            if (rout.size() < 3)
                return;
            auto best_len = utils::length_rout(rout, dist_mat);
            auto best_rout = rout;
            for (int i = 0; i < (int)rout.size() - 2; ++i)
                for (int j = i + 1; j < (int)rout.size() - 1; ++j)
                {
                    rout = vec_int_float();
                    for (int k = 0; k < i; ++k)
                        rout.push_back(best_rout[k]);
                    for (int k = j; k >= i; --k)
                        rout.push_back(best_rout[k]);
                    for (size_t k = j + 1; k < best_rout.size(); ++k)
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

        void TSP_3_opt(vec_int_float& rout, const matrix& dist_mat)
        {
            if (rout.size() < 4)
                return;
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
                    vec_int_float rout1;
                    vec_int_float rout2;

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
                    vec_int_float rout1;
                    vec_int_float rout2;

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
                        vec_int_float rout1;
                        vec_int_float rout2;

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

        void swap(vec_int_float& rout, const matrix& dist_mat)
        {
            if (rout.size() < 3)
                return;
            auto best_len = utils::length_rout(rout, dist_mat);
            auto best_rout = rout;

            for (int i = 0; i < (int)rout.size() - 1; ++i)
                for (int j = i + 1; j < (int)rout.size(); ++j)
                {
                    std::swap(rout[i], rout[j]);
                    auto len = utils::length_rout(rout, dist_mat);
                    if (best_len > len)
                    {
                        best_len = len;
                        best_rout = rout;
                    }
                    std::swap(rout[i], rout[j]);
                }
            rout = best_rout;
        }

        void shift(vec_int_float& rout, const matrix& dist_mat)
        {
            if (rout.size() < 3)
                return;
            auto best_len = utils::length_rout(rout, dist_mat);
            auto best_rout = rout;

            for (int i = 0; i < (int)rout.size() - 1; ++i)
            {
                rout = best_rout;
                for (int j = i; j < (int)rout.size() - 1; ++j)
                {
                    std::swap(rout[j], rout[j + 1]);
                    auto len = utils::length_rout(rout, dist_mat);
                    if (best_len > len)
                    {
                        best_len = len;
                        best_rout = rout;
                    }
                }
            }
            rout = best_rout;
        }

        void exchange(vec_int_float& rout, const size_t a, const size_t b)
        {
            if (rout[a].second < rout[b].second)
                std::swap(rout[a], rout[b]);

            size_t ind = rout[b].first;
            rout[b].first = rout[a].first;
            rout[a].second -= rout[b].second;

            rout.emplace(rout.begin() + a, std::pair<size_t, double>( ind, rout[b].second ));
        }
    }
}