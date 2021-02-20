#include "../headers/balancedVRP.hpp"
#include "../headers/TSP.hpp"
#include "../headers/utils.hpp"
#include <math.h>

namespace TSP
{
    void Lin_Kernighan_by_rout(std::vector<size_t>& rout, const matrix& dist_mat) {

        
        
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