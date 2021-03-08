#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <random>
#include "utils.hpp"

namespace TSP
{
    namespace Lin_Kernighan
    {
        std::vector<size_t> Lin_Kernighan(const double* const, const double* const, const size_t);

        void Lin_Kernighan_by_rout(std::vector<size_t>&, const matrix&, const sorted_matrix&);
    }

    namespace local_opt
    {
        void TSP_2_opt(std::vector<size_t>&, const matrix&);

        void TSP_3_opt(std::vector<size_t>&, const matrix&);

        void swap(std::vector<size_t>&, const matrix&);

        void shift(std::vector<size_t>&, const matrix&);

        void exchange(std::vector<size_t>&, const matrix&);
    }

    namespace local_opt_for_VND_STS
    {
        void TSP_2_opt(vec_int_float&, const matrix&);

        void TSP_3_opt(vec_int_float&, const matrix&);

        void swap(vec_int_float&, const matrix&);

        void shift(vec_int_float&, const matrix&);

        void exchange(vec_int_float&, const size_t, const size_t);
    }

}