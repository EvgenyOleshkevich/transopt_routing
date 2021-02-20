#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <random>


using matrix = std::vector<std::vector<double>>;
using int_matrix = std::vector<std::vector<size_t>>;

namespace TSP
{
    std::vector<size_t> Lin_Kernighan(const matrix&);

    namespace local_opt
    {
        void TSP_2_opt(std::vector<size_t>&, const matrix&);

        void TSP_3_opt(std::vector<size_t>&, const matrix&);
    }

}