#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <random>


using matrix = std::vector<std::vector<double>>;
using int_matrix = std::vector<std::vector<size_t>>;

namespace balancedVRP
{
	namespace clustering
	{
		std::pair<int, int> two_farthest_vertex(const matrix&, const std::vector<size_t>&);

		void update_dist_to_cluaster(const matrix& ,
			const std::vector<size_t>&, std::vector<double>&,
			const int, const int);

		int_matrix dichotomous_division(const matrix&, const int);

		int_matrix dichotomous_division(const matrix&, const std::vector<size_t>&, const int);
	}

}