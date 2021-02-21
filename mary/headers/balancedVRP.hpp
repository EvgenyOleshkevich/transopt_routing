#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <random>
#include "../headers/utils.hpp"

namespace balancedVRP
{
	namespace clustering
	{
		std::pair<int, int> two_farthest_vertex(const matrix&, const std::vector<size_t>&);

		void update_dist_to_cluaster(const matrix& ,
			const std::vector<size_t>&, std::vector<double>&,
			const int, const int);

		int_matrix dichotomous_division(const matrix&, const int);

		int_matrix dichotomous_division(const matrix&, const std::vector<size_t>&, const int, const int);
	
		std::vector<size_t> radian_sort(const double* const, const double* const, const size_t);
		
		int_matrix sweeping(const double* const, const double* const, const size_t, size_t);
	}

	int_matrix cutting_rout(const std::vector<size_t>&, const matrix&, size_t);

}