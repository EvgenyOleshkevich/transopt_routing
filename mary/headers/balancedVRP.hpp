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
	
		double demand_of_vertex(const std::vector<double>&, const std::vector<size_t>&);

		int_matrix dichotomous_division_weight(const matrix&,
			const std::vector<double>&, const double, const size_t);

		int_matrix dichotomous_division_weight(const matrix&, const std::vector<double>&,
			const std::vector<size_t>&, const double, const size_t, const bool);

		std::vector<size_t> radian_sort(const double* const, const double* const, const size_t);
		
		int_matrix sweeping(const double* const, const double* const, const size_t, size_t);
	}

	int_matrix cutting_rout(const std::vector<size_t>&, const matrix&, size_t);

	class VND_STS
	{
	public:
		VND_STS(const matrix& dist_mat, const size_t count_point, const size_t need_routs
			, const double capacity, const std::vector<double>& weight);

		// can be faster without utils::length_rout(rout, dist_mat);
		std::pair<int_matrix, matrix> dynamic_decode(const vec_int_float&);
		// can be faster without utils::length_rout(rout, dist_mat);
		std::pair<int_matrix, matrix> greedy_decode(const vec_int_float&);

		std::pair<int_matrix, matrix> calculate(vec_int_float&);
	private:
		matrix dist_mat;
		const size_t count_point;
		const size_t need_routs;
		const double capacity;
		const std::vector<double> weight;
		const size_t fict_vertex;

		double VND(vec_int_float&);

		double STS(vec_int_float&);

		// can be faster without utils::length_rout(rout, dist_mat);
		double dynamic_decode_fast(const vec_int_float&);

		// can be faster without utils::length_rout(rout, dist_mat);
		void dynamic_decode_add_fict(vec_int_float&);

		// can be faster without utils::length_rout(rout, dist_mat);
		double greedy_decode_fast(const vec_int_float&);
	};
}