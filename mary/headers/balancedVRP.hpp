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
		pair<int, int> two_farthest_vertex(const matrix&, const vector<size_t>&);

		void update_dist_to_cluaster(const matrix& ,
			const vector<size_t>&, vector<double>&,
			const int, const int);

		int_matrix dichotomous_division(const matrix&, const int);

		int_matrix dichotomous_division(const matrix&, const vector<size_t>&, const int, const int);
	
		double demand_of_vertex(const vector<double>&, const vector<size_t>&);

		int_matrix dichotomous_division_weight(const matrix&,
			const vector<double>&, const double, const size_t);

		int_matrix dichotomous_division_weight(const matrix&, const vector<double>&,
			const vector<size_t>&, const double, const size_t, const bool);

		vector<size_t> radian_sort(const vector<double>&, const vector<double>&, const size_t);
		
		int_matrix sweeping(const double* const, const double* const, const size_t, size_t);

		// получение матрицы расстояние для каждого кластера
		vector<matrix> get_dist_inner_cluster(const matrix&, const int_matrix&);

		// получение весов для каждого кластера
		matrix get_weight_inner_cluster(const vector<double>&, const int_matrix&);

		// получение номер кластера для каждой вершины
		vector<int> get_number_cluster_by_vertex(const int_matrix&);
	}

	namespace project
	{
		int_matrix clusters_with_multiple_duplicates(const int_matrix&, const vector<size_t>&);
	}

	// разрезание общего маршрута
	int_matrix cutting_rout(const vector<size_t>&, const matrix&, size_t);

	class VND_STS
	{
	public:
		VND_STS(const matrix& dist_mat, const size_t need_routs
			, const double capacity, const vector<double>& weight);

		// can be faster without utils::length_rout(rout, dist_mat);
		pair<int_matrix, matrix> dynamic_decode(const vec_int_float&);
		// can be faster without utils::length_rout(rout, dist_mat);
		pair<int_matrix, matrix> greedy_decode(const vec_int_float&);

		pair<int_matrix, matrix> calculate(vec_int_float&);
	private:
		matrix dist_mat;
		const size_t count_point;
		const size_t need_routs;
		const double capacity;
		const vector<double> weight;
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

	struct Transport
	{
		double volume;
		double cost_start;
		double cost_by_dist;
		size_t count;
	};

	class Ant
	{
	public:
		Ant(const matrix& dist_mat, const vector<Transport>& transports, const vector<double>& weights);

		pair<int_matrix, matrix> calculate(vec_int_float&);
	private:
		const matrix& dist_mat;
		const vector<Transport>& transports;
		const vector<double>& weights;
		vector<vector<pair<double, size_t>>> prioritet;
		vector <matrix> pheromone_mat;
		const double base_pheromone = 0.2;
		const double min_pheromone = 0.1;
		const double max_pheromone = 0.85;
		const double alpha = 1;
		const double beta = 4;

		const double evaporation_rate = 0.8; // коэфициент испарения
		const double addition_pheromone_rate = 0.8; // коэфициент добавленяи феромона
		const double Q_rate = 5; // коэфициент добавленяи феромона
		std::mt19937 mersenne_rand;

		struct Calc_data
		{
			vector<double> weights_used;
			vector<size_t> transport_used;
			vector<size_t> transport_current_pos;
			vector<size_t> vertex_used;

		};

		void fill_pheromone() {
			pheromone_mat = vector <matrix>(
				transports.size(),
				matrix(
					dist_mat.size(),
					vector<double>(dist_mat.size(), 0)
				));

			for (size_t k = 0; k < pheromone_mat.size(); ++k) {
				for (size_t i = 0; i < pheromone_mat[k].size(); ++i)
					for (size_t j = 1; j < pheromone_mat[k].size(); ++j)
						pheromone_mat[k][i][j] = base_pheromone;

				//for (size_t i = 0; i < pheromone_mat[k].size(); ++i)
				//{
				//	pheromone_mat[k][i][i] = 0;
				//}
			}
		}

		void calcalate_pheromone(const int_matrix& routs, matrix& pher_mat, const Transport& transport) {
			for (size_t i = 0; i < pher_mat.size(); ++i)
				for (size_t j = 1; j < pher_mat.size(); ++j)
					pher_mat[i][j] *= evaporation_rate;

			for (const vector<size_t>& rout : routs) {
				double addition = Q_rate * (rout.size() - 1)
					/ (utils::length_rout_0(rout, dist_mat) * transport.cost_by_dist + transport.cost_start);

				for (size_t i = 1; i < rout.size(); ++i)
					pher_mat[rout[i - 1]][rout[i]] += addition;
			}

			for (size_t i = 0; i < pher_mat.size(); ++i)
				for (size_t j = 1; j < pher_mat.size(); ++j)
					pher_mat[i][j] = std::max(min_pheromone, std::min(evaporation_rate, max_pheromone));
		}

		pair<size_t, size_t> get_next_step(const Calc_data& data)
		{
			double sum = 0;

			auto cost_move = vector<vector<double>>(transports.size(), vector<double>(data.vertex_used.size(), 0));
			for (size_t trans = 0; trans < transports.size(); ++trans)
			{
				size_t pos = data.transport_current_pos[trans];
				for (size_t i = 0; i < data.vertex_used.size(); ++i)
					if (data.vertex_used[i] == 0)
					{
						cost_move[trans][i] = pow(pheromone_mat[trans][pos][i], alpha)
							* pow(1 / dist_mat[pos][i] * transports[trans].cost_by_dist, beta);
						sum += cost_move[trans][i];
					}
			}

			auto rand_value = (double)mersenne_rand() * sum / UINT32_MAX;
			double sum_after = 0;

			for (size_t trans = 0; trans < transports.size(); ++trans)
				for (size_t i = 0; i < data.vertex_used.size(); ++i)
				{
					sum_after += cost_move[trans][i];
					if (sum_after > rand_value)
						return { trans , i };
				}
			
			throw new std::runtime_error("rand_value is bigger sum weight");
		}
	};
}