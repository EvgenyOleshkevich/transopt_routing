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

		// получение матрицы рассто€ние дл€ каждого кластера
		vector<matrix> get_dist_inner_cluster(const matrix&, const int_matrix&);

		// получение весов дл€ каждого кластера
		matrix get_weight_inner_cluster(const vector<double>&, const int_matrix&);

		// получение номер кластера дл€ каждой вершины
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
		double capacity;
		double cost_start;
		double cost_by_dist;
		size_t count;
	};

	class Ant_algorithm
	{
	public:
		Ant_algorithm(const matrix& dist_mat, const vector<Transport>& transports, const vector<double>& weights);

		struct Ant
		{
			vector<size_t> rout;
			double remain_volume;
			size_t transport_type;
			size_t current_pos;
			bool is_last; // дл€ определени€ нужно ли после его заполнени€ пускать нового
		};

		struct Calc_data
		{
			//vector<double> weights_used;
			//vector<size_t> transport_used;
			//vector<size_t> transport_current_pos;
			vector<Ant> ants;
			vector<size_t> transport_remain;
			vector<size_t> vertex_used;
		};

		pair<vector<int_matrix>, double> calculate()
		{
			double best_lenght = UINT32_MAX;
			vector<int_matrix> best_res;


			double min_weight = min_weight_vertex();
			double max_weight = max_weight_vertex();

			unsigned int start_time = clock();
			unsigned int cur_time = start_time;
			unsigned int prev_cur_time = start_time;
			while (true)
			{
				fill_pheromone();

				vector<int_matrix> routs(transports.size());
				for (size_t iter = 1; iter < count_iter_inner; ++iter) {

					Calc_data data = fill_data();
					vector<int_matrix> res(transports.size());

					for (size_t used = 1; used < weights.size(); ++used)
					{
						if (used == weights.size() - 1)
						{
							int y = 0;
						}
						auto trans_vertex = get_next_step(data);

						if (trans_vertex.second == 0)
						{
							get_next_step_sharing(data);
						}

						size_t trans = trans_vertex.first;
						size_t vertex = trans_vertex.second;
						data.ants[trans].current_pos = vertex;
						data.ants[trans].remain_volume -= weights[vertex];
						data.ants[trans].rout.push_back(vertex);
						data.vertex_used[vertex] = 1;

						size_t transport_type = data.ants[trans].transport_type;

						if (data.ants[trans].is_last
							&& data.transport_remain[transport_type] > 0
							&& data.ants[trans].remain_volume <= max_weight)
						{
							data.ants[trans].is_last = false;
							--data.transport_remain[transport_type];
							data.ants.push_back({
								vector<size_t>(1, 0),
								transports[transport_type].capacity,
								transport_type,
								0 ,
								true });
						}

						if (data.ants[trans].remain_volume < min_weight) {
							data.ants[trans].rout.push_back(0);
							routs[transport_type].push_back(data.ants[trans].rout);
							res[transport_type].push_back(data.ants[trans].rout);
							data.ants.erase(data.ants.begin() + trans);
						}
					}
					for (Ant& ant: data.ants)
					{
						ant.rout.push_back(0);
						routs[ant.transport_type].push_back(ant.rout);
						res[ant.transport_type].push_back(ant.rout);
					}

					double lenght = length_roust(res);

					if (best_lenght > lenght) {
						best_lenght = lenght;
						best_res = res;
					}
				}
				// 130-131
				calcalate_pheromone(routs);
				unsigned int cur_time = clock();
				unsigned int diff = cur_time - prev_cur_time;
				prev_cur_time = cur_time;
				std::cout << diff << std::endl;
				if (cur_time - start_time > 60000)
					break;
			}

			return { best_res,best_lenght };
		}

	private:
		const matrix& dist_mat;
		const vector<Transport>& transports;
		const vector<double>& weights;
		vector <matrix> pheromone_mat;
		const double base_pheromone = 0.2;
		const double min_pheromone = 0.1;
		const double max_pheromone = 0.85;
		//const double first_step_decrease = 0.05;
		const double alpha = 1;
		const double beta = 4;
		const size_t count_iter_inner = 50;

		const double evaporation_rate = 0.8; // коэфициент испарени€
		//const double addition_pheromone_rate = 0.8; // коэфициент добавлен€и феромона
		const double Q_rate = 5; // коэфициент добавлен€и феромона
		//const double criteria_start_ant = 0.3;
		//const double criteria_close_ant = 0.99;
		std::mt19937 mersenne_rand;

		double min_weight_vertex() const
		{
			double min = weights[1];
			for (size_t i = 2; i < weights.size(); ++i)
				if (min > weights[i])
					min = weights[i];

			return min;
		}

		double max_weight_vertex() const
		{
			double max = weights[1];
			for (size_t i = 2; i < weights.size(); ++i)
				if (max < weights[i])
					max = weights[i];

			return max;
		}

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

				/*for (size_t i = 0; i < pheromone_mat[k].size(); ++i)
					pheromone_mat[k][0][i] -= first_step_decrease;*/
			}
		}

		Calc_data fill_data() {
			Calc_data data{ 
				vector<Ant>(),
				vector<size_t>(transports.size()),
				vector<size_t>(dist_mat.size(), 0)
			};

			for (size_t i = 0; i < transports.size(); ++i)
			{
				data.transport_remain[i] = transports[i].count - 1;
				data.ants.push_back({ vector<size_t>(1, 0), transports[i].capacity, i, 0 , true });
			}

			data.vertex_used[0] = 1;
			return data;
		}

		void calcalate_pheromone(const vector<int_matrix>& routs) {

			for (size_t trans_type = 0; trans_type < routs.size(); ++trans_type)
			{
				for (size_t i = 0; i < pheromone_mat[trans_type].size(); ++i)
					for (size_t j = 1; j < pheromone_mat[trans_type].size(); ++j)
						pheromone_mat[trans_type][i][j] *= evaporation_rate;

				for (const vector<size_t>& rout : routs[trans_type]) {
					double addition = Q_rate * (rout.size() - 1)
						/ (utils::length_rout_0(rout, dist_mat) * transports[trans_type].cost_by_dist + transports[trans_type].cost_start);

					for (size_t i = 1; i < rout.size(); ++i)
						pheromone_mat[trans_type][rout[i - 1]][rout[i]] += addition;
				}

				/*for (size_t i = 0; i < pher_mat.size(); ++i)
					pher_mat[0][i] -= first_step_decrease;*/

				for (size_t i = 0; i < pheromone_mat[trans_type].size(); ++i)
					for (size_t j = 1; j < pheromone_mat[trans_type].size(); ++j)
						pheromone_mat[trans_type][i][j] = std::max(min_pheromone, std::min(evaporation_rate, max_pheromone));
			}
		}

		pair<size_t, size_t> get_next_step(const Calc_data& data)
		{
			double sum = 0;

			auto cost_move = vector<vector<double>>(data.ants.size(), vector<double>(data.vertex_used.size(), 0));
			for (size_t trans = 0; trans < data.ants.size(); ++trans)
			{
				const Ant& ant = data.ants[trans];
				size_t pos = ant.current_pos;
				for (size_t i = 1; i < data.vertex_used.size(); ++i)
					if (data.vertex_used[i] == 0
						&& ant.remain_volume > weights[i]) {

						if (dist_mat[pos][i] > EPS)
							cost_move[trans][i] = pow(pheromone_mat[ant.transport_type][pos][i], alpha)
							* pow(1 / (dist_mat[pos][i] * transports[ant.transport_type].cost_by_dist), beta);
						else
						{
							return { trans , i };

						}
						sum += cost_move[trans][i];
					}
			}

			auto rand_value = (double)mersenne_rand() * sum / UINT32_MAX;
			double sum_after = 0;

			for (size_t trans = 0; trans < data.ants.size(); ++trans)
				for (size_t i = 1; i < data.vertex_used.size(); ++i)
				{
					sum_after += cost_move[trans][i];
					if (sum_after > rand_value)
						return { trans , i };
				}

			return { 0 , 0 };
		}

		pair<vector<size_t>, size_t> get_next_step_sharing(const Calc_data& data)
		{
			vector<size_t> used_trans;
			size_t vertex = 0;
			for (size_t i = 1; i < data.vertex_used.size(); ++i)
				if (data.vertex_used[i] == 0)
				{
					vertex = i;
					break;
				}
			double weight = weights[vertex];
			vector<pair<double, size_t>> unused(data.ants.size());

			for (size_t trans = 0; trans < data.ants.size(); ++trans)
			{
				const Ant& ant = data.ants[trans];
				unused[trans].first = ant.remain_volume * pow(pheromone_mat[ant.transport_type][ant.current_pos][vertex], alpha)
					* pow(1 / (dist_mat[ant.current_pos][vertex] * transports[ant.transport_type].cost_by_dist), beta);
				unused[trans].second = trans;
			}
			std::sort(unused.begin(), unused.end());

			for (int i = unused.size() - 1; i > 0; --i)
			{
				weight -= data.ants[unused[i].second].remain_volume;
				used_trans.push_back(unused[i].second);
				if (weight <= 0)
					break;
			}
			return { used_trans , vertex };
		}

		double length_roust(vector<int_matrix> routs)
		{
			double lenght = 0;
			for (size_t trans_type = 0; trans_type < routs.size(); ++trans_type)
				for (const vector<size_t>& rout : routs[trans_type])
					lenght += utils::length_rout_0(rout, dist_mat) * transports[trans_type].cost_by_dist + transports[trans_type].cost_start;

			return lenght;
		}
	};
}