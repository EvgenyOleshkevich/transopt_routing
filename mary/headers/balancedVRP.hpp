#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <random>
#include "../headers/TSP.hpp"
//#include "../headers/utils.hpp"

namespace balancedVRP
{
	struct Transport
	{
		double capacity;
		double cost_start;
		double cost_by_dist;
		size_t count;
	};

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

		vector<size_t> radian_sort(const vector<double>&, const vector<double>&);
		
		int_matrix sweeping(const double* const, const double* const, const size_t, size_t);

		// получение матрицы расстояние для каждого кластера
		vector<matrix> get_dist_inner_cluster(const matrix&, const int_matrix&);

		// получение весов для каждого кластера
		matrix get_weight_inner_cluster(const vector<double>&, const int_matrix&);

		// получение координат для каждого кластера
		vector<int> get_coord_inner_cluster(const vector<double>&, const int_matrix&);

		// получение номер кластера для каждой вершины
		vector<int> get_number_cluster_by_vertex(const int_matrix&);
	}

	namespace project
	{
		int_matrix clusters_with_multiple_duplicates(const int_matrix&, const vector<size_t>&);

		vector<size_t> radian_sort(const vector<double>&, const vector<double>&);

		class Sweeping
		{
		public:
			Sweeping(const matrix& dist_mat, const vector<double>& x, const vector<double>& y,
				const vector<double>& weights, const vector<Transport>& transports)
				: dist_mat(dist_mat), transports(transports), x(x), y(y), weights(weights) {}

			void run()
			{
				vector<size_t> order = radian_sort(x, y);
				size_t order_i = 0;

				for (const Transport& transport : transports)
				{
					int_matrix routs;

					for (size_t i = 0; i < transport.count; ++i)
					{
						vector<size_t> rout(1, 0);
						double remain_weight = transport.capacity;
						while (remain_weight > weights[order[order_i]])
						{
							remain_weight -= weights[order[order_i]];
							rout.push_back(order[order_i]);
							++order_i;
							if (order_i == order.size())
								break;
						} 
						rout.push_back(0);
						routs.push_back(rout);
						if (order_i == order.size())
							break;
					}
					res.push_back(routs);
					if (order_i == order.size())
						break;
				}
			}

			vector<int_matrix> res;
			double length()
			{
					double length = 0;
					for (size_t trans_type = 0; trans_type < res.size(); ++trans_type)
						for (const vector<size_t>& rout : res[trans_type])
							length += utils::length_rout_0(rout, dist_mat) * transports[trans_type].cost_by_dist + transports[trans_type].cost_start;
					return length;
			}
		private:
			const matrix& dist_mat;
			const vector<Transport>& transports;
			const vector<double>& x;
			const vector<double>& y;
			const vector<double>& weights;

		};

		class GreadyBase
		{
		public:
			GreadyBase(const matrix& dist_mat, const sorted_matrix& sorted_dist_mat,
				const vector<double>& weights, const vector<Transport>& transports)
				: sorted_dist_mat(sorted_dist_mat), dist_mat(dist_mat),
				transports(transports), weights(weights) {}

			void run()
			{
				vector <size_t> used(dist_mat.size(), 0);
				used[0] = 1;
				size_t count_add = 1;

				for (const Transport& transport : transports)
				{
					int_matrix routs;

					for (size_t i = 0; i < transport.count; ++i)
					{
						vector<size_t> rout(1, 0);
						double remain_weight = transport.capacity;
						size_t vertex = 0;
						for (size_t order_i = 0; order_i < sorted_dist_mat.size() && used[vertex] == 1; ++order_i)
							vertex = sorted_dist_mat[0][order_i].second;

						while (true/*remain_weight >= weights[vertex]*/)
						{
							used[vertex] = 1;
							remain_weight -= weights[vertex];
							rout.push_back(vertex);
							++count_add;
							if (count_add == used.size())
								break;

							size_t new_vertex = vertex;
							size_t order_i = 0;
							for (;
								order_i < sorted_dist_mat.size()
								&& (used[new_vertex] == 1
								|| remain_weight < weights[new_vertex]);
								++order_i)
							{
								new_vertex = sorted_dist_mat[vertex][order_i].second;
							}
							if (used[new_vertex] == 1 || remain_weight < weights[new_vertex])
								break;
							vertex = new_vertex;
						}
						//std::cout << "remain_weight: " << remain_weight << std::endl;
						rout.push_back(0);
						routs.push_back(rout);
						if (count_add == used.size())
							break;
					}
					res.push_back(routs);
					if (count_add == used.size())
						break;
				}
			}

			

			vector<int_matrix> res;
			double length()
			{
				double length = 0;
				for (size_t trans_type = 0; trans_type < res.size(); ++trans_type)
					for (const vector<size_t>& rout : res[trans_type])
						length += utils::length_rout_0(rout, dist_mat) * transports[trans_type].cost_by_dist + transports[trans_type].cost_start;
				return length;
			}
		private:
			const sorted_matrix& sorted_dist_mat;
			const matrix& dist_mat;
			const vector<Transport>& transports;
			const vector<double>& weights;

			void run2()
			{
				vector <size_t> order_i(dist_mat.size(), 0);
				vector <size_t> used(dist_mat.size(), 0);
				used[0] = 1;
				size_t count_add = 1;

				for (const Transport& transport : transports)
				{
					int_matrix routs;

					for (size_t i = 0; i < transport.count; ++i)
					{
						vector<size_t> rout(1, 0);
						double remain_weight = transport.capacity;
						size_t vertex = 0;
						while (used[vertex] == 1)
						{
							vertex = sorted_dist_mat[0][order_i[0]].second;
							++order_i[0];
						}

						while (remain_weight > weights[vertex])
						{
							used[vertex] = 1;
							remain_weight -= weights[vertex];
							rout.push_back(vertex);
							++count_add;
							if (count_add == used.size())
								break;

							size_t new_vertex = vertex;
							while (used[new_vertex] == 1)
							{
								new_vertex = sorted_dist_mat[vertex][order_i[vertex]].second;
								++order_i[vertex];
							}
							vertex = new_vertex;
						}
						rout.push_back(0);
						routs.push_back(rout);
						if (count_add == used.size())
							break;
					}
					res.push_back(routs);
					if (count_add == used.size())
						break;
				}
			}
		};

		class ClarkRight
		{
		public:
			ClarkRight(const matrix& dist_mat, 
				const vector<double>& weights, const vector<Transport>& transports)
				: dist_mat(dist_mat),
				transports(transports), weights(weights) {}

			void run()
			{
				{
					auto sotred_save_matrix = sort_save_matrix();
					auto routs_start = init_routs();
					auto routs_end = init_routs();
					size_t count_routs = dist_mat.size() - 1;
					while (count_routs > 1)
					{
						auto save_rout = sotred_save_matrix.back();
						sotred_save_matrix.pop_back();
						size_t i = (size_t)save_rout[1];
						size_t j = (size_t)save_rout[2];
						if (routs_start[j].empty() || routs_end[i].empty() || routs_end[i][0] == j)
							continue; // проверка, что они из разных маршрутов


						size_t m = routs_end[i][0];
						size_t n = routs_start[j].back();
						vector<size_t> vec;

						for (auto vertex : routs_end[i])
							vec.push_back(vertex);
						for (auto vertex : routs_start[j])
						{
							routs_start[m].push_back(vertex);
							vec.push_back(vertex);
						}
						routs_end[n] = vec;
						routs_end[i].clear();
						routs_start[j].clear();
						--count_routs;
					}


					vector<size_t> res_rout;
					for (const vector<size_t>& rout : routs_start)
					{
						if (rout.empty())
							continue;
						cut_rout(rout);
						break;
					}
				}
			}



			vector<int_matrix> res;
			double length()
			{
				double length = 0;
				for (size_t trans_type = 0; trans_type < res.size(); ++trans_type)
					for (const vector<size_t>& rout : res[trans_type])
						length += utils::length_rout_0(rout, dist_mat) * transports[trans_type].cost_by_dist + transports[trans_type].cost_start;
				return length;
			}
		private:
			const matrix& dist_mat;
			const vector<Transport>& transports;
			const vector<double>& weights;


			void cut_rout(const vector<size_t>& rout)
			{
				double cut_coef = get_cut_coef(rout);
				size_t rout_i = 0;
				for (size_t i = 0; i < transports.size(); ++i)
				{
					int_matrix routs;
					for (size_t j = 0; j < transports[i].count; ++j)
					{
						double remain_weight = transports[i].capacity;
						vector<size_t> rout_trans(1, 0);
						while (remain_weight - weights[rout[rout_i]] >= 0)
						{
							rout_trans.push_back(rout[rout_i]);
							remain_weight -= weights[rout[rout_i]];
							++rout_i;
							if (rout_i >= rout.size())
								break;
							if (dist_mat[rout[rout_i - 1]][rout[rout_i]] > cut_coef)
								break;

						}
						rout_trans.push_back(0);
						routs.push_back(rout_trans);
						if (rout_i >= rout.size())
							break;
					}
					
					res.push_back(routs);
					if (rout_i >= rout.size())
						break;
				}
			}

			double get_cut_coef(const vector<size_t>& rout)
			{
				const double speed_search = 1.1;

				double mean_edge = 0;
				size_t count_edge = 0;
				for (size_t i = 1; i < dist_mat.size(); ++i)
					for (size_t j = i + 1; j < dist_mat.size(); ++j)
					{
						++count_edge;
						mean_edge += dist_mat[i][j];
					}
				mean_edge /= count_edge;

				double cut_coef = 4;
				
				bool is_ok = false;
				while (!is_ok)
				{
					size_t rout_i = 0;
					for (size_t i = 0; i < transports.size(); ++i)
					{
						for (size_t j = 0; j < transports[i].count; ++j)
						{
							double remain_weight = transports[i].capacity;
							while (remain_weight - weights[rout[rout_i]] >= 0)
							{
								++rout_i;
								if (rout_i >= rout.size())
									break;
								if (dist_mat[rout[rout_i - 1]][rout[rout_i]] > cut_coef)
									break;
							}
							if (rout_i >= rout.size())
								break;
						}
						if (rout_i >= rout.size())
							break;
					}
					size_t c = 0;
					for (size_t i = 1; i < rout.size(); ++i)
						if (dist_mat[rout[i - 1]][rout[i]]) {
							++c;
						}
					if (rout_i < rout.size())
						cut_coef *= speed_search;
					else
						break;
				}

				return cut_coef;
			}

			matrix get_save_matrix()
			{
				matrix save_matrix(dist_mat.size(), vector<double>(dist_mat.size()));
				for (size_t i = 0; i < dist_mat.size(); ++i)
					for (size_t j = 0; j < dist_mat.size(); ++j)
						save_matrix[i][j] = dist_mat[i][0] + dist_mat[0][j] - dist_mat[i][j];
				return save_matrix;
			}

			matrix sort_save_matrix()
			{
				matrix save_matrix = get_save_matrix();
				matrix sorted_save_matrix;
				for (size_t i = 1; i < dist_mat.size(); ++i)
					for (size_t j = 1; j < dist_mat.size(); ++j)
					{
						if (i == j)
							continue;
						auto save = vector<double>(3);
						save[0] = save_matrix[i][j];
						save[1] = i;
						save[2] = j;
						sorted_save_matrix.push_back(save);
					}
				sort(sorted_save_matrix.begin(), sorted_save_matrix.end());
				//for (size_t i = 1; i < sorted_save_matrix.size() / 2; ++i)
				//	swap(sorted_save_matrix[i], sorted_save_matrix[sorted_save_matrix.size() - i]);
				//for (auto save :sorted_save_matrix)
				//    cout << save[0] << "; i=" << save[1] << "; j=" << save[2] << endl;
				return sorted_save_matrix;
			}

			int_matrix init_routs()
			{
				int_matrix routs(dist_mat.size());
				for (size_t i = 1; i < dist_mat.size(); ++i)
					routs[i].push_back(i);
				return routs;
			}
		};

		class Osman_new
		{
		public:
			Osman_new(const matrix& dist_mat, const vector<double>& weights,
				const vector<Transport>& transports, int_matrix routs,
				const vector<size_t> transport_id, const size_t count_rout) :
				res(routs), transport_id(transport_id),
				dist_mat(dist_mat), transports(transports), weights(weights),
				ts((size_t)(std::max(7.0, 9.6 * log(dist_mat.size() * count_rout) - 40))),
				count_rout(count_rout)// 13
			{}

			void run()
			{
				init_remain_weight();

				//std::unique_ptr<checker_change> checker(new checker_change_widht_local(this));
				auto routs = res;
				auto best_length = length();
				move_table = get_move_table();

				size_t max_iter = (size_t)(340 + 0.000353 * 5 * utils::sqr((double)dist_mat.size() * count_rout));
				fill_BSTM_RECM(routs);
				max_iter = 10;
				unsigned int start_time = clock();
				for (size_t i = 0; clock() - start_time < time_limit; ++i)
				{
					auto index = get_max_free_BSTM(i, routs);
					size_t p = index.first;
					size_t q = index.second;
					if (BSTM[p][q] - EPS < 0)
						break;

					auto a = RECM[p][q].first.first;
					auto b = RECM[p][q].first.second;
					auto way = RECM[p][q].second;
					switch (way)
					{
					case 0:
					{
						move_table[p][routs[p][a]] = i;
						move_table[q][routs[q][b]] = i;
						remain_weight[p] += weights[routs[p][a]] - weights[routs[q][b]];
						remain_weight[q] += weights[routs[q][b]] - weights[routs[p][a]];
						swap(routs[p][a], routs[q][b]);
						break;
					}
					case 1:
					{
						move_table[p][routs[p][a]] = i;
						remain_weight[p] += weights[routs[p][a]];
						remain_weight[q] -= weights[routs[p][a]];
						routs[q].emplace(routs[q].begin() + b, routs[p][a]);
						routs[p].erase(routs[p].begin() + a);
						break;
					}
					case 2:
					{
						move_table[q][routs[q][b]] = i;
						remain_weight[p] -= weights[routs[q][b]];
						remain_weight[q] += weights[routs[q][b]];
						routs[p].emplace(routs[p].begin() + a, routs[q][b]);
						routs[q].erase(routs[q].begin() + b);
						break;
					}
					case 3:
					{
						move_table[p][routs[p][a]] = i;
						remain_weight[p] += weights[routs[p][a]];
						remain_weight[q] -= weights[routs[p][a]];
						routs[q].push_back(routs[p][a]);
						routs[p].erase(routs[p].begin() + a);
						break;
					}
					case 4:
					{
						remain_weight[p] -= weights[routs[q][b]];
						remain_weight[q] += weights[routs[q][b]];
						move_table[q][routs[q][b]] = i;
						routs[p].push_back(routs[q][b]);
						routs[q].erase(routs[q].begin() + b);
						break;
					}
					default:
						break;
					}
					auto len = length(routs);
					if (best_length > len)
					{
						best_length = len;
						res = routs;
					}
					recalculate_BSTM_RECM(routs, p, q);
				}
			}

			double length()
			{
				double length = 0;
				for (size_t i = 0; i < res.size(); ++i)
				{
					double len = utils::length_rout(res[i], dist_mat) * transports[transport_id[i]].cost_by_dist;
					if (len >= EPS)
						length += len + transports[transport_id[i]].cost_start;
				}
				return length;
			}

			double length(const int_matrix& routs)
			{
				double length = 0;
				for (size_t i = 0; i < routs.size(); ++i)
				{
					double len = utils::length_rout(routs[i], dist_mat) * transports[transport_id[i]].cost_by_dist;
					if (len >= EPS)
						length += len + transports[transport_id[i]].cost_start;
				}
				return length;
			}

			static pair<int_matrix, vector<size_t>> transform_data(const vector<int_matrix>& routs)
			{
				int_matrix res;
				vector<size_t> transport_id;

				for (size_t trans_id = 0; trans_id < routs.size(); trans_id++)
					for (size_t i = 0; i < routs[trans_id].size(); i++)
					{
						transport_id.push_back(trans_id);
						res.push_back(routs[trans_id][i]);
					}
				//remove_zero_vertex(res);

				return { res , transport_id };
			}

			static void remove_zero_vertex(int_matrix& routs)
			{
				for (size_t i = 0; i < routs.size(); i++)
				{
					routs[i].pop_back();
					routs[i].erase(routs[i].begin());
				}
			}

			static void add_zero_vertex(int_matrix& routs)
			{
				for (size_t i = 0; i < routs.size(); i++)
				{
					routs[i].push_back(0);
					routs[i].emplace(routs[i].begin(), 0);
				}
			}

			int_matrix res;
			const vector<size_t> transport_id;
		private:
			const matrix& dist_mat;
			const vector<Transport>& transports;
			const vector<double>& weights;
			const size_t ts;// длина списка запрета
			const size_t count_rout;
			const double coef_access = 1.1;
			const double time_limit = 40000;

			matrix BSTM;
			vector<vector<pair<pair<size_t, size_t>, size_t>>> RECM;
			vector<double> remain_weight;
			vector<vector<int>> move_table; // rout*vertex
			// вершины и способ обмена
			// 0 - обмен, 1 - левая в правую перед указанной, 2 - правая в левую перед указанной
			// 3 - левая в правую в конец, 4 - правая в левую в конец, 5 - ничего

			void init_remain_weight()
			{
				remain_weight = vector<double>(res.size());
				for (size_t i = 0; i < res.size(); i++)
				{
					double weight = 0;
					for (const size_t vertex : res[i])
						weight += weights[vertex];

					remain_weight[i] = transports[transport_id[i]].capacity - weight;
				}
			}

			vector<vector<int>> get_move_table()
			{
				return vector<vector<int>>(count_rout, vector<int>(dist_mat.size(), -1000000000));
			}

			void check_change(const int_matrix& routs, size_t p, size_t q)
			{
				if (p > q)
					swap(p, q);

				const vector<size_t>& rout1 = routs[p];
				const vector<size_t>& rout2 = routs[q];
				pair<pair<size_t, size_t>, size_t> decision = { {0,0}, 5 };


				double cost1 = transports[transport_id[p]].cost_by_dist;
				double cost2 = transports[transport_id[q]].cost_by_dist;
				double cur_len = utils::length_rout(rout1, dist_mat)
					* cost1
					+ utils::length_rout(rout2, dist_mat)
					* cost2;
				double best_len = cur_len * coef_access;

				// обмен
				for (size_t i = 1; i < rout1.size() - 1; i++)
				{
					size_t A = rout1[i];
					double len = cur_len
						- (dist_mat[rout1[i - 1]][A]
							+ dist_mat[A][rout1[i + 1]]) * cost1;

					for (size_t j = 1; j < rout2.size() - 1; j++)
					{
						size_t B = rout2[j];
						if (remain_weight[p] + weights[A] - weights[B] < 0
							|| remain_weight[q] - weights[A] + weights[B] < 0)
							continue;

						double add_len = (dist_mat[rout2[j - 1]][A]
							+ dist_mat[A][rout2[j + 1]]) * cost2
							- (dist_mat[rout2[j - 1]][B]
								+ dist_mat[B][rout2[j + 1]]) * cost2
							+ (dist_mat[rout1[i - 1]][B]
								+ dist_mat[B][rout1[i + 1]]) * cost1;

						if (best_len > len + add_len)
						{
							best_len = len + add_len;
							decision = { {i,j}, 0 };
						}
					}
				}


				// 1 идет в 2
				if (rout1.size() > 2)
					for (size_t i = 1; i < rout1.size() - 1; i++)
					{
						size_t A = rout1[i];
						if (remain_weight[q] < weights[A]);
							continue;
						double len = cur_len
							- (dist_mat[rout1[i - 1]][A]
								+ dist_mat[A][rout1[i + 1]]) * cost1;

						for (size_t j = 1; j < rout2.size(); j++)
						{
							double add_len = (dist_mat[rout2[j - 1]][A]
								+ dist_mat[A][rout2[j]]) * cost2;

							if (best_len > len + add_len)
							{
								best_len = len + add_len;
								decision = { {i,j}, 1 };
							}
						}
					}

				// 2 идет в 1
				if (rout2.size() > 2)
					for (size_t i = 1; i < rout2.size() - 1; i++)
					{
						size_t A = rout2[i];
						if (remain_weight[p] < weights[A]);
							continue;
						double len = cur_len
							- (dist_mat[rout2[i - 1]][A]
								+ dist_mat[A][rout2[i + 1]]) * cost2;

						for (size_t j = 1; j < rout1.size(); j++)
						{
							double add_len = (dist_mat[rout1[j - 1]][A]
								+ dist_mat[A][rout1[j]]) * cost1;

							if (best_len > len + add_len)
							{
								best_len = len + add_len;
								decision = { {j,i}, 2 };
							}
						}
					}

				BSTM[p][q] = cur_len * coef_access - best_len;
				BSTM[q][p] = 0;
				RECM[p][q] = decision;
				RECM[q][p] = { {0, 0}, 5 };
			}

			void fill_BSTM_RECM(const vector<vector<size_t>>& routs)
			{
				BSTM = matrix(count_rout, vector < double>(count_rout, 0));
				RECM = vector<vector<pair<pair<size_t, size_t>, size_t>>>(count_rout,
					vector<pair<pair<size_t, size_t>, size_t>>(count_rout, { {0,0}, 5 }));
				for (size_t p = 0; p < count_rout; ++p)
					for (size_t q = p + 1; q < count_rout; ++q)
						check_change(routs, p, q);
			}

			void recalculate_BSTM_RECM(const vector<vector<size_t>>& routs, size_t p, size_t q)
			{
				for (size_t i = 0; i < count_rout; ++i)
					if (i != p)
						check_change(routs, p, i);

				for (size_t i = 0; i < count_rout; ++i)
					if (i != p && i != q)
						check_change(routs, i, q);
			}

			pair<size_t, size_t> get_max_free_BSTM(size_t iter, const vector<vector<size_t>>& routs)
			{
				auto res = pair<size_t, size_t>(0, 0);
				double max = 0;
				for (size_t p = 0; p < count_rout; ++p)
					for (size_t q = p + 1; q < count_rout; ++q)
						if (max < BSTM[p][q])
						{
							auto iter_num1 = iter - move_table[q][routs[p][RECM[p][q].first.first]];
							auto iter_num2 = iter - move_table[p][routs[q][RECM[p][q].first.second]];
							if (iter_num1 > ts && iter_num2 > ts)
							{
								max = BSTM[p][q];
								res = { p, q };
							}
						}
				return res;
			}
		};

		class Osman
		{
		public:
			Osman(const matrix& dist_mat, const vector<double>& weights,
				const vector<Transport>& transports, int_matrix routs,
				const vector<size_t> transport_id, const size_t count_rout) :
				res(routs), transport_id(transport_id),
				dist_mat(dist_mat), transports(transports), weights(weights),
				ts((size_t)(std::max(7.0, 9.6 * log(dist_mat.size() * count_rout) - 40))),
				count_rout(count_rout)// 13
			{}

			void run()
			{
				init_remain_weight();

				std::unique_ptr<checker_change> checker(new checker_change_widht_local(this));
				auto routs = res;
				auto best_length = length();
				move_table = get_move_table();

				size_t max_iter = (size_t)(340 + 0.000353 * 5 * utils::sqr((double)dist_mat.size() * count_rout));
				fill_BSTM_RECM(routs, checker);
				max_iter = 10;
				unsigned int start_time = clock();
				for (size_t i = 0; clock() - start_time < time_limit; ++i)
				{
					auto index = get_max_free_BSTM(i, routs);
					if (BSTM[index.first][index.second] == 0)
						break;
					size_t p = index.first;
					size_t q = index.second;
					auto a = RECM[p][q].first.first;
					auto b = RECM[p][q].first.second;
					auto way = RECM[p][q].second;
					switch (way)
					{
					case 0:
					{
						move_table[p][routs[p][a]] = i;
						move_table[q][routs[q][b]] = i;
						remain_weight[p] += weights[routs[p][a]] - weights[routs[q][b]];
						remain_weight[q] += weights[routs[q][b]] - weights[routs[p][a]];
						swap(routs[p][a], routs[q][b]);
						break;
					}
					case 1:
					{
						move_table[p][routs[p][a]] = i;
						remain_weight[p] += weights[routs[p][a]];
						remain_weight[q] -= weights[routs[p][a]];
						routs[q].emplace(routs[q].begin() + b, routs[p][a]);
						routs[p].erase(routs[p].begin() + a);
						break;
					}
					case 2:
					{
						move_table[q][routs[q][b]] = i;
						remain_weight[p] -= weights[routs[q][b]];
						remain_weight[q] += weights[routs[q][b]];
						routs[p].emplace(routs[p].begin() + a, routs[q][b]);
						routs[q].erase(routs[q].begin() + b);
						break;
					}
					case 3:
					{
						move_table[p][routs[p][a]] = i;
						remain_weight[p] += weights[routs[p][a]];
						remain_weight[q] -= weights[routs[p][a]];
						routs[q].push_back(routs[p][a]);
						routs[p].erase(routs[p].begin() + a);
						break;
					}
					case 4:
					{
						remain_weight[p] -= weights[routs[q][b]];
						remain_weight[q] += weights[routs[q][b]];
						move_table[q][routs[q][b]] = i;
						routs[p].push_back(routs[q][b]);
						routs[q].erase(routs[q].begin() + b);
						break;
					}
					default:
						break;
					}
					auto len = length(routs);
					if (best_length > len)
					{
						best_length = len;
						res = routs;
					}
					recalculate_BSTM_RECM(routs, p, q, checker);
				}
			}

			double length()
			{
				double length = 0;
				for (size_t i = 0; i < res.size(); ++i)
				{
					double len = utils::length_rout(res[i], dist_mat) * transports[transport_id[i]].cost_by_dist;
					if (len >= EPS)
						length += len + transports[transport_id[i]].cost_start;
				}
				return length;
			}

			double length(const int_matrix& routs)
			{
				double length = 0;
				for (size_t i = 0; i < routs.size(); ++i)
				{
					double len = utils::length_rout(routs[i], dist_mat) * transports[transport_id[i]].cost_by_dist;
					if (len >= EPS)
						length += len + transports[transport_id[i]].cost_start;
				}
				return length;
			}

			static pair<int_matrix, vector<size_t>> transform_data(const vector<int_matrix>& routs)
			{
				int_matrix res;
				vector<size_t> transport_id;

				for (size_t trans_id = 0; trans_id < routs.size(); trans_id++)
					for (size_t i = 0; i < routs[trans_id].size(); i++)
					{
						transport_id.push_back(trans_id);
						res.push_back(routs[trans_id][i]);
					}
				remove_zero_vertex(res);

				return { res , transport_id };
			}

			static void remove_zero_vertex(int_matrix& routs)
			{
				for (size_t i = 0; i < routs.size(); i++)
				{
					routs[i].pop_back();
					routs[i].erase(routs[i].begin());
				}
			}

			void add_zero_vertex()
			{
				for (size_t i = 0; i < res.size(); i++)
				{
					res[i].push_back(0);
					res[i].emplace(res[i].begin(), 0);
				}
			}

			static void add_zero_vertex(int_matrix& routs)
			{
				for (size_t i = 0; i < routs.size(); i++)
				{
					routs[i].push_back(0);
					routs[i].emplace(routs[i].begin(), 0);
				}
			}

			int_matrix res;
			const vector<size_t> transport_id;
		private:
			const matrix& dist_mat;
			const vector<Transport>& transports;
			const vector<double>& weights;
			const size_t ts;// длина списка запрета
			const size_t count_rout;
			const double coef_access = 1.05;
			const double time_limit = 40000;

			vector<vector<size_t>> BSTM;
			vector<vector<pair<pair<size_t, size_t>, size_t>>> RECM;
			vector<double> remain_weight;
			vector<vector<int>> move_table; // rout*vertex
			// вершины и способ обмена
			// 0 - обмен, 1 - левая в правую перед указанной, 2 - правая в левую перед указанной
			// 3 - левая в правую в конец, 4 - правая в левую в конец, 5 - ничего

			class checker_change
			{
			public:
				checker_change(Osman* osman) : osman(osman) {}
				void virtual check_change(const vector<vector<size_t>>& routs, size_t p, size_t q) const = 0;
				Osman* const osman;
			};

			class checker_change_widht_local : public checker_change
			{
			public:
				checker_change_widht_local(Osman* osman) : checker_change(osman) {}

				void check_change(const vector<vector<size_t>>& routs, size_t p, size_t q) const override
				{
					if (p > q)
						swap(p, q);
					double base_length = utils::length_rout(routs[p], osman->dist_mat);
					base_length += utils::length_rout(routs[q], osman->dist_mat);

					double best_length = base_length * osman->coef_access;

					auto rout1 = routs[p];
					auto rout2 = routs[q];
					pair<pair<size_t, size_t>, size_t> decision = { {0,0}, 5 };

					// обмен
					for (size_t i = 0; i < rout1.size(); ++i)
						for (size_t j = 0; j < rout2.size(); ++j)
						{
							if (osman->remain_weight[p] + osman->weights[rout1[i]] - osman->weights[rout2[j]] < 0
								|| osman->remain_weight[q] - osman->weights[rout1[i]] + osman->weights[rout2[j]] < 0)
								continue;
							swap(rout1[i], rout2[j]);
							double length = utils::length_rout(rout1, osman->dist_mat);
							length += utils::length_rout(rout2, osman->dist_mat);
							if (best_length > length)
							{
								best_length = length;
								decision = { {i,j}, 0 };
							}
							swap(rout1[i], rout2[j]);
						}

					// 1 идет в 2
					if (rout1.size() > 1)
					{
						for (size_t i = 0; i < rout1.size(); ++i)
							if (osman->remain_weight[q] - osman->weights[rout1[i]] >= 0)
								for (size_t j = 0; j < rout2.size(); ++j)
								{
									rout2.emplace(rout2.begin() + j, rout1[i]);
									rout1.erase(rout1.begin() + i);
									double length = utils::length_rout(rout1, osman->dist_mat);
									length += utils::length_rout(rout2, osman->dist_mat);
									if (best_length > length)
									{
										best_length = length;
										decision = { {i,j}, 1 };
									}
									rout1.emplace(rout1.begin() + i, rout2[j]);
									rout2.erase(rout2.begin() + j);
								}

						// вставка в конец
						for (size_t i = 0; i < rout1.size(); ++i)
						{
							if (osman->remain_weight[q] - osman->weights[rout1[i]] < 0)
								continue;
							rout2.push_back(rout1[i]);
							rout1.erase(rout1.begin() + i);
							double length = utils::length_rout(rout1, osman->dist_mat);
							length += utils::length_rout(rout2, osman->dist_mat);
							if (best_length > length)
							{
								best_length = length;
								decision = { {i,0}, 3 };
							}
							rout1.emplace(rout1.begin() + i, rout2.back());
							rout2.pop_back();
						}
					}

					// 2 идет в 1
					if (rout2.size() > 1)
					{
						for (size_t i = 0; i < rout2.size(); ++i)
							if (osman->remain_weight[p] - osman->weights[rout2[i]] >= 0)
								for (size_t j = 0; j < rout1.size(); ++j)
								{

									rout1.emplace(rout1.begin() + j, rout2[i]);
									rout2.erase(rout2.begin() + i);
									double length = utils::length_rout(rout1, osman->dist_mat);
									length += utils::length_rout(rout2, osman->dist_mat);
									if (best_length > length)
									{
										best_length = length;
										decision = { {j,i}, 2 };
									}
									rout2.emplace(rout2.begin() + i, rout1[j]);
									rout1.erase(rout1.begin() + j);
								}

						// вставка в конец
						for (size_t i = 0; i < rout2.size(); ++i)
						{
							if (osman->remain_weight[p] - osman->weights[rout2[i]] < 0)
								continue;
							rout1.push_back(rout2[i]);
							rout2.erase(rout2.begin() + i);
							double length = utils::length_rout(rout1, osman->dist_mat);
							length += utils::length_rout(rout2, osman->dist_mat);
							if (best_length > length)
							{
								best_length = length;
								decision = { {0,i}, 4 };
							}
							rout2.emplace(rout2.begin() + i, rout1.back());
							rout1.pop_back();
						}
					}

					osman->BSTM[p][q] = (size_t)(best_length - base_length);
					osman->BSTM[q][p] = 0;
					osman->RECM[p][q] = decision;
					osman->RECM[q][p] = { {0, 0}, 5 };
				}
			};

			void init_remain_weight()
			{
				remain_weight = vector<double>(res.size());
				for (size_t i = 0; i < res.size(); i++)
				{
					double weight = 0;
					for (const size_t vertex : res[i])
						weight += weights[vertex];

					remain_weight[i] = transports[transport_id[i]].capacity - weight;
				}
			}

			vector<vector<int>> get_move_table()
			{
				return vector<vector<int>>(count_rout, vector<int>(dist_mat.size(), -1000000000));
			}

			void check_change(const vector<vector<size_t>>& routs, size_t p, size_t q)
			{
				if (p > q)
					swap(p, q);
				double base_length = utils::length_rout(routs[p], dist_mat);
				base_length += utils::length_rout(routs[q], dist_mat);

				double best_length = base_length * coef_access;

				auto rout1 = routs[p];
				auto rout2 = routs[q];
				pair<pair<size_t, size_t>, size_t> decision = { {0,0}, 5 };

				// обмен
				for (size_t i = 0; i < rout1.size(); ++i)
					for (size_t j = 0; j < rout2.size(); ++j)
					{
						if (remain_weight[p] + weights[rout1[i]] - weights[rout2[j]] < 0
							|| remain_weight[q] - weights[rout1[i]] + weights[rout2[j]] < 0)
							continue;
						swap(rout1[i], rout2[j]);
						double length = utils::length_rout(rout1, dist_mat);
						length += utils::length_rout(rout2, dist_mat);
						if (best_length > length)
						{
							best_length = length;
							decision = { {i,j}, 0 };
						}
						swap(rout1[i], rout2[j]);
					}

				// 1 идет в 2
				if (rout1.size() > 1)
				{
					for (size_t i = 0; i < rout1.size(); ++i)
						if (remain_weight[q] - weights[rout1[i]] >= 0)
							for (size_t j = 0; j < rout2.size(); ++j)
							{
								rout2.emplace(rout2.begin() + j, rout1[i]);
								rout1.erase(rout1.begin() + i);
								double length = utils::length_rout(rout1, dist_mat);
								length += utils::length_rout(rout2, dist_mat);
								if (best_length > length)
								{
									best_length = length;
									decision = { {i,j}, 1 };
								}
								rout1.emplace(rout1.begin() + i, rout2[j]);
								rout2.erase(rout2.begin() + j);
							}

					// вставка в конец
					for (size_t i = 0; i < rout1.size(); ++i)
					{
						if (remain_weight[q] - weights[rout1[i]] < 0)
							continue;
						rout2.push_back(rout1[i]);
						rout1.erase(rout1.begin() + i);
						double length = utils::length_rout(rout1, dist_mat);
						length += utils::length_rout(rout2, dist_mat);
						if (best_length > length)
						{
							best_length = length;
							decision = { {i,0}, 3 };
						}
						rout1.emplace(rout1.begin() + i, rout2.back());
						rout2.pop_back();
					}
				}

				// 2 идет в 1
				if (rout2.size() > 1)
				{
					for (size_t i = 0; i < rout2.size(); ++i)
						if (remain_weight[p] - weights[rout2[i]] >= 0)
							for (size_t j = 0; j < rout1.size(); ++j)
							{

								rout1.emplace(rout1.begin() + j, rout2[i]);
								rout2.erase(rout2.begin() + i);
								double length = utils::length_rout(rout1, dist_mat);
								length += utils::length_rout(rout2, dist_mat);
								if (best_length > length)
								{
									best_length = length;
									decision = { {j,i}, 2 };
								}
								rout2.emplace(rout2.begin() + i, rout1[j]);
								rout1.erase(rout1.begin() + j);
							}

					// вставка в конец
					for (size_t i = 0; i < rout2.size(); ++i)
					{
						if (remain_weight[p] - weights[rout2[i]] < 0)
							continue;
						rout1.push_back(rout2[i]);
						rout2.erase(rout2.begin() + i);
						double length = utils::length_rout(rout1, dist_mat);
						length += utils::length_rout(rout2, dist_mat);
						if (best_length > length)
						{
							best_length = length;
							decision = { {0,i}, 4 };
						}
						rout2.emplace(rout2.begin() + i, rout1.back());
						rout1.pop_back();
					}
				}

				BSTM[p][q] = (size_t)(best_length - base_length);
				BSTM[q][p] = 0;
				RECM[p][q] = decision;
				RECM[q][p] = { {0, 0}, 5 };
			}

			void fill_BSTM_RECM(const vector<vector<size_t>>& routs, const std::unique_ptr<checker_change>& checker)
			{
				BSTM = vector<vector<size_t>>(count_rout, vector < size_t>(count_rout, 0));
				RECM = vector<vector<pair<pair<size_t, size_t>, size_t>>>(count_rout,
					vector<pair<pair<size_t, size_t>, size_t>>(count_rout, { {0,0}, 5 }));
				for (size_t p = 0; p < count_rout; ++p)
					for (size_t q = p + 1; q < count_rout; ++q)
						checker->check_change(routs, p, q);
			}

			void recalculate_BSTM_RECM(const vector<vector<size_t>>& routs, size_t p, size_t q, const std::unique_ptr<checker_change>& checker)
			{
				for (size_t i = 0; i < count_rout; ++i)
					if (i != p)
						checker->check_change(routs, p, i);

				for (size_t i = 0; i < count_rout; ++i)
					if (i != p && i != q)
						checker->check_change(routs, i, q);
			}

			pair<size_t, size_t> get_max_free_BSTM(size_t iter, const vector<vector<size_t>>& routs)
			{
				auto res = pair<size_t, size_t>(0, 0);
				double max = 0;
				for (size_t p = 0; p < count_rout; ++p)
					for (size_t q = p + 1; q < count_rout; ++q)
						if (max < BSTM[p][q])
						{
							auto iter_num1 = iter - move_table[q][routs[p][RECM[p][q].first.first]];
							auto iter_num2 = iter - move_table[p][routs[q][RECM[p][q].first.second]];
							if (iter_num1 > ts && iter_num2 > ts)
							{
								max = BSTM[p][q];
								res = { p, q };
							}
						}
				return res;
			}
		};

		class WidhtNeighborhoodSearch
		{
		public:
			WidhtNeighborhoodSearch(const matrix& dist_mat, const vector<double>& weights,
				const vector<Transport>& transports, int_matrix routs,
				const vector<size_t>& transport_id) :
				res(routs), transport_id(transport_id),
				dist_mat(dist_mat), transports(transports), weights(weights)
			{}

			void run() {
				init_remain_weight();
				double best_len = length();
				auto routs = res;

				std::random_device rd;
				std::mt19937 mersenne_rand = std::mt19937(rd());

				using ptr = std::unique_ptr<Neighborhood>;
				vector<ptr> methods;
				methods.push_back(ptr(new Opt2(this)));
				methods.push_back(ptr(new VertexShift(this)));
				methods.push_back(ptr(new EdgeShift(this)));
				methods.push_back(ptr(new VertexSwap(this)));
				methods.push_back(ptr(new EdgeSwap(this)));
				methods.push_back(ptr(new VertexRoutSwap(this)));
				methods.push_back(ptr(new VertexRoutEject(this)));

				vector<size_t> used(methods.size(), 1);
				size_t used_count = methods.size();
				used[0] += 5;
				used_count += 5;

				unsigned int start_time = clock();
				while (clock() - start_time < time_limit)
				{
					size_t rand_value = mersenne_rand() % used_count;
					size_t sum_value = 0;
					size_t meth = 0;
					for (; meth < methods.size(); ++meth)
					{
						sum_value += used[meth];
						if (sum_value > rand_value)
							break;
					}

					double length = methods[meth]->check(routs);

					if (best_len > length)
					{
						best_len = length;
						res = routs;
						++used[meth];
						++used_count;
					}
				}
				for (size_t f : used)
					std::cout << f << " ";
				std::cout << std::endl;
			}

			double length()
			{
				double length = 0;
				for (size_t i = 0; i < res.size(); ++i)
				{
					double len = utils::length_rout_0(res[i], dist_mat) * transports[transport_id[i]].cost_by_dist;
					if (len >= EPS)
						length += len + transports[transport_id[i]].cost_start;
				}
				return length;
			}

			double length(const int_matrix& routs)
			{
				double length = 0;
				for (size_t i = 0; i < routs.size(); ++i)
				{
					double len = utils::length_rout_0(routs[i], dist_mat) * transports[transport_id[i]].cost_by_dist;
					if (len >= EPS)
						length += len + transports[transport_id[i]].cost_start;
				}
				return length;
			}

			int_matrix res;
			const vector<size_t> transport_id;
		private:
			const matrix& dist_mat;
			const vector<Transport>& transports;
			const vector<double>& weights;
			const double coef_access = 1.09;
			const double time_limit = 20000;
			vector<double> remain_weight;

			void init_remain_weight()
			{
				remain_weight = vector<double>(res.size());
				for (size_t i = 0; i < res.size(); i++)
				{
					double weight = 0;
					for (const size_t vertex : res[i])
						weight += weights[vertex];

					remain_weight[i] = transports[transport_id[i]].capacity - weight;
				}
			}

			class Neighborhood
			{
			public:
				Neighborhood(WidhtNeighborhoodSearch* data) : data(data) {}
				virtual double check(int_matrix& routs) = 0;
				WidhtNeighborhoodSearch* const data;
			};

			class Opt2 : public Neighborhood
			{
			public:
				Opt2(WidhtNeighborhoodSearch* data) : Neighborhood(data) {}
				double check(int_matrix& routs) override {
					double length = 0;

					for (size_t i = 0; i < routs.size(); i++)
						TSP::local_opt::opt_2_fast2(
							routs[i],
							data->dist_mat);
					return data->length(routs);
				}
			};

			class VertexShift : public Neighborhood
			{
			public:
				VertexShift(WidhtNeighborhoodSearch* data) : Neighborhood(data) {}
				double check(int_matrix& routs) override {

					for (size_t r = 0; r < routs.size(); r++)
					{
						vector<size_t>& rout = routs[r];

						double cur_len = utils::length_rout(rout, data->dist_mat);
						double best_len = cur_len * data->coef_access;

						if (rout.size() < 4)
							continue;
						// intial and final position are fixed (initial/final node remains 0)
						for (size_t i = 1; i < rout.size() - 2; i++)
						{
							size_t A = rout[i];
							double len = cur_len
								- data->dist_mat[rout[i - 1]][A]
								- data->dist_mat[A][rout[i + 1]];
							size_t to = i;

							for (size_t j = i + 1; j < rout.size() - 1; j++)
							{
								double add_len = data->dist_mat[rout[j]][A]
									+ data->dist_mat[A][rout[j + 1]];

								if (best_len > len + add_len)
								{
									best_len = len + add_len;
									to = j;
								}
							}

							if (to != i)
							{
								rout.erase(rout.begin() + i);
								rout.emplace(rout.begin() + to, A);
								cur_len = best_len;
							}
						}
					}

					return data->length(routs);
				}
			};

			class EdgeShift : public Neighborhood
			{
			public:
				EdgeShift(WidhtNeighborhoodSearch* data) : Neighborhood(data) {}
				double check(int_matrix& routs) override {

					for (size_t r = 0; r < routs.size(); r++)
					{
						vector<size_t>& rout = routs[r];

						double cur_len = utils::length_rout(rout, data->dist_mat);
						double best_len = cur_len * data->coef_access;

						if (rout.size() < 5)
							continue;
						// intial and final position are fixed (initial/final node remains 0)
						for (size_t i = 1; i < rout.size() - 3; i++)
						{
							size_t A1 = rout[i];
							size_t A2 = rout[i + 1];
							double len = cur_len
								- data->dist_mat[rout[i - 1]][A1]
								- data->dist_mat[A2][rout[i + 2]];
							size_t to = i;

							for (size_t j = i + 2; j < rout.size() - 1; j++)
							{
								double add_len = data->dist_mat[rout[j]][A1]
									+ data->dist_mat[A2][rout[j + 1]];

								if (best_len > len + add_len)
								{
									best_len = len + add_len;
									to = j;
								}
							}

							if (to != i)
							{
								rout.erase(rout.begin() + i);
								rout.erase(rout.begin() + i);
								rout.emplace(rout.begin() + to - 1, A2);
								rout.emplace(rout.begin() + to - 1, A1);
								cur_len = best_len;
							}
						}
					}

					return data->length(routs);
				}
			};

			class VertexSwap : public Neighborhood
			{
			public:
				VertexSwap(WidhtNeighborhoodSearch* data) : Neighborhood(data) {}
				double check(int_matrix& routs) override {

					for (size_t r = 0; r < routs.size(); r++)
					{
						vector<size_t>& rout = routs[r];

						double cur_len = utils::length_rout(rout, data->dist_mat);
						double best_len = cur_len * data->coef_access;

						if (rout.size() < 4)
							continue;
						// intial and final position are fixed (initial/final node remains 0)
						for (size_t i = 1; i < rout.size() - 2; i++)
						{
							size_t A = rout[i];
							double len = cur_len
								- data->dist_mat[rout[i - 1]][A]
								- data->dist_mat[A][rout[i + 1]];
							size_t to = i;

							for (size_t j = i + 2; j < rout.size() - 1; j++)
							{

								size_t B = rout[j];
								double add_len = data->dist_mat[rout[j - 1]][A]
									+ data->dist_mat[A][rout[j + 1]]
									- data->dist_mat[rout[j - 1]][B]
									- data->dist_mat[B][rout[j + 1]]
									+ data->dist_mat[rout[i - 1]][B]
									+ data->dist_mat[B][rout[i + 1]];

								if (best_len > len + add_len)
								{
									best_len = len + add_len;
									to = j;
								}
							}

							if (to != i)
							{
								swap(rout[i], rout[to]);
								cur_len = best_len;
							}
						}
					}

					return data->length(routs);
				}
			};

			class EdgeSwap : public Neighborhood
			{
			public:
				EdgeSwap(WidhtNeighborhoodSearch* data) : Neighborhood(data) {}
				double check(int_matrix& routs) override {

					for (size_t r = 0; r < routs.size(); r++)
					{
						vector<size_t>& rout = routs[r];

						double cur_len = utils::length_rout(rout, data->dist_mat);
						double best_len = cur_len * data->coef_access;

						if (rout.size() < 6)
							continue;
						// intial and final position are fixed (initial/final node remains 0)
						for (size_t i = 1; i < rout.size() - 3; i++)
						{
							size_t A1 = rout[i];
							size_t A2 = rout[i + 1];
							double len = cur_len
								- data->dist_mat[rout[i - 1]][A1]
								- data->dist_mat[A2][rout[i + 2]];
							size_t to = i;

							for (size_t j = i + 3; j < rout.size() - 2; j++)
							{

								size_t B1 = rout[j];
								size_t B2 = rout[j + 1];
								double add_len = data->dist_mat[rout[j - 1]][A1]
									+ data->dist_mat[A2][rout[j + 1]]
									- data->dist_mat[rout[j - 1]][B1]
									- data->dist_mat[B2][rout[j + 1]]
									+ data->dist_mat[rout[i - 1]][B1]
									+ data->dist_mat[B2][rout[i + 1]];

								if (best_len > len + add_len)
								{
									best_len = len + add_len;
									to = j;
								}
							}

							if (to != i)
							{
								swap(rout[i], rout[to]);
								swap(rout[i + 1], rout[to + 1]);
								cur_len = best_len;
							}
						}
					}

					return data->length(routs);
				}
			};

			class VertexRoutSwap : public Neighborhood
			{
			public:
				VertexRoutSwap(WidhtNeighborhoodSearch* data) : Neighborhood(data) {}
				double check(int_matrix& routs) override {
					for (size_t r1 = 0; r1 < routs.size(); r1++)
						for (size_t r2 = r1 + 1; r2 < routs.size(); r2++)
						{
							vector<size_t>& rout1 = routs[r1];
							vector<size_t>& rout2 = routs[r2];
							double cost1 = data->transports[data->transport_id[r1]].cost_by_dist;
							double cost2 = data->transports[data->transport_id[r2]].cost_by_dist;
							double cur_len = utils::length_rout(rout1, data->dist_mat)
								* cost1
								+ utils::length_rout(rout2, data->dist_mat)
								* cost2;
							double best_len = cur_len * data->coef_access;
							// intial and final position are fixed (initial/final node remains 0)
							for (size_t i = 1; i < rout1.size() - 1; i++)
							{
								size_t A = rout1[i];
								double len = cur_len
									- (data->dist_mat[rout1[i - 1]][A]
									+ data->dist_mat[A][rout1[i + 1]]) * cost1;
								size_t to = 0;

								for (size_t j = 1; j < rout2.size() - 1; j++)
								{
									size_t B = rout2[j];
									if (data->remain_weight[r1] + data->weights[A] - data->weights[B] < 0
										|| data->remain_weight[r2] - data->weights[A] + data->weights[B] < 0)
										continue;
									double add_len = (data->dist_mat[rout2[j - 1]][A]
										+ data->dist_mat[A][rout2[j + 1]]) * cost2
										- (data->dist_mat[rout2[j - 1]][B]
										+ data->dist_mat[B][rout2[j + 1]]) * cost2
										+ (data->dist_mat[rout1[i - 1]][B]
										+ data->dist_mat[B][rout1[i + 1]]) * cost1;

									if (best_len > len + add_len)
									{
										best_len = len + add_len;
										to = j;
									}
								}

								if (to != 0)
								{
									data->remain_weight[r1] += data->weights[A] - data->weights[rout2[to]] ;
									data->remain_weight[r2] += data->weights[rout2[to]] - data->weights[A];
									swap(rout1[i], rout2[to]);
									cur_len = best_len;
								}
							}

						}

					return data->length(routs);
				}
			};

			class VertexRoutEject : public Neighborhood
			{
			public:
				VertexRoutEject(WidhtNeighborhoodSearch* data) : Neighborhood(data) {}
				double check(int_matrix& routs) override {
					for (size_t r1 = 0; r1 < routs.size(); r1++)
						for (size_t r2 = 0; r2 < routs.size(); r2++)
						{
							if (r1 == r2)
								continue;
							vector<size_t>& rout1 = routs[r1];
							vector<size_t>& rout2 = routs[r2];

							if (rout2.size() < 3)
								continue;

							double cost1 = data->transports[data->transport_id[r1]].cost_by_dist;
							double cost2 = data->transports[data->transport_id[r2]].cost_by_dist;
							double cur_len = utils::length_rout(rout1, data->dist_mat)
								* cost1
								+ utils::length_rout(rout2, data->dist_mat)
								* cost2;
							double best_len = cur_len * data->coef_access;
							// intial and final position are fixed (initial/final node remains 0)
							for (size_t i = 1; i < rout1.size() - 1; i++)
							{
								size_t A = rout1[i];
								if (data->remain_weight[r2] < data->weights[A]);
									continue;
								double len = cur_len
									- (data->dist_mat[rout1[i - 1]][A]
										+ data->dist_mat[A][rout1[i + 1]]) * cost1;
								size_t to = 0;

								for (size_t j = 1; j < rout2.size(); j++)
								{

									double add_len = (data->dist_mat[rout2[j - 1]][A]
										+ data->dist_mat[A][rout2[j]]) * cost2;

									if (best_len > len + add_len)
									{
										best_len = len + add_len;
										to = j;
									}
								}

								if (to != 0)
								{
									data->remain_weight[r1] += data->weights[A];
									data->remain_weight[r2] -= data->weights[A];
									rout1.erase(rout1.begin() + i);
									rout2.emplace(rout2.begin() + to, A);
									cur_len = best_len;
								}
							}

						}

					return data->length(routs);
				}
			};
		};

		class GlobalWidhtNeighborhoodSearch
		{
		public:
			GlobalWidhtNeighborhoodSearch(const matrix& dist_mat, const vector<double>& weights,
				const vector<Transport>& transports, int_matrix routs,
				const vector<size_t>& transport_id,
				const vector<size_t>& cluster_id, const vector<size_t>& frequence) :
				res(routs), transport_id(transport_id),
				cluster_id(cluster_id), frequence(frequence),
				dist_mat(dist_mat), transports(transports), weights(weights)
			{}

			void run() {
				init_remain_weight();
				double best_len = length();
				auto routs = res;

				std::random_device rd;
				std::mt19937 mersenne_rand = std::mt19937(rd());

				using ptr = std::unique_ptr<Neighborhood>;
				vector<ptr> methods;
				methods.push_back(ptr(new Opt2(this)));
				methods.push_back(ptr(new VertexShift(this)));
				methods.push_back(ptr(new EdgeShift(this)));
				methods.push_back(ptr(new VertexSwap(this)));
				methods.push_back(ptr(new EdgeSwap(this)));
				methods.push_back(ptr(new VertexRoutSwap(this)));
				methods.push_back(ptr(new VertexRoutEject(this)));

				vector<size_t> used{5, 1, 1, 1, 1, 3, 3 };
				size_t used_count = 15;

				unsigned int start_time = clock();
				while (clock() - start_time < time_limit)
				{
					size_t rand_value = mersenne_rand() % used_count;
					size_t sum_value = 0;
					size_t meth = 0;
					for (; meth < methods.size(); ++meth)
					{
						sum_value += used[meth];
						if (sum_value > rand_value)
							break;
					}

					double length = methods[meth]->check(routs);

					if (best_len > length)
					{
						best_len = length;
						res = routs;
						++used[meth];
						++used_count;
					}
				}
				for (size_t f : used)
					std::cout << f << " ";
				std::cout << std::endl;
			}

			double length()
			{
				double length = 0;
				for (size_t i = 0; i < res.size(); ++i)
				{
					double len = utils::length_rout_0(res[i], dist_mat) * transports[transport_id[i]].cost_by_dist;
					if (len >= EPS)
						length += len + transports[transport_id[i]].cost_start;
				}
				return length;
			}

			double length(const int_matrix& routs)
			{
				double length = 0;
				for (size_t i = 0; i < routs.size(); ++i)
				{
					double len = utils::length_rout_0(routs[i], dist_mat) * transports[transport_id[i]].cost_by_dist;
					if (len >= EPS)
						length += len + transports[transport_id[i]].cost_start;
				}
				return length;
			}

			int_matrix res;
			const vector<size_t> transport_id;
			const vector<size_t>& cluster_id;
			const vector<size_t>& frequence;
		private:
			const matrix& dist_mat;
			const vector<Transport>& transports;
			const vector<double>& weights;
			const double coef_access = 1.09;
			const double time_limit = 120000;
			vector<double> remain_weight;

			void init_remain_weight()
			{
				remain_weight = vector<double>(res.size());
				for (size_t i = 0; i < res.size(); i++)
				{
					double weight = 0;
					for (const size_t vertex : res[i])
						weight += weights[vertex];

					remain_weight[i] = transports[transport_id[i]].capacity - weight;
				}
			}

			class Neighborhood
			{
			public:
				Neighborhood(GlobalWidhtNeighborhoodSearch* data) : data(data) {}
				virtual double check(int_matrix& routs) = 0;
				GlobalWidhtNeighborhoodSearch* const data;
			};

			class Opt2 : public Neighborhood
			{
			public:
				Opt2(GlobalWidhtNeighborhoodSearch* data) : Neighborhood(data) {}
				double check(int_matrix& routs) override {
					double length = 0;

					for (size_t i = 0; i < routs.size(); i++)
						TSP::local_opt::opt_2_fast2(
							routs[i],
							data->dist_mat);
					return data->length(routs);
				}
			};

			class VertexShift : public Neighborhood
			{
			public:
				VertexShift(GlobalWidhtNeighborhoodSearch* data) : Neighborhood(data) {}
				double check(int_matrix& routs) override {

					for (size_t r = 0; r < routs.size(); r++)
					{
						vector<size_t>& rout = routs[r];

						double cur_len = utils::length_rout(rout, data->dist_mat);
						double best_len = cur_len * data->coef_access;

						if (rout.size() < 4)
							continue;
						// intial and final position are fixed (initial/final node remains 0)
						for (size_t i = 1; i < rout.size() - 2; i++)
						{
							size_t A = rout[i];
							double len = cur_len
								- data->dist_mat[rout[i - 1]][A]
								- data->dist_mat[A][rout[i + 1]];
							size_t to = i;

							for (size_t j = i + 1; j < rout.size() - 1; j++)
							{
								double add_len = data->dist_mat[rout[j]][A]
									+ data->dist_mat[A][rout[j + 1]];

								if (best_len > len + add_len)
								{
									best_len = len + add_len;
									to = j;
								}
							}

							if (to != i)
							{
								rout.erase(rout.begin() + i);
								rout.emplace(rout.begin() + to, A);
								cur_len = best_len;
							}
						}
					}

					return data->length(routs);
				}
			};

			class EdgeShift : public Neighborhood
			{
			public:
				EdgeShift(GlobalWidhtNeighborhoodSearch* data) : Neighborhood(data) {}
				double check(int_matrix& routs) override {

					for (size_t r = 0; r < routs.size(); r++)
					{
						vector<size_t>& rout = routs[r];

						double cur_len = utils::length_rout(rout, data->dist_mat);
						double best_len = cur_len * data->coef_access;

						if (rout.size() < 5)
							continue;
						// intial and final position are fixed (initial/final node remains 0)
						for (size_t i = 1; i < rout.size() - 3; i++)
						{
							size_t A1 = rout[i];
							size_t A2 = rout[i + 1];
							double len = cur_len
								- data->dist_mat[rout[i - 1]][A1]
								- data->dist_mat[A2][rout[i + 2]];
							size_t to = i;

							for (size_t j = i + 2; j < rout.size() - 1; j++)
							{
								double add_len = data->dist_mat[rout[j]][A1]
									+ data->dist_mat[A2][rout[j + 1]];

								if (best_len > len + add_len)
								{
									best_len = len + add_len;
									to = j;
								}
							}

							if (to != i)
							{
								rout.erase(rout.begin() + i);
								rout.erase(rout.begin() + i);
								rout.emplace(rout.begin() + to - 1, A2);
								rout.emplace(rout.begin() + to - 1, A1);
								cur_len = best_len;
							}
						}
					}

					return data->length(routs);
				}
			};

			class VertexSwap : public Neighborhood
			{
			public:
				VertexSwap(GlobalWidhtNeighborhoodSearch* data) : Neighborhood(data) {}
				double check(int_matrix& routs) override {

					for (size_t r = 0; r < routs.size(); r++)
					{
						vector<size_t>& rout = routs[r];

						double cur_len = utils::length_rout(rout, data->dist_mat);
						double best_len = cur_len * data->coef_access;

						if (rout.size() < 4)
							continue;
						// intial and final position are fixed (initial/final node remains 0)
						for (size_t i = 1; i < rout.size() - 2; i++)
						{
							size_t A = rout[i];
							double len = cur_len
								- data->dist_mat[rout[i - 1]][A]
								- data->dist_mat[A][rout[i + 1]];
							size_t to = i;

							for (size_t j = i + 2; j < rout.size() - 1; j++)
							{

								size_t B = rout[j];
								double add_len = data->dist_mat[rout[j - 1]][A]
									+ data->dist_mat[A][rout[j + 1]]
									- data->dist_mat[rout[j - 1]][B]
									- data->dist_mat[B][rout[j + 1]]
									+ data->dist_mat[rout[i - 1]][B]
									+ data->dist_mat[B][rout[i + 1]];

								if (best_len > len + add_len)
								{
									best_len = len + add_len;
									to = j;
								}
							}

							if (to != i)
							{
								swap(rout[i], rout[to]);
								cur_len = best_len;
							}
						}
					}

					return data->length(routs);
				}
			};

			class EdgeSwap : public Neighborhood
			{
			public:
				EdgeSwap(GlobalWidhtNeighborhoodSearch* data) : Neighborhood(data) {}
				double check(int_matrix& routs) override {

					for (size_t r = 0; r < routs.size(); r++)
					{
						vector<size_t>& rout = routs[r];

						double cur_len = utils::length_rout(rout, data->dist_mat);
						double best_len = cur_len * data->coef_access;

						if (rout.size() < 6)
							continue;
						// intial and final position are fixed (initial/final node remains 0)
						for (size_t i = 1; i < rout.size() - 3; i++)
						{
							size_t A1 = rout[i];
							size_t A2 = rout[i + 1];
							double len = cur_len
								- data->dist_mat[rout[i - 1]][A1]
								- data->dist_mat[A2][rout[i + 2]];
							size_t to = i;

							for (size_t j = i + 3; j < rout.size() - 2; j++)
							{

								size_t B1 = rout[j];
								size_t B2 = rout[j + 1];
								double add_len = data->dist_mat[rout[j - 1]][A1]
									+ data->dist_mat[A2][rout[j + 1]]
									- data->dist_mat[rout[j - 1]][B1]
									- data->dist_mat[B2][rout[j + 1]]
									+ data->dist_mat[rout[i - 1]][B1]
									+ data->dist_mat[B2][rout[i + 1]];

								if (best_len > len + add_len)
								{
									best_len = len + add_len;
									to = j;
								}
							}

							if (to != i)
							{
								swap(rout[i], rout[to]);
								swap(rout[i + 1], rout[to + 1]);
								cur_len = best_len;
							}
						}
					}

					return data->length(routs);
				}
			};

			class VertexRoutSwap : public Neighborhood
			{
			public:
				VertexRoutSwap(GlobalWidhtNeighborhoodSearch* data) : Neighborhood(data) {}
				double check(int_matrix& routs) override {
					for (size_t r1 = 0; r1 < routs.size(); r1++)
						for (size_t r2 = r1 + 1; r2 < routs.size(); r2++)
						{
							if (data->cluster_id[r1] != data->cluster_id[r2])
							{
								check_different_clust(routs, r1, r2);
								continue;
							}
							vector<size_t>& rout1 = routs[r1];
							vector<size_t>& rout2 = routs[r2];
							double cost1 = data->transports[data->transport_id[r1]].cost_by_dist;
							double cost2 = data->transports[data->transport_id[r2]].cost_by_dist;
							double cur_len = utils::length_rout(rout1, data->dist_mat)
								* cost1
								+ utils::length_rout(rout2, data->dist_mat)
								* cost2;
							double best_len = cur_len * data->coef_access;
							// intial and final position are fixed (initial/final node remains 0)
							for (size_t i = 1; i < rout1.size() - 1; i++)
							{
								size_t A = rout1[i];
								double len = cur_len
									- (data->dist_mat[rout1[i - 1]][A]
										+ data->dist_mat[A][rout1[i + 1]]) * cost1;
								size_t to = 0;

								for (size_t j = 1; j < rout2.size() - 1; j++)
								{
									size_t B = rout2[j];
									if (data->remain_weight[r1] + data->weights[A] - data->weights[B] < 0
										|| data->remain_weight[r2] - data->weights[A] + data->weights[B] < 0)
										continue;
									double add_len = (data->dist_mat[rout2[j - 1]][A]
										+ data->dist_mat[A][rout2[j + 1]]) * cost2
										- (data->dist_mat[rout2[j - 1]][B]
											+ data->dist_mat[B][rout2[j + 1]]) * cost2
										+ (data->dist_mat[rout1[i - 1]][B]
											+ data->dist_mat[B][rout1[i + 1]]) * cost1;

									if (best_len > len + add_len)
									{
										best_len = len + add_len;
										to = j;
									}
								}

								if (to != 0)
								{
									data->remain_weight[r1] += data->weights[A] - data->weights[rout2[to]];
									data->remain_weight[r2] += data->weights[rout2[to]] - data->weights[A];
									swap(rout1[i], rout2[to]);
									cur_len = best_len;
								}
							}

						}

					return data->length();;
				}

				void check_different_clust(int_matrix& routs, const size_t r1, const size_t r2)
				{
					vector<size_t>& rout1 = routs[r1];
					vector<size_t>& rout2 = routs[r2];
					double cost1 = data->transports[data->transport_id[r1]].cost_by_dist;
					double cost2 = data->transports[data->transport_id[r2]].cost_by_dist;
					double cur_len = utils::length_rout(rout1, data->dist_mat)
						* cost1
						+ utils::length_rout(rout2, data->dist_mat)
						* cost2;
					double best_len = cur_len * data->coef_access;
					// intial and final position are fixed (initial/final node remains 0)
					for (size_t i = 1; i < rout1.size() - 1; i++)
					{
						size_t A = rout1[i];
						if (data->frequence[A] != 1)
							continue;
						double len = cur_len
							- (data->dist_mat[rout1[i - 1]][A]
								+ data->dist_mat[A][rout1[i + 1]]) * cost1;
						size_t to = 0;

						for (size_t j = 1; j < rout2.size() - 1; j++)
						{
							size_t B = rout2[j];
							if (data->frequence[B] != 1
								|| data->remain_weight[r1] + data->weights[A] - data->weights[B] < 0
								|| data->remain_weight[r2] - data->weights[A] + data->weights[B] < 0)
								continue;
							double add_len = (data->dist_mat[rout2[j - 1]][A]
								+ data->dist_mat[A][rout2[j + 1]]) * cost2
								- (data->dist_mat[rout2[j - 1]][B]
									+ data->dist_mat[B][rout2[j + 1]]) * cost2
								+ (data->dist_mat[rout1[i - 1]][B]
									+ data->dist_mat[B][rout1[i + 1]]) * cost1;

							if (best_len > len + add_len)
							{
								best_len = len + add_len;
								to = j;
							}
						}

						if (to != 0)
						{
							data->remain_weight[r1] += data->weights[A] - data->weights[rout2[to]];
							data->remain_weight[r2] += data->weights[rout2[to]] - data->weights[A];
							swap(rout1[i], rout2[to]);
							cur_len = best_len;
						}
					}
				}
			};

			class VertexRoutEject : public Neighborhood
			{
			public:
				VertexRoutEject(GlobalWidhtNeighborhoodSearch* data) : Neighborhood(data) {}
				double check(int_matrix& routs) override {
					for (size_t r1 = 0; r1 < routs.size(); r1++)
						for (size_t r2 = 0; r2 < routs.size(); r2++)
						{
							if (r1 == r2)
								continue;
							if (data->cluster_id[r1] != data->cluster_id[r2])
							{
								check_different_clust(routs, r1, r2);
								continue;
							}
							vector<size_t>& rout1 = routs[r1];
							vector<size_t>& rout2 = routs[r2];

							if (rout2.size() < 3)
								continue;

							double cost1 = data->transports[data->transport_id[r1]].cost_by_dist;
							double cost2 = data->transports[data->transport_id[r2]].cost_by_dist;
							double cur_len = utils::length_rout(rout1, data->dist_mat)
								* cost1
								+ utils::length_rout(rout2, data->dist_mat)
								* cost2;
							double best_len = cur_len * data->coef_access;
							// intial and final position are fixed (initial/final node remains 0)
							for (size_t i = 1; i < rout1.size() - 1; i++)
							{
								size_t A = rout1[i];
								if (data->remain_weight[r2] < data->weights[A]);
									continue;
								double len = cur_len
									- (data->dist_mat[rout1[i - 1]][A]
										+ data->dist_mat[A][rout1[i + 1]]) * cost1;
								size_t to = 0;

								for (size_t j = 1; j < rout2.size(); j++)
								{

									double add_len = (data->dist_mat[rout2[j - 1]][A]
										+ data->dist_mat[A][rout2[j]]) * cost2;

									if (best_len > len + add_len)
									{
										best_len = len + add_len;
										to = j;
									}
								}

								if (to != 0)
								{
									data->remain_weight[r1] += data->weights[A];
									data->remain_weight[r2] -= data->weights[A];
									rout1.erase(rout1.begin() + i);
									rout2.emplace(rout2.begin() + to, A);
									cur_len = best_len;
								}
							}

						}

					return data->length();;
				}

				void check_different_clust(int_matrix& routs, const size_t r1, const size_t r2)
				{
					vector<size_t>& rout1 = routs[r1];
					vector<size_t>& rout2 = routs[r2];

					if (rout2.size() < 3)
						return;

					double cost1 = data->transports[data->transport_id[r1]].cost_by_dist;
					double cost2 = data->transports[data->transport_id[r2]].cost_by_dist;
					double cur_len = utils::length_rout(rout1, data->dist_mat)
						* cost1
						+ utils::length_rout(rout2, data->dist_mat)
						* cost2;
					double best_len = cur_len * data->coef_access;
					// intial and final position are fixed (initial/final node remains 0)
					for (size_t i = 1; i < rout1.size() - 1; i++)
					{
						size_t A = rout1[i];
						if (data->frequence[A] != 1
							|| data->remain_weight[r2] < data->weights[A])
							continue;
						double len = cur_len
							- (data->dist_mat[rout1[i - 1]][A]
								+ data->dist_mat[A][rout1[i + 1]]) * cost1;
						size_t to = 0;

						for (size_t j = 1; j < rout2.size(); j++)
						{

							double add_len = (data->dist_mat[rout2[j - 1]][A]
								+ data->dist_mat[A][rout2[j]]) * cost2;

							if (best_len > len + add_len)
							{
								best_len = len + add_len;
								to = j;
							}
						}

						if (to != 0)
						{
							data->remain_weight[r1] -= data->weights[A];
							data->remain_weight[r2] += data->weights[A];
							rout1.erase(rout1.begin() + i);
							rout2.emplace(rout2.begin() + to, A);
							cur_len = best_len;
						}
					}
				}
			};
		};
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

	class Ant_algorithm
	{
	public:
		Ant_algorithm(const matrix& dist_mat, const vector<double>& weights, const vector<Transport>& transports);

		struct Ant
		{
			vector<size_t> rout;
			double remain_volume;
			size_t transport_type;
			size_t current_pos;
			bool is_last; // для определения нужно ли после его заполнения пускать нового
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

		void run()
		{
			double best_length = UINT32_MAX;
			vector<int_matrix> best_res;

			double min_weight = min_weight_vertex();
			double max_weight = max_weight_vertex();

			unsigned int start_time = clock();
			//unsigned int cur_time = start_time;
			//unsigned int prev_cur_time = start_time;
			while (true)
			{
				fill_pheromone();

				vector<int_matrix> routs(transports.size());
				for (size_t iter = 1; iter < count_iter_inner; ++iter) {

					Calc_data data = fill_data();
					vector<int_matrix> cur_res(transports.size());

					for (size_t used = 1; used < weights.size(); ++used)
					{
						auto trans_vertex = get_next_step(data);

						if (trans_vertex.second == UINT32_MAX)
						{
							get_next_step_sharing(data);
						}

						size_t trans = trans_vertex.first;
						size_t vertex = data.vertex_used[trans_vertex.second];
						data.ants[trans].current_pos = vertex;
						data.ants[trans].remain_volume -= weights[vertex];
						data.ants[trans].rout.push_back(vertex);
						data.vertex_used.erase(data.vertex_used.begin() + trans_vertex.second);

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
							cur_res[transport_type].push_back(data.ants[trans].rout);
							data.ants.erase(data.ants.begin() + trans);
						}
					}
					for (Ant& ant: data.ants)
					{
						ant.rout.push_back(0);
						routs[ant.transport_type].push_back(ant.rout);
						cur_res[ant.transport_type].push_back(ant.rout);
					}

					double len= length(cur_res);

					if (best_length > len) {
						best_length = len;
						best_res = cur_res;
					}
				}
				// 130-131
				calcalate_pheromone(routs);
				//unsigned int diff = cur_time - prev_cur_time;
				//prev_cur_time = cur_time;
				//std::cout << diff << std::endl;
				if (clock() - start_time > time_limit)
					break;
			}

			for (size_t i = 0; i < best_res.size(); ++i)
				for (size_t j = 0; j < best_res[i].size(); ++j)
					TSP::local_opt::opt_2_fast2(best_res[i][j], dist_mat);

			res = best_res;
		}

		double length(const vector<int_matrix>& routs)
		{
			double length = 0;
			for (size_t trans_type = 0; trans_type < routs.size(); ++trans_type)
				for (const vector<size_t>& rout : routs[trans_type])
					length += utils::length_rout_0(rout, dist_mat) * transports[trans_type].cost_by_dist + transports[trans_type].cost_start;

			return length;
		}

		double length()
		{
			double length = 0;
			for (size_t trans_type = 0; trans_type < res.size(); ++trans_type)
				for (const vector<size_t>& rout : res[trans_type])
					length += utils::length_rout_0(rout, dist_mat) * transports[trans_type].cost_by_dist + transports[trans_type].cost_start;

			return length;
		}

		vector<int_matrix> res;
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
		const size_t count_iter_inner = 10;
		const double time_limit = 40000;

		const double evaporation_rate = 0.8; // коэфициент испарения
		//const double addition_pheromone_rate = 0.8; // коэфициент добавленяи феромона
		const double Q_rate = 5; // коэфициент добавленяи феромона
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
				vector<size_t>(dist_mat.size() - 1, 0)
			};

			for (size_t i = 0; i < transports.size(); ++i)
			{
				data.transport_remain[i] = transports[i].count - 1;
				data.ants.push_back({ vector<size_t>(1, 0), transports[i].capacity, i, 0 , true });
			}

			for (size_t i = 0; i < data.vertex_used.size(); ++i)
			{
				data.vertex_used[i] = i + 1;
			}

			//data.vertex_used[0] = 1;
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
				for (size_t i = 0; i < data.vertex_used.size(); ++i)
				{
					size_t vertex = data.vertex_used[i];
					if (ant.remain_volume > weights[vertex]) {

						

						if (dist_mat[pos][vertex] > EPS)
							cost_move[trans][i] = pow(pheromone_mat[ant.transport_type][pos][vertex], alpha)
							* pow(1 / (dist_mat[pos][vertex] * transports[ant.transport_type].cost_by_dist), beta);
						else
						{
							return { trans , i };
						}
						sum += cost_move[trans][i];
					}
				}
			}

			auto rand_value = (double)mersenne_rand() * sum / UINT32_MAX;
			double sum_after = 0;

			for (size_t trans = 0; trans < data.ants.size(); ++trans)
				for (size_t i = 0; i < data.vertex_used.size(); ++i)
				{
					sum_after += cost_move[trans][i];
					if (sum_after > rand_value)
						return { trans , i };
				}

			return { UINT32_MAX , UINT32_MAX };
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
	};
}