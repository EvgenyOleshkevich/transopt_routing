#include "../headers/balancedVRP.hpp"
#include "../headers/TSP.hpp"
#include <math.h>

namespace balancedVRP
{
	namespace clustering
	{
		std::pair<int, int> two_farthest_vertex(const matrix& dist_mat,
			const std::vector<size_t>& vertexes) {
			double max_dist = 0;
			std::pair<int, int> res = { 0, 1 };
			for (size_t i = 0; i < vertexes.size(); ++i)
				for (size_t j = i + 1; j < vertexes.size(); ++j)
					if (max_dist < dist_mat[vertexes[i]][vertexes[j]])
					{
						max_dist = dist_mat[vertexes[i]][vertexes[j]];
						res.first = i;
						res.second = j;
					}
			return res;
		}

		void update_dist_to_cluaster(const matrix& dist_mat,
			const std::vector<size_t>& vertexes, std::vector<double>& dist_to_cluaster,
			const int first, const int second)
		{
			for (int i = 0; i < vertexes.size(); ++i)
				if (dist_to_cluaster[i] < INT_MAX)
					dist_to_cluaster[i] = dist_mat[vertexes[first]][vertexes[i]]
					- dist_mat[vertexes[second]][vertexes[i]];
		}

		int_matrix dichotomous_division(const matrix& dist_mat, const int clusters)
		{
			std::vector<size_t> vertexes(dist_mat.size() - 1);
			for (size_t i = 1; i < dist_mat.size(); i++)
				vertexes[i - 1] = i;
			return dichotomous_division(dist_mat, vertexes, clusters, vertexes.size() / clusters);
		}


		int_matrix dichotomous_division(const matrix& dist_mat,
			const std::vector<size_t>& vertexes, const int clusters, const int cluster_size)
		{
			int k1 = clusters / 2;
			int k2 = clusters / 2 + clusters % 2;

			std::vector<double> dist_to_cluaster(vertexes.size());
			for (size_t i = 0; i < dist_to_cluaster.size(); i++)
				dist_to_cluaster[i] = (double)INT_MAX - 1;
			auto pair = two_farthest_vertex(dist_mat, vertexes);
			int first = pair.first;
			int second = pair.second;

			dist_to_cluaster[first] = (double)INT_MAX + 1;
			dist_to_cluaster[second] = (double)INT_MAX + 1;

			std::vector<size_t> first_cluster = { vertexes[first] };
			std::vector<size_t> second_cluster = { vertexes[second] };

			while (vertexes.size() - first_cluster.size() > k2 * (cluster_size + 1) ||
				first_cluster.size() < k1 * cluster_size)
			{
				/* inline body of : update_dist_to_cluaster
				for (int i = 0; i < vertexes.size(); ++i)
					if (dist_to_cluaster[i] < INT_MAX)
						dist_to_cluaster[i] = dist_mat[vertexes[first]][vertexes[i]]
						- dist_mat[vertexes[second]][vertexes[i]];
				//update_dist_to_cluaster(dist_mat, vertexes,
					//dist_to_cluaster, first, second);
				*/
				int min = INT_MAX;
				int max = INT_MIN;
				int min_ind = -1;
				int max_ind = -1;

				for (int i = 0; i < dist_to_cluaster.size(); ++i)
				{
					// inline body of : update_dist_to_cluaster
					if (dist_to_cluaster[i] < INT_MAX)
					{
						dist_to_cluaster[i] = dist_mat[vertexes[first]][vertexes[i]]
							- dist_mat[vertexes[second]][vertexes[i]];
						//

						if (min > dist_to_cluaster[i])
						{
							min = dist_to_cluaster[i];
							min_ind = i;
						}
						if (max < dist_to_cluaster[i] && dist_to_cluaster[i] < INT_MAX)
						{
							max = dist_to_cluaster[i];
							max_ind = i;
						}
					}
				}

				if (min_ind != -1)
				{
					first_cluster.push_back(vertexes[min_ind]);
					dist_to_cluaster[min_ind] = (double)INT_MAX + 1;
				}
				if (min_ind != max_ind && max_ind != -1)
				{
					second_cluster.push_back(vertexes[max_ind]);
					dist_to_cluaster[max_ind] = (double)INT_MAX + 1;
				}
			}

			for (int i = 0; i < dist_to_cluaster.size(); ++i)
				if (dist_to_cluaster[i] < INT_MAX)
					second_cluster.push_back(vertexes[i]);


			if (clusters == 2)
				return int_matrix{ first_cluster , second_cluster };
			else if (clusters == 3)
			{
				auto res = dichotomous_division(dist_mat, second_cluster, k2, cluster_size);
				res.push_back(first_cluster);
				return res;
			}
			auto res1 = dichotomous_division(dist_mat, first_cluster, k1, cluster_size);
			auto res2 = dichotomous_division(dist_mat, second_cluster, k2, cluster_size);
			for (const auto& cluster : res2)
				res1.push_back(cluster);
			return res1;
		}

		double demand_of_vertex(const std::vector<double>& weights,
			const std::vector<size_t>& vertexes)
		{
			double sum = 0;
			for (const size_t v : vertexes)
				sum += weights[v];
			return sum;
		}

		int_matrix dichotomous_division_weight(const matrix& dist_mat,
			const std::vector<double>& weights, const double capacity,
			const size_t count_cluster)
		{
			std::vector<size_t> vertexes(dist_mat.size() - 1);
			for (size_t i = 1; i < dist_mat.size(); i++)
				vertexes[i - 1] = i;
			return dichotomous_division_weight(dist_mat, weights, vertexes, capacity, count_cluster, true);
		}

		int_matrix dichotomous_division_weight(const matrix& dist_mat,
			const std::vector<double>& weights,
			const std::vector<size_t>& vertexes,
			const double capacity, const size_t count_cluster, const bool can_reduce)
		{
			double demand = demand_of_vertex(weights, vertexes);
			size_t k1 = count_cluster / 2;
			size_t k2 = count_cluster - k1;
			double D1 = k1 * capacity;
			double D2 = k2 * capacity;

			double d1 = demand - D2;
			double d2 = demand - D1;
			if (can_reduce && (d1 < 0 || d2 < 0))
				return dichotomous_division_weight(dist_mat, weights,
					vertexes, capacity, count_cluster - 1, can_reduce);

			std::vector<double> dist_to_cluaster(vertexes.size());
			for (size_t i = 0; i < dist_to_cluaster.size(); i++)
				dist_to_cluaster[i] = (double)INT_MAX - 1;
			auto pair = two_farthest_vertex(dist_mat, vertexes);
			int first = pair.first;
			int second = pair.second;

			dist_to_cluaster[first] = (double)INT_MAX + 1;
			dist_to_cluaster[second] = (double)INT_MAX + 1;

			std::vector<size_t> first_cluster = { vertexes[first] };
			std::vector<size_t> second_cluster = { vertexes[second] };
			double demand1 = weights[vertexes[first]];
			double demand2 = weights[vertexes[second]];

			while (/*demand1 / k1 < demand  / count_cluster &&
				demand2 / k2 < demand / count_cluster*/
				first_cluster.size() + second_cluster.size() < vertexes.size())
			{
				int min = INT_MAX;
				int max = INT_MIN;
				int min_ind = -1;
				int max_ind = -1;

				for (int i = 0; i < dist_to_cluaster.size(); ++i)
				{
					// inline body of : update_dist_to_cluaster
					if (dist_to_cluaster[i] < INT_MAX)
					{
						dist_to_cluaster[i] = std::min(dist_to_cluaster[i], 
							dist_mat[vertexes[first]][vertexes[i]]
							- dist_mat[vertexes[second]][vertexes[i]]);
						//

						if (min > dist_to_cluaster[i] &&
							weights[vertexes[i]] + demand1 <= D1)
						{
							min = dist_to_cluaster[i];
							min_ind = i;
						}
						if (max < dist_to_cluaster[i] && dist_to_cluaster[i] < INT_MAX
							&& weights[vertexes[i]] + demand2 <= D2)
						{
							max = dist_to_cluaster[i];
							max_ind = i;
						}
					}
				}

				if (min_ind != -1)
				{
					first_cluster.push_back(vertexes[min_ind]);
					dist_to_cluaster[min_ind] = (double)INT_MAX + 1;
					first = min_ind;
					demand1 += weights[vertexes[first]];
				}
				if (min_ind != max_ind && max_ind != -1)
				{
					second_cluster.push_back(vertexes[max_ind]);
					dist_to_cluaster[max_ind] = (double)INT_MAX + 1;
					second = max_ind;
					demand2 += weights[vertexes[second]];
				}
				if (min_ind == -1 && max_ind == -1)
					break;
			}

			if (first_cluster.size() + second_cluster.size() < vertexes.size())
			{
				size_t ind = 0;
				for (; ind < dist_to_cluaster.size(); ++ind)
					if (dist_to_cluaster[ind] < INT_MAX)
						break;
				if (demand1 > demand2)
				{
					second_cluster.push_back(vertexes[ind]);
					demand2 += weights[vertexes[ind]];
				} else {
					first_cluster.push_back(vertexes[ind]);
					demand1 += weights[vertexes[ind]];
				}
				if (demand2 > D2) {// переброс всегда из 1 во 2
					std::swap(first_cluster, second_cluster);
					std::swap(demand1, demand2);
					std::swap(D1, D2);
					std::swap(k1, k2);
				}

				bool is_ok = false;
				for (int i = first_cluster.size() - 1; i > 0; --i)
				{
					if (demand2 + weights[first_cluster[i]] <= D2 &&
						demand1 - weights[first_cluster[i]] <= D1)
					{
						second_cluster.push_back(first_cluster[i]);
						first_cluster.erase(first_cluster.begin() + i);
						is_ok - true;
						break;
					}
				}

				if (!is_ok)
					return dichotomous_division_weight(dist_mat, weights,
						vertexes, capacity, count_cluster + 1, false);
			}

			if (count_cluster == 2)
				return int_matrix{ first_cluster , second_cluster };
			else if (count_cluster == 3)
			{
				if (k1 == 2)
					swap(first_cluster, second_cluster);
				auto res = dichotomous_division_weight(dist_mat, weights,
					second_cluster, capacity, 2, true);
				res.push_back(first_cluster);
				return res;
			}
			auto res1 = dichotomous_division_weight(dist_mat, weights,
				first_cluster, capacity, k1, true);
			auto res2 = dichotomous_division_weight(dist_mat, weights,
				second_cluster, capacity, k2, true);
			for (const auto& cluster : res2)
				res1.push_back(cluster);
			return res1;
		}


		std::vector<size_t> radian_sort(const double* const x, const double* const y, size_t size)
		{
			--size;
			matrix points(size, std::vector<double>(4));
			for (size_t i = 0; i < size; ++i)
			{
				if (abs(x[i + 1]) < 0.000001 && abs(y[i + 1]) < 0.000001)
				{
					points[i][0] = -5;
					points[i][1] = x[i + 1];
					points[i][2] = y[i + 1];
					points[i][3] = i + 1;
					continue;
				}
				points[i][0] = std::atan2(x[i + 1], y[i + 1]);
				points[i][1] = x[i + 1];
				points[i][2] = y[i + 1];
				points[i][3] = i + 1;
			}

			std::sort(points.begin(), points.end());

			auto res = std::vector<size_t>(size);
			for (size_t i = 0; i < size; ++i)
				res[i] = points[i][3];
			return res;
		}
		
		int_matrix sweeping(const double* const x, const double* const y, const size_t size, size_t count_clusters)
		{
			matrix points(size, std::vector<double>(4));
			for (size_t i = 0; i < size; ++i)
			{
				points[i][0] = std::atan2(x[i + 1], y[i + 1]);
				points[i][1] = x[i + 1];
				points[i][2] = y[i + 1];
				points[i][3] = i + 1;
			}

			std::sort(points.begin(), points.end());

			int_matrix clusters;

			while (points.size() > 0)
			{
				std::vector<size_t> cluster;
				size_t size_cluster = points.size() / count_clusters;
				--count_clusters;
				for (size_t i = 0; i < size_cluster; ++i)
				{
					cluster.push_back((int)points.back()[3]);
					points.pop_back();
				}
				clusters.push_back(cluster);
			}
			return clusters;
		}
	}

	int_matrix cutting_rout(const std::vector<size_t>& rout,
		const matrix& dist_mat, size_t count_rout)
	{
		using namespace utils;
		size_t count_points = rout.size();
		int_matrix routs;
		routs.push_back(rout);
		for (size_t i = count_rout; i > 1; --i)
		{
			size_t size_rout = count_points / i;
			std::vector<size_t> new_rout1;
			for (size_t j = 0; j < size_rout; ++j)
			{
				new_rout1.push_back(routs[0].back());
				routs[0].pop_back();
			}
			if (count_points % i > 0)
			{
				std::vector<size_t> new_rout2(new_rout1);
				new_rout2.push_back(routs[0].back());
				std::vector<size_t> rout0(routs[0]);
				rout0.pop_back();
				auto len1 = length_rout(new_rout1, dist_mat) + length_rout(routs[0], dist_mat);
				auto len2 = length_rout(new_rout2, dist_mat) + length_rout(rout0, dist_mat);

				if (len1 > len2)
				{
					new_rout1 = new_rout2;
					routs[0] = rout0;
					++size_rout;
				}
			}
			count_points -= size_rout;
			routs.push_back(new_rout1);

		}

		for (size_t i = 0; i < routs.size(); ++i)
		{
			TSP::local_opt::TSP_2_opt(routs[i], dist_mat);
			TSP::local_opt::TSP_3_opt(routs[i], dist_mat);
			TSP::local_opt::TSP_2_opt(routs[i], dist_mat);
			TSP::local_opt::TSP_3_opt(routs[i], dist_mat);
		}
		return routs;
	}

}