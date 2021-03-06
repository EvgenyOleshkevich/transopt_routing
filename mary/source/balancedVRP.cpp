#include "../headers/balancedVRP.hpp"

#include <math.h>

namespace balancedVRP
{
	namespace clustering
	{
		pair<int, int> two_farthest_vertex(const matrix& dist_mat,
			const vector<size_t>& vertexes) {
			double max_dist = 0;
			pair<int, int> res = { 0, 1 };
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
			const vector<size_t>& vertexes, vector<double>& dist_to_cluaster,
			const int first, const int second)
		{
			for (int i = 0; i < vertexes.size(); ++i)
				if (dist_to_cluaster[i] < INT_MAX)
					dist_to_cluaster[i] = dist_mat[vertexes[first]][vertexes[i]]
					- dist_mat[vertexes[second]][vertexes[i]];
		}

		int_matrix dichotomous_division(const matrix& dist_mat, const int clusters)
		{
			vector<size_t> vertexes(dist_mat.size() - 1);
			for (size_t i = 1; i < dist_mat.size(); i++)
				vertexes[i - 1] = i;
			return dichotomous_division(dist_mat, vertexes, clusters, vertexes.size() / clusters);
		}


		int_matrix dichotomous_division(const matrix& dist_mat,
			const vector<size_t>& vertexes, const int clusters, const int cluster_size)
		{
			int k1 = clusters / 2;
			int k2 = clusters / 2 + clusters % 2;

			vector<double> dist_to_cluaster(vertexes.size());
			for (size_t i = 0; i < dist_to_cluaster.size(); i++)
				dist_to_cluaster[i] = (double)INT_MAX - 1;
			auto pair = two_farthest_vertex(dist_mat, vertexes);
			int first = pair.first;
			int second = pair.second;

			dist_to_cluaster[first] = (double)INT_MAX + 1;
			dist_to_cluaster[second] = (double)INT_MAX + 1;

			vector<size_t> first_cluster = { vertexes[first] };
			vector<size_t> second_cluster = { vertexes[second] };

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

		double demand_of_vertex(const vector<double>& weights,
			const vector<size_t>& vertexes)
		{
			double sum = 0;
			for (const size_t v : vertexes)
				sum += weights[v];
			return sum;
		}

		int_matrix dichotomous_division_weight(const matrix& dist_mat,
			const vector<double>& weights, const double capacity,
			const size_t count_cluster)
		{
			vector<size_t> vertexes(dist_mat.size() - 1);
			for (size_t i = 1; i < dist_mat.size(); i++)
				vertexes[i - 1] = i;
			auto res = dichotomous_division_weight(dist_mat, weights, vertexes, capacity, count_cluster, true);
			for (size_t i = 0; i < res.size(); i++)
				res[i].emplace(res[i].begin(), 0);
			return res;
		}

		int_matrix dichotomous_division_weight(const matrix& dist_mat,
			const vector<double>& weights,
			const vector<size_t>& vertexes,
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

			vector<double> dist_to_cluaster(vertexes.size());
			for (size_t i = 0; i < dist_to_cluaster.size(); i++)
				dist_to_cluaster[i] = (double)INT_MAX - 1;
			auto pair = two_farthest_vertex(dist_mat, vertexes);
			int first = pair.first;
			int second = pair.second;

			dist_to_cluaster[first] = (double)INT_MAX + 1;
			dist_to_cluaster[second] = (double)INT_MAX + 1;

			vector<size_t> first_cluster = { vertexes[first] };
			vector<size_t> second_cluster = { vertexes[second] };
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
				if (demand2 > D2) {// �������� ������ �� 1 �� 2
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


		vector<size_t> radian_sort(const vector<double>& x, const vector<double>& y)
		{
			if (x.size() != y.size())
				throw std::runtime_error("size dont match");

			size_t size = x.size();
			vector<pair<double, size_t>> points(size);
			for (size_t i = 0; i < size; ++i)
			{
				if (abs(x[i]) < 0.000001 && abs(y[i]) < 0.000001)
				{
					points[i].first = -5;
					points[i].second = i;
					continue;
				}
				points[i].first = std::atan2(x[i + 1], y[i + 1]);
				points[i].second = i;
			}

			std::sort(points.begin(), points.end());

			auto res = vector<size_t>(size);
			for (size_t i = 0; i < size; ++i)
				res[i] = points[i].second;
			return res;
		}
		
		int_matrix sweeping(const double* const x, const double* const y, const size_t size, size_t count_clusters)
		{
			matrix points(size, vector<double>(4));
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
				vector<size_t> cluster;
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
	
		vector<matrix> get_dist_inner_cluster(const matrix& dist_mat, const int_matrix& clusters)
		{
			vector<matrix> res;
			for (size_t i = 0; i < clusters.size(); ++i) {
				matrix cluster_dist_mat =
					matrix(clusters[i].size(),
						vector<double>(clusters[i].size()));

				for (size_t j = 0; j < clusters[i].size(); ++j)
					for (size_t k = 0; k < clusters[i].size(); ++k)
						cluster_dist_mat[j][k] = dist_mat[clusters[i][j]][clusters[i][k]];

				res.push_back(cluster_dist_mat);
			}
			return res;
		}

		matrix get_weight_inner_cluster(const vector<double>& weight, const int_matrix& clusters)
		{
			matrix weights;
			for (size_t i = 0; i < clusters.size(); ++i) {
				auto cluster_weight = vector<double>(clusters[i].size());

				for (size_t j = 0; j < clusters[i].size(); ++j)
					cluster_weight[j] = weight[clusters[i][j]];

				weights.push_back(cluster_weight);
			}
			return weights;
		}

		vector<int> get_number_cluster_by_vertex(const int_matrix& clusters)
		{
			size_t size = 1;
			for (const auto& cluster : clusters)
				size += cluster.size();

			vector<int> vertex_x_cluster(size, -1);

			for (size_t i = 0; i < clusters.size(); ++i)
				for (const auto vertex : clusters[i])
					//if (vertex_x_cluster[vertex] == -1)
					vertex_x_cluster[vertex] = i;
					/*else
					{
						int y = 0;
						int e = 0;
					}

			for (size_t i = 1; i < vertex_x_cluster.size(); ++i)
				if (vertex_x_cluster[i] < 0 || vertex_x_cluster[i] > 6)
				{
					int y = 0;
					int e = 0;
				}*/

			return vertex_x_cluster;
		}
	}

	namespace project
	{
		int_matrix clusters_with_multiple_duplicates(const int_matrix& clusters, const vector<size_t>& frequence)
		{
			int_matrix res = clusters;

			for (size_t i = 0; i < clusters.size(); ++i)
				for (size_t j = 0; j < clusters[i].size(); ++j)
					if (frequence[clusters[i][j]] == 4) // clusters.size() - max period
					{
						res[(i + 2) % clusters.size()].push_back(clusters[i][j]);
						res[(i + 4) % clusters.size()].push_back(clusters[i][j]);
						res[(i + 6) % clusters.size()].push_back(clusters[i][j]);
					}
					else if (frequence[clusters[i][j]] == 3) {
						res[(i + 3) % clusters.size()].push_back(clusters[i][j]);
						res[(i + 6) % clusters.size()].push_back(clusters[i][j]);
					}
					else if (frequence[clusters[i][j]] == 2) {
						res[(i + 4) % clusters.size()].push_back(clusters[i][j]);
					}

			return res;
		}

		vector<size_t> radian_sort(const vector<double>& x, const vector<double>& y)
		{
			if (x.size() != y.size())
				throw std::runtime_error("size dont match");

			size_t size = x.size() - 1;
			vector<pair<double, size_t>> points(size);

			double cx = x[0];
			double cy = y[0];
			for (size_t i = 0; i < size; ++i)
			{
				if (abs(x[i + 1] - cx) < 0.000001 && abs(y[i + 1] - cy) < 0.000001)
				{
					points[i].first = -5;
					points[i].second = i + 1;
					continue;
				}
				points[i].first = std::atan2(x[i + 1] - cx, y[i + 1] - cy);
				points[i].second = i + 1;
			}

			std::sort(points.begin(), points.end());

			auto res = vector<size_t>(size);
			for (size_t i = 0; i < size; ++i)
				res[i] = points[i].second;
			return res;
		}
	}

	int_matrix cutting_rout(const vector<size_t>& rout,
		const matrix& dist_mat, size_t count_rout)
	{
		using namespace utils;
		size_t count_points = rout.size();
		int_matrix routs;
		routs.push_back(rout);
		for (size_t i = count_rout; i > 1; --i)
		{
			size_t size_rout = count_points / i;
			vector<size_t> new_rout1;
			for (size_t j = 0; j < size_rout; ++j)
			{
				new_rout1.push_back(routs[0].back());
				routs[0].pop_back();
			}
			if (count_points % i > 0)
			{
				vector<size_t> new_rout2(new_rout1);
				new_rout2.push_back(routs[0].back());
				vector<size_t> rout0(routs[0]);
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
			TSP::local_opt::opt_2(routs[i], dist_mat);
			TSP::local_opt::opt_3(routs[i], dist_mat);
			TSP::local_opt::opt_2(routs[i], dist_mat);
			TSP::local_opt::opt_3(routs[i], dist_mat);
		}
		return routs;
	}

#pragma region VND_STS
	VND_STS::VND_STS(const matrix& dist_mat, const size_t need_routs
		, const double capacity, const vector<double>& weight) :
		dist_mat(dist_mat), count_point(dist_mat.size()), need_routs(need_routs),
		capacity(capacity), weight(weight), fict_vertex(count_point)
	{
		for (auto& vec : this->dist_mat)
		{
			vec.push_back(0);
		}
		this->dist_mat.push_back(vector<double>(fict_vertex + 1, 0));
	}

	pair<int_matrix, matrix> VND_STS::dynamic_decode(const vec_int_float& vertexes)
	{
		vector < size_t> init_rout = { vertexes[0].first };
		double length = utils::length_rout_fict(init_rout, dist_mat, fict_vertex);
		vector<pair<double, int>> F = { {length , -1} }; // F(0);
		// length of routs, index for restore routs

		for (size_t i = 1; i < vertexes.size(); ++i)
		{
			double weight = 0;
			vector <size_t> rout;
			F.push_back({ INT_MAX, i - 1 });
			for (int j = i; j > 0; --j)
			{
				rout.push_back(vertexes[j].first);
				weight += vertexes[j].second;
				if (weight > capacity)
					break;
				length = utils::length_rout_fict(rout, dist_mat, fict_vertex);
				if (j > 0)
					length += F[j - 1].first;

				if (F[i].first > length)
				{
					F[i].first = length;
					F[i].second = j - 1;
				}
			}
		}

		int_matrix routs;
		matrix routs_weight;
		int i = F.size() - 1;
		do
		{
			vector <size_t> rout;
			vector <double> rout_weight;
			for (int j = i; j > F[i].second; --j)
			{
				rout.push_back(vertexes[j].first);
				rout_weight.push_back(vertexes[j].second);
			}

			routs.push_back(rout);
			routs_weight.push_back(rout_weight);
			i = F[i].second;
		} while (i > 0);

		return { routs, routs_weight };
	}

	pair<int_matrix, matrix> VND_STS::greedy_decode(const vec_int_float& vertexes)
	{
		int_matrix routs;
		matrix routs_weight;

		vector<size_t> rout;
		vector<double> rout_weight;
		double w = 0;
		for (size_t i = 0; i < vertexes.size(); ++i)
		{
			rout.push_back(vertexes[i].first);
			w += vertexes[i].second;

			if (w > capacity)
			{
				w -= vertexes[i].second + capacity;
				rout_weight.push_back(-w);
				--i;
				routs.push_back(rout);
				routs_weight.push_back(rout_weight);
				rout = vector<size_t>();
				rout_weight = vector<double>();
			}
			else {
				rout_weight.push_back(vertexes[i].second);
			}
		}

		routs.push_back(rout);
		routs_weight.push_back(rout_weight);

		return { routs , routs_weight };
	}

	double VND_STS::dynamic_decode_fast(const vec_int_float& vertexes)
	{
		vector < size_t> init_rout = { vertexes[0].first };
		double length = utils::length_rout_fict(init_rout, dist_mat, fict_vertex);
		vector<pair<double, int>> F = { {length , -1} }; // F(0);
		// length of routs, index for restore routs

		for (size_t i = 1; i < vertexes.size(); ++i)
		{
			double weight = 0;
			vector <size_t> rout;
			F.push_back({ INT_MAX, i - 1 });
			for (int j = i; j >= 0; --j)
			{
				rout.push_back(vertexes[j].first);
				weight += vertexes[j].second;
				if (weight > capacity)
					break;
				length = utils::length_rout_fict(rout, dist_mat, fict_vertex);
				if (j > 0)
					length += F[j - 1].first;

				if (F[i].first > length)
				{
					F[i].first = length;
					F[i].second = j - 1;
				}
			}
		}

		return F.back().first;
	}

	double VND_STS::greedy_decode_fast(const vec_int_float& vertexes)
	{
		double length = 0;
		vector<size_t> rout;
		double w = 0;
		for (size_t i = 0; i < vertexes.size(); ++i)
		{
			rout.push_back(vertexes[i].first);
			w += vertexes[i].second;

			if (w > capacity)
			{
				w -= vertexes[i].second + capacity;
				--i;
				length += utils::length_rout_fict(rout, dist_mat, fict_vertex);
				rout = vector<size_t>();
			}
		}
		length += utils::length_rout_fict(rout, dist_mat, fict_vertex);
		return length;
	}

	void VND_STS::dynamic_decode_add_fict(vec_int_float& vertexes)
	{
		vector < size_t> init_rout = { vertexes[0].first };
		double length = utils::length_rout_fict(init_rout, dist_mat, fict_vertex);
		vector<pair<double, int>> F = { {length , -1} }; // F(0);
		// length of routs, index for restore routs

		for (size_t i = 1; i < vertexes.size(); ++i)
		{
			double weight = 0;
			vector <size_t> rout;
			F.push_back({ INT_MAX, i - 1 });
			for (int j = i; j > 0; --j)
			{
				rout.push_back(vertexes[j].first);
				weight += vertexes[j].second;
				if (weight > capacity)
					break;
				length = utils::length_rout_fict(rout, dist_mat, fict_vertex);
				if (j > 0)
					length += F[j - 1].first;

				if (F[i].first > length)
				{
					F[i].first = length;
					F[i].second = j - 1;
				}
			}
		}

		int_matrix routs;
		int i = F.size() - 1;
		do
		{
			double w = 0;
			int j = i;
			for (; j > F[i].second; --j)
				w += vertexes[j].second;
			if (w < capacity && w > 3 * capacity / 4)
			{
				++j;
				vertexes.emplace(vertexes.begin() + j,
					pair<size_t, double>(fict_vertex, capacity - w));
			}
			i = F[i].second;
		} while (i > 0);
	}

	double VND_STS::VND(vec_int_float& vertexes)
	{
		vec_int_float best_vertexes = vertexes;
		double best_length = dynamic_decode_fast(vertexes);
		while (true)
		{
			TSP::local_opt_for_VND_STS::TSP_2_opt(vertexes, dist_mat);
			double length = dynamic_decode_fast(vertexes);
			bool is_break = true;
			if (best_length > length)
			{
				best_length = length;
				best_vertexes = vertexes;
				is_break = false;
			}
			TSP::local_opt_for_VND_STS::swap(vertexes, dist_mat);
			length = dynamic_decode_fast(vertexes);
			if (best_length > length)
			{
				best_length = length;
				best_vertexes = vertexes;
				is_break = false;
			}
			TSP::local_opt_for_VND_STS::shift(vertexes, dist_mat);
			length = dynamic_decode_fast(vertexes);
			if (best_length > length)
			{
				best_length = length;
				best_vertexes = vertexes;
				is_break = false;
			}
			if (is_break)
			{
				break;
			}
		}

		vertexes = best_vertexes;
		return best_length;
	}
	
	double VND_STS::STS(vec_int_float& vertexes)
	{
		vec_int_float best_vertexes = vertexes;
		double best_length = greedy_decode_fast(vertexes);

		for (size_t i = 0; i < 30; ++i)
		{
			size_t a = rand() % vertexes.size();
			size_t b = rand() % vertexes.size();
			while (a == b)
				b = rand() % vertexes.size();
			if (a > b)
				std::swap(a, b);

			size_t type = rand() % 4;
			switch (type)
			{
			case 0:
			{
				TSP::local_opt_for_VND_STS::TSP_2_opt(vertexes, dist_mat);
				double length = greedy_decode_fast(vertexes);
				if (best_length > length)
				{
					best_length = length;
					best_vertexes = vertexes;
				}
				break;
			}
			case 1:
			{
				std::swap(vertexes[a], vertexes[b]);
				double length = greedy_decode_fast(vertexes);
				if (best_length > length)
				{
					best_length = length;
					best_vertexes = vertexes;
				}
				break;
			}
			case 2:
			{
				for (size_t j = a; j < b - 1; ++j)
					std::swap(vertexes[j], vertexes[j + 1]);
				double length = greedy_decode_fast(vertexes);
				if (best_length > length)
				{
					best_length = length;
					best_vertexes = vertexes;
				}
				break;
			}
			case 3:
			{
				TSP::local_opt_for_VND_STS::exchange(vertexes, a, b);
				double length = greedy_decode_fast(vertexes);
				if (best_length > length)
				{
					best_length = length;
					best_vertexes = vertexes;
				}
				break;
			}
			default:
				break;
			}
		}
		vertexes = best_vertexes;
		return best_length;
	}

	pair<int_matrix, matrix> VND_STS::calculate(vec_int_float& vertexes)
	{
		vec_int_float best_vertexes = vertexes;
		double best_length = dynamic_decode_fast(vertexes);

		for (size_t i = 0; i < 10; ++i)
		{
			double length = VND(vertexes);
			if (best_length > length)
			{
				best_length = length;
				best_vertexes = vertexes;
			}
			//dynamic_decode_add_fict(vertexes);
			length = STS(vertexes);
			if (best_length > length)
			{
				best_length = length;
				best_vertexes = vertexes;
			}

			for (size_t j = 1; j < vertexes.size(); ++j)
				if (vertexes[j].first == vertexes[j - 1].first)
				{
					vertexes[j].second += vertexes[j - 1].second;
					--j;
					vertexes.erase(vertexes.begin() + j);
				}
		}
		vertexes = best_vertexes;
		/*for (size_t i = 0; i < vertexes.size(); ++i)
			if (vertexes[i].first == fict_vertex)
			{
				vertexes.erase(vertexes.begin() + i);
				--i;
			}*/

		return dynamic_decode(vertexes);
	}
#pragma endregion

#pragma region Ant
	Ant_algorithm::Ant_algorithm(const matrix& dist_mat, const vector<double>& weights, const vector<Transport>& transports) :
		dist_mat(dist_mat), transports(transports), weights(weights) 
	{
		std::random_device rd;
		mersenne_rand = std::mt19937(rd());
	}

#pragma endregion

}