#include "../headers/balancedVRP.hpp"

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
			std::vector<size_t> vertexes(dist_mat.size());
			for (size_t i = 0; i < dist_mat.size(); i++)
				vertexes[i] = i;
			return dichotomous_division(dist_mat, vertexes, clusters);
		}


		int_matrix dichotomous_division(const matrix& dist_mat,
			const std::vector<size_t>& vertexes, const int clusters)
		{
			int cluster_size = vertexes.size() / clusters;
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

			while (vertexes.size() - first_cluster.size() > k2 * (cluster_size + 1))
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
				auto res = dichotomous_division(dist_mat, second_cluster, k2);
				res.push_back(first_cluster);
				return res;
			}
			auto res1 = dichotomous_division(dist_mat, first_cluster, k1);
			auto res2 = dichotomous_division(dist_mat, second_cluster, k2);
			for (const auto& cluster : res2)
				res1.push_back(cluster);
			return res1;
		}
	}

}