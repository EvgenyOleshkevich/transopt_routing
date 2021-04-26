#include "headers/algorithms.hpp"
#include <fstream>
#include <thread>
#include <windows.h>
#include <chrono>

using namespace std;
using namespace algorithms;

vector<double> x;
vector<double> y;
vector<double> weights;
vector<size_t> frequence;
vector<vector<double>> dist_mat;


void read_file()
{
    x.clear();
    y.clear();
    weights.clear();

    x.push_back(0);
    y.push_back(0);
    weights.push_back(0);
    frequence.push_back(0);
    std::ifstream in("cpp_data.csv");
    string line;
    getline(in, line);
    double s = 0;
    while (!in.eof())
    {
        double freq, w, X, Y;
        in >> freq >> w >> X >> Y;
        frequence.push_back(freq);
        x.push_back(X);
        y.push_back(Y);
        weights.push_back(w);
        s += w;
    }
    cout << "sum volume: " << s << endl;
    in.close();
}

void cluster_test1()
{
    dist_mat = utils::fill_matrix_with_end_point(x, y);
    double capacity = 700;
    auto clusters = balancedVRP::clustering::dichotomous_division_weight(dist_mat, weights, capacity, 7);

    double acc = 0;
    double acc_weight = 0;
    for (const auto& cluster : clusters)
        for (auto vertex : cluster)
        {
            if (frequence[vertex] == 2)
            {
                acc += 3;
                acc_weight += 3 * weights[vertex];
            }
            else if (frequence[vertex] == 3)
            {
                acc += 2;
                acc_weight += 2 * weights[vertex];
            }
            else if (frequence[vertex] == 4)
            {
                acc += 1;
                acc_weight += weights[vertex];
            }
        }

    cout << "frequence: " << acc << endl;
    cout << "acc_weight: " << acc_weight << endl;


    auto vertex_x_clust = balancedVRP::clustering::get_number_cluster_by_vertex(clusters);

    cout << "vertex_x_clust.size(): " << vertex_x_clust.size() << endl;
    //    sum volume : 4215.56
    //    frequence : 396
    //    acc_weight : 1225.62
}

void cluster_test2()
{
    dist_mat = utils::fill_matrix_with_end_point(x, y);
    double capacity = 700;
    auto clusters = balancedVRP::clustering::dichotomous_division_weight(dist_mat, weights, capacity, 7);

    double acc = 0;
    double acc_weight = 0;
    for (const auto& cluster : clusters)
        for (auto vertex : cluster)
        {
            if (frequence[vertex] == 2)
            {
                acc += 3;
                acc_weight += 3 * weights[vertex];
            }
            else if (frequence[vertex] == 3)
            {
                acc += 2;
                acc_weight += 2 * weights[vertex];
            }
            else if (frequence[vertex] == 4)
            {
                acc += 1;
                acc_weight += weights[vertex];
            }
        }

    cout << "frequence: " << acc << endl;
    cout << "acc_weight: " << acc_weight << endl;


    auto clusters2 = balancedVRP::project::clusters_with_multiple_duplicates(clusters, frequence);

    acc_weight = 0;
    for (const auto& cluster : clusters2)
        for (auto vertex : cluster)
        {
            acc_weight += weights[vertex];
        }

    cout << "acc_weight: " << acc_weight << endl;
    //    sum volume : 4215.56
    //    frequence : 396
    //    acc_weight : 1225.62
    //    5441.18
}

int main()
{
    read_file();
    cluster_test2();

    return 0;
}

