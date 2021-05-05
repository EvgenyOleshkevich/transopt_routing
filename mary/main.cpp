#include "headers/algorithms.hpp"
#include <fstream>
#include <thread>
#include <windows.h>




using namespace std;
using namespace algorithms;
using namespace balancedVRP;

vector<double> x;
vector<double> y;
vector<double> weights;
vector<size_t> frequence;
vector<vector<double>> dist_mat;
vector<Transport> transports;

void read_transports()
{
    std::ifstream in("data_transport.csv");
    string line;
    getline(in, line);
    while (!in.eof())
    {
        double volume, cost_start, cost_by_dist;
        size_t count;
        in >> volume >> cost_start >> cost_by_dist >> count;
        transports.push_back({ volume, cost_start, cost_by_dist, count });
    }
    in.close();

    sort(transports.begin(), transports.end(), [](const Transport& a, const Transport& b)
        {
            return a.capacity < b.capacity;
            //return a.second.sep.size() < b.second.sep.size();
        });
}

void read_file()
{
    x.clear();
    y.clear();
    weights.clear();

    std::ifstream in("data_vertex.csv");
    string line;
    getline(in, line);
    double s = 0;
    //double X_mean = 0;
    //double Y_mean = 0;
    
    while (!in.eof())
    {
        double freq, w, X, Y;
        in >> freq >> w >> X >> Y;
        frequence.push_back(freq);
        x.push_back(X);
        y.push_back(Y);
        //X_mean += X;
        //Y_mean += Y;
        weights.push_back(w);
        s += w;
    }
    //x[0] = X_mean / (x.size() - 1);
    //y[0] = Y_mean / (y.size() - 1);
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

void Test(vector<double>& vec) {
    vec.clear();
}

void print(const vector<int_matrix>& routs)
{
    for (const int_matrix& rout_mat : routs)
        for (const vector<size_t>& rout : rout_mat)
        {
            cout << "[";
            for (const size_t vertex : rout)
                cout << vertex << ", ";
            cout << 0 << "]," << endl;
        }
}

void print(const int_matrix& routs)
{
    for (const vector<size_t>& rout : routs)
    {
        cout << "[";
        for (const size_t vertex : rout)
            cout << vertex << ", ";
        cout << 0 << "]," << endl;
    }
}

bool check_matrix(const vector<int_matrix>& routs, const size_t size)
{
    vector<size_t> used(size, 0);
    for (const int_matrix& rout_mat : routs)
    {
        //cout << "new type" << endl;
        for (const vector<size_t>& rout : rout_mat)
        {
            //cout << "[" ;
            for (const size_t vertex : rout)
            {
                ++used[vertex];
                //cout << vertex << ", ";
            }
            //cout << 0 << "]," << endl;
        }
    }
    bool is_twice = false;
    bool is_null = false;
    for (size_t i = 1; i < used.size(); ++i)
    {
        if (used[i] > 1)
            is_twice = true;
        if (used[i] == 0)
            is_null = true;
    }
    cout << "check is_twice: " << is_twice << endl;
    cout << "check is_null: " << is_null << endl;
    return !is_twice && !is_null;
}

bool check_matrix(const int_matrix& routs, const size_t size)
{
    vector<size_t> used(size, 0);

    for (const vector<size_t>& rout : routs)
        for (const size_t vertex : rout)
            ++used[vertex];
    bool is_twice = false;
    bool is_null = false;
    for (size_t i = 1; i < used.size(); ++i)
    {
        if (used[i] > 1)
            is_twice = true;
        if (used[i] == 0)
            is_null = true;
    }
    cout << "check is_twice: " << is_twice << endl;
    cout << "check is_null: " << is_null << endl;
    return !is_twice && !is_null;
}

void Ant_test()
{
    read_file();
    read_transports();
    dist_mat = utils::fill_matrix_with_end_point(x, y);

    balancedVRP::Ant_algorithm ant(dist_mat, transports, weights);
    auto mat_len = ant.run();
    cout << "lenght: " << mat_len.second << endl;
    cout << "check: " << check_matrix(mat_len.first, dist_mat.size()) << endl;
}

void sweep_test()
{
    read_file();
    read_transports();
    dist_mat = utils::fill_matrix_with_end_point(x, y);

    balancedVRP::project::Sweeping sweep(dist_mat, x, y, weights, transports);
    sweep.run();
    cout << "lenght: " << sweep.lenght() << endl;
    print(sweep.res);
    cout << "check: " << check_matrix(sweep.res, dist_mat.size()) << endl;
}

void greedy_test()
{
    read_file();
    read_transports();
    auto P = utils::fill_matrix_and_sort(x, y);
    dist_mat = P.first;

    balancedVRP::project::GreadyBase greedy(dist_mat, P.second, weights, transports);
    greedy.run();
    cout << "lenght: " << greedy.lenght() << endl;
    print(greedy.res);
    cout << "check: " << check_matrix(greedy.res, dist_mat.size()) << endl;
}

void greedy_2opt_test()
{
    read_file();
    read_transports();
    auto P = utils::fill_matrix_and_sort(x, y);
    dist_mat = P.first;

    balancedVRP::project::GreadyBase greedy(dist_mat, P.second, weights, transports);
    greedy.run();
    for (int_matrix& rout_mat : greedy.res)
        for (vector<size_t>& rout : rout_mat)
        {
            TSP::local_opt::opt_2_fast2(rout, dist_mat);
            //TSP::local_opt::opt_3(rout, dist_mat);
        }
    cout << "lenght: " << greedy.lenght() << endl;
    print(greedy.res);
    cout << "check: " << check_matrix(greedy.res, dist_mat.size()) << endl;
}

void clurk_test()
{
    read_file();
    read_transports();
    auto P = utils::fill_matrix_and_sort(x, y);
    dist_mat = P.first;

    balancedVRP::project::ClarkRight clurk(dist_mat, weights, transports);
    clurk.run();

    cout << "lenght: " << clurk.lenght() << endl;
    print(clurk.res);
    cout << "check: " << check_matrix(clurk.res, dist_mat.size()) << endl;
}

void osman_test()
{
    read_file();
    read_transports();
    auto P = utils::fill_matrix_and_sort(x, y);
    dist_mat = P.first;

    balancedVRP::project::GreadyBase greedy(dist_mat, P.second, weights, transports);
    greedy.run();
    cout << "lenght greedy: " << greedy.lenght() << endl;
    cout << "check: " << check_matrix(greedy.res, dist_mat.size()) << endl;

    auto data = balancedVRP::project::Osman::transform_data(greedy.res);

    balancedVRP::project::Osman osman(dist_mat, weights, transports, data.first,
        data.second, data.first.size());

    osman.run();
    osman.add_zero_vertex();

    print(osman.res);
    cout << "lenght osman: " << osman.lenght() << endl;
    cout << "check: " << check_matrix(osman.res, dist_mat.size()) << endl;
}

int main()
{
    osman_test();
    return 0;
    read_file();
    read_transports();
    dist_mat = utils::fill_matrix_with_end_point(x, y);
    Test(dist_mat[0]);
    int y = 0;
    return 0;
    read_transports();
    //srand((unsigned int)time(0));
    srand(0);
    rand();

    std::random_device rd;
    std::mt19937 mersenne(rd());

    for (int count = 0; count < 48; ++count)
    {
        std::cout << mersenne() << "\t";

        // Если вывели 5 чисел, то вставляем символ новой строки
        if ((count + 1) % 5 == 0)
            std::cout << "\n";
    }
    //read_file();
    //cluster_test2();

    return 0;
}

