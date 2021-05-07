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
const double clust_capacity = 700;

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

void print(const vector<int_matrix>& routs, const int_matrix& clusters, const size_t clust_id)
{
    cout << "[" << endl;
    for (const int_matrix& rout_mat : routs)
        for (const vector<size_t>& rout : rout_mat)
        {
            cout << "[";
            for (const size_t vertex : rout)
                cout << clusters[clust_id][vertex] << ", ";
            cout << 0 << "]," << endl;
        }
    cout << "]," << endl;
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

void print_file(const vector<int_matrix>& routs)
{
    std::ofstream out("res.txt");
    for (const int_matrix& rout_mat : routs)
        for (const vector<size_t>& rout : rout_mat)
        {
            for (const size_t vertex : rout)
                out << vertex << " ";
            out << ";";
        }
    out.close();
}

void print_file(const int_matrix& routs)
{
    std::ofstream out("res.txt");
    for (const vector<size_t>& rout : routs)
    {
        for (const size_t vertex : rout)
            out << vertex << " ";
        out << ";";
    }
    out.close();
}

void print_file(const vector<int_matrix>& routs, const int_matrix& clusters, const size_t clust_id)
{
    std::ofstream out("res_clust.txt", std::ios::app);
    for (const int_matrix& rout_mat : routs)
        for (const vector<size_t>& rout : rout_mat)
        {
            for (const size_t vertex : rout)
                out << clusters[clust_id][vertex] << " ";
            out << ";";
        }
    out << ":";
    out.close();
}

void print_file(const int_matrix& routs, const int_matrix& clusters, const size_t clust_id)
{
    std::ofstream out("res_clust.txt", std::ios::app);
    for (const vector<size_t>& rout : routs)
    {
        for (const size_t vertex : rout)
            out << clusters[clust_id][vertex] << " ";
        out << ";";
    }
    out << ":";
    out.close();
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
    //cout << "check is_twice: " << is_twice << endl;
    //cout << "check is_null: " << is_null << endl;
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

    Ant_algorithm ant(dist_mat, transports, weights);
    auto mat_len = ant.run();
    cout << "lenght: " << mat_len.second << endl;
    cout << "check: " << check_matrix(mat_len.first, dist_mat.size()) << endl;
}

void sweep_test()
{
    read_file();
    read_transports();
    dist_mat = utils::fill_matrix_with_end_point(x, y);

    project::Sweeping sweep(dist_mat, x, y, weights, transports);
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

    project::GreadyBase greedy(dist_mat, P.second, weights, transports);
    greedy.run();
    cout << "lenght: " << greedy.lenght() << endl;
    print_file(greedy.res);
    cout << "check: " << check_matrix(greedy.res, dist_mat.size()) << endl;
}

void greedy_2opt_test()
{
    read_file();
    read_transports();
    auto P = utils::fill_matrix_and_sort(x, y);
    dist_mat = P.first;

    project::GreadyBase greedy(dist_mat, P.second, weights, transports);
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

    project::ClarkRight clurk(dist_mat, weights, transports);
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

    project::GreadyBase greedy(dist_mat, P.second, weights, transports);
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

void save_cluster()
{
    read_file();
    read_transports();
    dist_mat = utils::fill_matrix(x, y);
    auto clusters = balancedVRP::clustering::dichotomous_division_weight(dist_mat, weights, clust_capacity, 7);

    std::ofstream out("cluster.txt");
    for (size_t i = 0; i < clusters.size(); i++)
    {
        double weight = 0;
        for (size_t j = 0; j < clusters[i].size(); j++)
        {
            weight += weights[clusters[i][j]];
            out << clusters[i][j] << ' ';
        }
        out << ";";
        cout << i << " claut weight: " << weight << " size: " << clusters[i].size() << endl;
    }
    out.close();

    //    sum volume : 4215.56
    //    frequence : 396
    //    acc_weight : 1225.62
    //    5441.18
}

void greedy_cluster_test()
{
    std::ofstream out("res_clust.txt");
    out.close();
    read_file();
    read_transports();
    dist_mat = utils::fill_matrix(x, y);
    auto clusters = balancedVRP::clustering::dichotomous_division_weight(dist_mat, weights, clust_capacity, 7);
    
    auto dist_inner_cluster = balancedVRP::clustering::get_dist_inner_cluster(dist_mat, clusters);
    auto weight_inner_cluster = balancedVRP::clustering::get_weight_inner_cluster(weights, clusters);


    double lenght = 0;
    vector<size_t> checks;
    for (size_t i = 0; i < clusters.size(); i++)
    {
        auto sort_matrix = utils::fill_sort_matrix(dist_inner_cluster[i]);

        balancedVRP::project::GreadyBase greedy(
            dist_inner_cluster[i],
            sort_matrix,
            weight_inner_cluster[i],
            transports);
        greedy.run();
        lenght += greedy.lenght();
        print_file(greedy.res, clusters, i);
        checks.push_back(check_matrix(greedy.res, dist_inner_cluster[i].size()));
    }

    cout << "lenght: " << lenght << endl;
    for (size_t check : checks)
        cout << "check: " << check << endl;
    
}

void osman_cluster_test()
{
    std::ofstream out("res_clust.txt");
    out.close();
    read_file();
    read_transports();
    dist_mat = utils::fill_matrix(x, y);
    auto clusters = balancedVRP::clustering::dichotomous_division_weight(dist_mat, weights, clust_capacity, 7);

    auto dist_inner_cluster = balancedVRP::clustering::get_dist_inner_cluster(dist_mat, clusters);
    auto weight_inner_cluster = balancedVRP::clustering::get_weight_inner_cluster(weights, clusters);


    double lenght = 0;
    vector<size_t> checks;
    for (size_t i = 0; i < clusters.size(); i++)
    {
        auto sort_matrix = utils::fill_sort_matrix(dist_inner_cluster[i]);

        balancedVRP::project::GreadyBase greedy(
            dist_inner_cluster[i],
            sort_matrix,
            weight_inner_cluster[i],
            transports);
        greedy.run();

        auto data = balancedVRP::project::Osman::transform_data(greedy.res);

        balancedVRP::project::Osman osman(dist_inner_cluster[i],
            weight_inner_cluster[i],
            transports, data.first,
            data.second, data.first.size());

        osman.run();
        osman.add_zero_vertex();

        lenght += osman.lenght();


        print_file(osman.res, clusters, i);
        checks.push_back(check_matrix(osman.res, dist_inner_cluster[i].size()));
    }

    cout << "lenght: " << lenght << endl;
    for (size_t check : checks)
        cout << "check: " << check << endl;

}

void greedy_cluster_test2()
{
    std::ofstream out("res_clust.txt");
    out.close();
    read_file();
    read_transports();
    dist_mat = utils::fill_matrix(x, y);
    auto clusters = balancedVRP::clustering::dichotomous_division_weight(dist_mat, weights, clust_capacity, 7);

    auto dist_inner_cluster = balancedVRP::clustering::get_dist_inner_cluster(dist_mat, clusters);
    auto weight_inner_cluster = balancedVRP::clustering::get_weight_inner_cluster(weights, clusters);


    double lenght = 0;
    vector<size_t> checks;
    for (size_t i = 0; i < clusters.size(); i++)
    {
        auto sort_matrix = utils::fill_sort_matrix(dist_inner_cluster[i]);

        balancedVRP::project::GreadyBase greedy(
            dist_inner_cluster[i],
            sort_matrix,
            weight_inner_cluster[i],
            transports);
        greedy.run();
        lenght += greedy.lenght();

        for (size_t i = 0; i < greedy.res.size(); i++)
            for (size_t j = 0; j < greedy.res[i].size(); j++)
                TSP::local_opt::opt_3(greedy.res[i][j], dist_inner_cluster[i]);

        print(greedy.res, clusters, i);
        checks.push_back(check_matrix(greedy.res, dist_inner_cluster[i].size()));
    }

    cout << "lenght: " << lenght << endl;
    for (size_t check : checks)
        cout << "check: " << check << endl;

}

int main()
{
    //int_matrix t(2, vector<size_t>(2, 0));

    //vector<size_t>& p = t[0];
    //p[0] = 1;

    osman_cluster_test();
    return 0;
}

