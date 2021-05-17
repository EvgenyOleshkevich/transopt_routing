#include "headers/algorithms.hpp"
#include <fstream>
#include <thread>
#include <boost/thread.hpp>
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
const size_t count_clust = 7;

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

namespace print {
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


void print_file(const vector<int_matrix>& routs, const char* name)
{
    std::ofstream out(name);
    for (const int_matrix& rout_mat : routs)
        for (const vector<size_t>& rout : rout_mat)
        {
            for (const size_t vertex : rout)
                out << vertex << " ";
            out << ";";
        }
    out.close();
}

void print_file(const int_matrix& routs, const char* name)
{
    std::ofstream out(name);
    for (const vector<size_t>& rout : routs)
    {
        for (const size_t vertex : rout)
            out << vertex << " ";
        out << ";";
    }
    out.close();
}

void print_file(const vector<int_matrix>& routs, const int_matrix& clusters,
    const size_t clust_id, const char* name)
{
    std::ofstream out(name, std::ios::app);
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

void print_file(const int_matrix& routs, const int_matrix& clusters,
    const size_t clust_id, const char* name)
{
    std::ofstream out(name, std::ios::app);
    for (const vector<size_t>& rout : routs)
    {
        for (const size_t vertex : rout)
            out << clusters[clust_id][vertex] << " ";
        out << ";";
    }
    out << ":";
    out.close();
}
}

namespace checkers {
    bool check_matrix(const vector<int_matrix>& routs, const size_t size)
    {
        vector<size_t> used(size, 0);
        for (const int_matrix& rout_mat : routs)
            for (const vector<size_t>& rout : rout_mat)
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
        //cout << "check is_twice: " << is_twice << endl;
        //cout << "check is_null: " << is_null << endl;
        return !is_twice && !is_null;
    }

    bool check_transport(const int_matrix& routs, const vector<size_t>& transports_id)
    {
        bool is_ok = true;

        for (size_t i = 0; i < routs.size(); ++i)
        {
            double w = 0;
            for (const size_t vertex : routs[i])
                w += weights[vertex];
            is_ok &= transports[transports_id[i]].capacity >= w;
        }

        return is_ok;
    }

    bool check_transport(const int_matrix& routs, const vector<size_t>& transports_id, const vector<double>& W)
    {
        bool is_ok = true;

        for (size_t i = 0; i < routs.size(); ++i)
        {
            double w = 0;
            for (const size_t vertex : routs[i])
                w += W[vertex];
            if (transports[transports_id[i]].capacity < w)
            {
                is_ok &= transports[transports_id[i]].capacity >= w;
            }

        }

        return is_ok;
    }

    bool check_transport(const int_matrix& routs, const Transport& trans, const vector<double>& W)
    {
        bool is_ok = true;

        for (size_t i = 0; i < routs.size(); ++i)
        {
            double w = 0;
            for (const size_t vertex : routs[i])
                w += W[vertex];
            if (trans.capacity < w)
                is_ok &= trans.capacity >= w;
        }

        return is_ok;
    }
}

void save_cluster()
{
    read_file();
    read_transports();
    dist_mat = utils::fill_matrix(x, y);
    auto clusters = clustering::dichotomous_division_weight(dist_mat, weights, clust_capacity, count_clust);

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

void greedy_parallel_test()
{
    read_file();
    read_transports();
    dist_mat = utils::fill_matrix(x, y);
    auto start_time = clock();
    std::ofstream out("res_clust.txt");
    out.close();

    auto clusters = clustering::dichotomous_division_weight(dist_mat, weights, clust_capacity, count_clust);

    auto dist_inner_cluster = clustering::get_dist_inner_cluster(dist_mat, clusters);
    auto weight_inner_cluster = balancedVRP::clustering::get_weight_inner_cluster(weights, clusters);

    vector<std::shared_ptr<boost::thread>> threads;
    vector< project::GreadyBase> greedies;

    project::GreadyBase k0(
        dist_inner_cluster[0],
        utils::fill_sort_matrix(dist_inner_cluster[0]),
        weight_inner_cluster[0],
        transports);
    project::GreadyBase k1(
        dist_inner_cluster[1],
        utils::fill_sort_matrix(dist_inner_cluster[1]),
        weight_inner_cluster[1],
        transports);
    project::GreadyBase k2(
        dist_inner_cluster[2],
        utils::fill_sort_matrix(dist_inner_cluster[2]),
        weight_inner_cluster[2],
        transports);
    project::GreadyBase k3(
        dist_inner_cluster[3],
        utils::fill_sort_matrix(dist_inner_cluster[3]),
        weight_inner_cluster[3],
        transports);
    project::GreadyBase k4(
        dist_inner_cluster[4],
        utils::fill_sort_matrix(dist_inner_cluster[4]),
        weight_inner_cluster[4],
        transports);
    project::GreadyBase k5(
        dist_inner_cluster[5],
        utils::fill_sort_matrix(dist_inner_cluster[5]),
        weight_inner_cluster[5],
        transports);
    project::GreadyBase k6(
        dist_inner_cluster[6],
        utils::fill_sort_matrix(dist_inner_cluster[6]),
        weight_inner_cluster[6],
        transports);


    threads.push_back(std::make_shared<boost::thread>(&project::GreadyBase::run, &k0));
    threads.push_back(std::make_shared<boost::thread>(&project::GreadyBase::run, &k1));
    threads.push_back(std::make_shared<boost::thread>(&project::GreadyBase::run, &k2));
    threads.push_back(std::make_shared<boost::thread>(&project::GreadyBase::run, &k3));
    threads.push_back(std::make_shared<boost::thread>(&project::GreadyBase::run, &k4));
    threads.push_back(std::make_shared<boost::thread>(&project::GreadyBase::run, &k5));
    threads.push_back(std::make_shared<boost::thread>(&project::GreadyBase::run, &k6));

    for (auto t : threads) {
        t->join();
    }
    // time: 1692
    cout << "time: " << clock() - start_time;
}

void clurk_test()
{
    const char* name = "tests/clurk_res.txt";
    std::ofstream out(name);
    out.close();
    read_file();
    read_transports();
    dist_mat = utils::fill_matrix(x, y);

    unsigned int start_time = clock();
    auto clusters = clustering::dichotomous_division_weight(dist_mat, weights, clust_capacity, count_clust);

    auto dist_inner_cluster = clustering::get_dist_inner_cluster(dist_mat, clusters);
    auto weight_inner_cluster = clustering::get_weight_inner_cluster(weights, clusters);


    double length = 0;
    vector<size_t> checks;
    vector<project::ClarkRight> reults;
    for (size_t i = 0; i < clusters.size(); i++)
    {
        project::ClarkRight clurk(dist_inner_cluster[i], weight_inner_cluster[i], transports);
        clurk.run();
        reults.push_back(clurk);
    }

    cout << "# clurk" << endl;
    cout << "# time: " << clock() - start_time << " millisec " << endl;

    for (size_t i = 0; i < clusters.size(); i++)
    {

        length += reults[i].length();
        print::print_file(reults[i].res, clusters, i, name);
        checks.push_back(checkers::check_matrix(reults[i].res, dist_inner_cluster[i].size()));
        for (size_t j = 0; j < reults[i].res.size(); j++)
            checks.push_back(checkers::check_transport(reults[i].res[j], transports[j], weight_inner_cluster[i]));
    }

    cout << "# length: " << length << endl;
    for (size_t check : checks)
        cout << "check: " << check << endl;
}

void sweep_test()
{
    const char* name = "tests/sweeping_res.txt";
    std::ofstream out(name);
    out.close();
    read_file();
    read_transports();
    dist_mat = utils::fill_matrix(x, y);

    unsigned int start_time = clock();
    auto clusters = clustering::dichotomous_division_weight(dist_mat, weights, clust_capacity, count_clust);

    auto dist_inner_cluster = clustering::get_dist_inner_cluster(dist_mat, clusters);
    auto weight_inner_cluster = clustering::get_weight_inner_cluster(weights, clusters);
    auto x_coord = clustering::get_weight_inner_cluster(x, clusters);
    auto y_coord = clustering::get_weight_inner_cluster(y, clusters);


    double length = 0;
    vector<size_t> checks;
    vector<project::Sweeping> reults;
    for (size_t i = 0; i < clusters.size(); i++)
    {
        project::Sweeping sweep(dist_inner_cluster[i], x_coord[i], y_coord[i], weight_inner_cluster[i], transports);
        sweep.run();
        reults.push_back(sweep);
    }

    cout << "# sweeping" << endl;
    cout << "# time: " << clock() - start_time << " millisec " << endl;

    for (size_t i = 0; i < clusters.size(); i++)
    {

        length += reults[i].length();
        print::print_file(reults[i].res, clusters, i, name);
        checks.push_back(checkers::check_matrix(reults[i].res, dist_inner_cluster[i].size()));
        for (size_t j = 0; j < reults[i].res.size(); j++)
            checks.push_back(checkers::check_transport(reults[i].res[j], transports[j], weight_inner_cluster[i]));
    }

    cout << "# length: " << length << endl;
    for (size_t check : checks)
        cout << "check: " << check << endl;
}

void greedy_test()
{
    const char* name = "tests/greedy_res.txt";
    std::ofstream out(name);
    out.close();
    read_file();
    read_transports();
    dist_mat = utils::fill_matrix(x, y);

    unsigned int start_time = clock();
    auto clusters = clustering::dichotomous_division_weight(dist_mat, weights, clust_capacity, count_clust);
    
    auto dist_inner_cluster = clustering::get_dist_inner_cluster(dist_mat, clusters);
    auto weight_inner_cluster = clustering::get_weight_inner_cluster(weights, clusters);


    double length = 0;
    vector<size_t> checks;
    for (size_t i = 0; i < clusters.size(); i++)
    {
        auto sort_matrix = utils::fill_sort_matrix(dist_inner_cluster[i]);

        project::GreadyBase greedy(
            dist_inner_cluster[i],
            sort_matrix,
            weight_inner_cluster[i],
            transports);
        greedy.run();
        length += greedy.length();
        print::print_file(greedy.res, clusters, i, name);
        checks.push_back(checkers::check_matrix(greedy.res, dist_inner_cluster[i].size()));
    }

    cout << "# greedy" << endl;
    cout << "# time: " << clock() - start_time << " millisec "<< endl;
    cout << "# length: " << length << endl;
    for (size_t check : checks)
        cout << "check: " << check << endl;
    
}

void ant_test()
{
    const char* name = "tests/ant_res.txt";
    std::ofstream out(name);
    out.close();
    read_file();
    read_transports();
    dist_mat = utils::fill_matrix(x, y);

    unsigned int start_time = clock();
    auto clusters = clustering::dichotomous_division_weight(dist_mat, weights, clust_capacity, count_clust);

    auto dist_inner_cluster = clustering::get_dist_inner_cluster(dist_mat, clusters);
    auto weight_inner_cluster = clustering::get_weight_inner_cluster(weights, clusters);


    double length = 0;
    vector<size_t> checks;
    vector<Ant_algorithm> reults;
    for (size_t i = 0; i < clusters.size(); i++)
    {
        Ant_algorithm ant(dist_inner_cluster[i], weight_inner_cluster[i], transports);
        ant.run();
        reults.push_back(ant);
    }

    cout << "# ant" << endl;
    cout << "# time: " << clock() - start_time << " millisec " << endl;

    for (size_t i = 0; i < clusters.size(); i++)
    {
        length += reults[i].length();
        print::print_file(reults[i].res, clusters, i, name);
        checks.push_back(checkers::check_matrix(reults[i].res, dist_inner_cluster[i].size()));
        for (size_t j = 0; j < reults[i].res.size(); j++)
            checks.push_back(checkers::check_transport(reults[i].res[j], transports[j], weight_inner_cluster[i]));
    }

    cout << "# length: " << length << endl;
    for (size_t check : checks)
        cout << "check: " << check << endl;
}

void osman_test()
{
    const char* name = "tests/osman_res.txt";
    std::ofstream out(name);
    out.close();
    read_file();
    read_transports();
    dist_mat = utils::fill_matrix(x, y);

    unsigned int start_time = clock();
    auto clusters = clustering::dichotomous_division_weight(dist_mat, weights, clust_capacity, count_clust);

    auto dist_inner_cluster = clustering::get_dist_inner_cluster(dist_mat, clusters);
    auto weight_inner_cluster = clustering::get_weight_inner_cluster(weights, clusters);

    double length = 0;
    vector<size_t> checks;
    vector<project::Osman> reults;
    for (size_t i = 0; i < clusters.size(); i++)
    {
        Ant_algorithm ant(dist_inner_cluster[i], weight_inner_cluster[i], transports);
        ant.run();

        auto data = project::Osman::transform_data(ant.res);

        project::Osman osman(dist_inner_cluster[i],
            weight_inner_cluster[i],
            transports, data.first,
            data.second, data.first.size());

        osman.run();
        osman.add_zero_vertex();
        reults.push_back(osman);
    }

    cout << "# osman" << endl;
    cout << "# time: " << clock() - start_time << " millisec " << endl;

    for (size_t i = 0; i < clusters.size(); i++)
    {
        length += reults[i].length();
        print::print_file(reults[i].res, clusters, i, name);
        checks.push_back(checkers::check_matrix(reults[i].res, dist_inner_cluster[i].size()));
        checks.push_back(checkers::check_transport(reults[i].res, reults[i].transport_id, weight_inner_cluster[i]));
    }

    cout << "# length: " << length << endl;
    for (size_t check : checks)
        cout << "check: " << check << endl;
}

void global_local_search_cluster_test(const vector<int_matrix>& vec_routs,
    const vector<vector<size_t>>& vec_trans, const int_matrix& clusters)
{
    int_matrix routs;
    vector<size_t> trans_id;
    vector<size_t> clust_id;

    for (size_t i = 0; i < clusters.size(); i++)
        for (size_t t = 0; t < vec_routs[i].size(); t++)
        {
            const vector<size_t>& rout_i = vec_routs[i][t];
            vector<size_t> rout(rout_i.size());
            trans_id.push_back(vec_trans[i][t]);
            clust_id.push_back(i);

            for (size_t v = 0; v < rout_i.size(); v++)
                rout[v] = clusters[i][rout_i[v]];
            routs.push_back(rout);
        }

    project::GlobalWidhtNeighborhoodSearch local(dist_mat, weights,
        transports, routs, trans_id, clust_id, frequence);

    local.run();
    cout << "length: " << local.length() << endl;
    cout << "check: " << checkers::check_transport(local.res, trans_id) << endl;

    print::print_file(local.res);
}

void local_search_cluster_test_old()
{
    std::ofstream out("res_clust.txt");
    out.close();
    read_file();
    read_transports();
    dist_mat = utils::fill_matrix(x, y);
    auto clusters = balancedVRP::clustering::dichotomous_division_weight(dist_mat, weights, clust_capacity, count_clust);
    clusters = project::clusters_with_multiple_duplicates(clusters, frequence);
    auto dist_inner_cluster = clustering::get_dist_inner_cluster(dist_mat, clusters);
    auto weight_inner_cluster = clustering::get_weight_inner_cluster(weights, clusters);

    vector<size_t> checks;

    vector<int_matrix> vec_routs;
    vector<vector<size_t>> vec_trans;
    double sum_gready = 0;
    double sum_osman = 0;
    double sum_local= 0;
    
    for (size_t i = 0; i < clusters.size(); i++)
    {
        auto sort_matrix = utils::fill_sort_matrix(dist_inner_cluster[i]);

        project::GreadyBase greedy(
            dist_inner_cluster[i],
            sort_matrix,
            weight_inner_cluster[i],
            transports);
        greedy.run();
        double len = greedy.length();
        cout << "greedy length: " << len << endl;
        sum_gready += len;

        auto data = project::Osman::transform_data(greedy.res);
        cout << "check transport greedy: " << checkers::check_transport(data.first, data.second, weight_inner_cluster[i]) << endl;


        project::Osman osman(dist_inner_cluster[i],
            weight_inner_cluster[i],
            transports, data.first,
            data.second, data.first.size());
        osman.run();
        osman.add_zero_vertex();
        len = osman.length();
        cout << "osman length: " << len << endl;
        sum_osman += len;
        cout << "check transport osman: " << checkers::check_transport(osman.res, osman.transport_id, weight_inner_cluster[i]) << endl;

        project::WidhtNeighborhoodSearch local(dist_inner_cluster[i],
            weight_inner_cluster[i],
            transports, osman.res,
            osman.transport_id);

        local.run();
        len = local.length();
        cout << "local length: " << len << endl;
     

        sum_local += len;
        print::print_file(local.res, clusters, i);
        cout << "check transport local: " << checkers::check_transport(local.res, local.transport_id, weight_inner_cluster[i]) << endl;

        checks.push_back(checkers::check_matrix(local.res, dist_inner_cluster[i].size()));

        vec_routs.push_back(local.res);
        vec_trans.push_back(local.transport_id);
    }

    cout << "sum_gready length: " << sum_gready << endl;
    cout << "sum_osman length: " << sum_osman << endl;
    cout << "sum_local length: " << sum_local << endl;
    for (size_t check : checks)
        cout << "check: " << check << endl;

    global_local_search_cluster_test(vec_routs, vec_trans, clusters);
}

void final_test()
{
    const char* name = "tests/final_res.txt";
    std::ofstream out(name);
    out.close();
    read_file();
    read_transports();
    dist_mat = utils::fill_matrix(x, y);

    unsigned int start_time = clock();
    auto clusters = clustering::dichotomous_division_weight(dist_mat, weights, clust_capacity, count_clust);

    auto dist_inner_cluster = clustering::get_dist_inner_cluster(dist_mat, clusters);
    auto weight_inner_cluster = clustering::get_weight_inner_cluster(weights, clusters);

    double length = 0;
    vector<size_t> checks;

    vector<int_matrix> vec_routs;
    vector<vector<size_t>> vec_trans;
    for (size_t i = 0; i < clusters.size(); i++)
    {
        Ant_algorithm ant(dist_inner_cluster[i], weight_inner_cluster[i], transports);
        ant.run();

        auto data = project::Osman::transform_data(ant.res);

        project::Osman osman(dist_inner_cluster[i],
            weight_inner_cluster[i],
            transports, data.first,
            data.second, data.first.size());

        osman.run();
        osman.add_zero_vertex();

        project::WidhtNeighborhoodSearch local(dist_inner_cluster[i],
            weight_inner_cluster[i],
            transports, osman.res,
            osman.transport_id);

        local.run();

        vec_routs.push_back(local.res);
        vec_trans.push_back(local.transport_id);
    }

    int_matrix routs;
    vector<size_t> trans_id;
    vector<size_t> clust_id;

    for (size_t i = 0; i < clusters.size(); i++)
        for (size_t t = 0; t < vec_routs[i].size(); t++)
        {
            const vector<size_t>& rout_i = vec_routs[i][t];
            vector<size_t> rout(rout_i.size());
            trans_id.push_back(vec_trans[i][t]);
            clust_id.push_back(i);

            for (size_t v = 0; v < rout_i.size(); v++)
                rout[v] = clusters[i][rout_i[v]];
            routs.push_back(rout);
        }

    project::GlobalWidhtNeighborhoodSearch local(dist_mat, weights,
        transports, routs, trans_id, clust_id, frequence);

    local.run();

    cout << "# final" << endl;
    cout << "# time: " << clock() - start_time << " millisec " << endl;
    cout << "# length: " << local.length() << endl;
    cout << "check: " << checkers::check_transport(local.res, trans_id) << endl;

    print::print_file(local.res, name);
}

int main()
{
    final_test();
    return 0;
}
