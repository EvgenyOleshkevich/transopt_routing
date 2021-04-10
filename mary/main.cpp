#include "headers/algorithms.hpp"
#include <fstream>


using namespace std;
using namespace algorithms;

const size_t count_point = 51;
const size_t need_routs = 5;
vector<double> x { 0, 28, -7, 78, 23, -99, -96, -39, 45, -56, -2, -77, -25, 43, -72,
-31, 78, -49, -76, 23, -33, -22, -58, -79, 24, -6, -92, 41, -42, 92, -72, -25, -3,
74, 44, 46, -16, -2, 36, -21, -69, 19, 10, 39, 40, -84, 46, 15, -86, -65, 9 };
vector<double> y { 0, 62, 80, 60, 83, 36, 62, -93, 95, -22, 42, 77, -34, -66, -43,
-69, 98, 94, -15, -34, -51, 3, -66, -69, 90, -53, 75, -28, 53, 46, -7, 44, 62,
-23, 48, -22, -90, -45, -50, 46, -18, 30, -28, 10, 43, 83, -84, 64, 45, 7, -63 };
vector<double> weight{ 0, 49, 68, 34, 27, 51, 3, 8, 76, 18, 76, 79, 41, 68,
79, 87, 14, 17, 39, 47, 30, 2, 56, 97, 37, 54, 39, 45, 88, 31, 65, 49, 54, 92, 8,
83, 31, 65, 80, 66, 40, 85, 10, 64, 40, 54, 86, 80, 63, 31, 24 };

vector<vector<double>> dist_mat;
// smth
void c()
{
    vector<size_t> g1 = { 42 };
    vector<size_t> g2 = { 37, 25, 50, 46, 13, 38, 33, 35, 27, 19 };
    vector<size_t> g3 = { 41, 44, 34, 29, 3, 16, 8, 24, 4, 1, 47, 32, 2, 17, 11, 45, 26, 6, 5, 48, 28, 31, 39, 10 };
    vector<size_t> g4 = { 12, 20, 15, 36, 7, 22, 23, 14, 9, 40, 18, 30, 49, 21 };
    vector<size_t> g5 = { 43 };
    vector<vector<size_t>> o = { g1, g2, g3, g4, g5 };
    cout << utils::length_routs(o, dist_mat);
}

void c2()
{
    srand((unsigned int)time(0));
    dist_mat = utils::fill_matrix(x, y, count_point);
    auto rout = TSP::Lin_Kernighan::Lin_Kernighan(x, y, 51);
    TSP::local_opt::TSP_2_opt(rout, dist_mat);
    TSP::local_opt::TSP_3_opt(rout, dist_mat);
    TSP::local_opt::TSP_2_opt(rout, dist_mat);
    TSP::local_opt::TSP_3_opt(rout, dist_mat);
    balancedVRP::VND_STS a(dist_mat, need_routs, 600, weight);

    std::vector<std::pair<size_t, double>> t;
    for (size_t v : rout)
    {
        t.push_back({ v, weight[v] });
    }

    auto routs = a.calculate(t).first;

    cout << "dynamic: lenght= " << utils::length_routs(routs, dist_mat) << endl;
    for (const vector<size_t>& rout : routs)
    {
        double w = 0;
        cout << "[" << 0 << ", ";
        for (auto vertex : rout)
        {
            w += weight[vertex];
            cout << vertex << ", ";
        }
        cout << 0 << "]," << endl;
    }
}

void read_file()
{
    x.clear();
    y.clear();
    weight.clear();

    x.push_back(0);
    y.push_back(0);
    weight.push_back(0);
    std::ifstream in("file.csv");
    while (!in.eof())
    {
        string line;
        getline(in, line);
        double w, X, Y;
        in >> w >> X >> Y;
        x.push_back(X);
        y.push_back(Y);
        weight.push_back(w);
    }

    in.close();

    auto dist = utils::fill_matrix_with_end_point(x, y);

    auto clusters  = balancedVRP::clustering::dichotomous_division(dist, 4);
    auto dist_matrixes = balancedVRP::clustering::get_dist_inner_cluster(dist, clusters);
    auto weights = balancedVRP::clustering::get_weight_inner_cluster(weight, clusters);
    double lenght = 0;
    for (size_t i = 0; i < dist_matrixes.size(); ++i) {
        

        balancedVRP::VND_STS a(dist_matrixes[i], 60, 24, weight);

        vector<size_t> rout(clusters[i].size());
        for (size_t i = 0; i < rout.size(); ++i)
            rout[i] = i + 1;

        unsigned int start_time = clock(); // начальное время
        TSP::local_opt::TSP_3_opt_fast(rout, dist_matrixes[i]);
        unsigned int end_time = clock(); // конечное время
        unsigned int search_time = end_time - start_time; // искомое время
        cout << "time: " << search_time << endl;

        TSP::local_opt::TSP_3_opt(rout, dist_matrixes[i]);
        TSP::local_opt::TSP_2_opt(rout, dist_matrixes[i]);
        TSP::local_opt::TSP_3_opt(rout, dist_matrixes[i]);

        std::vector<std::pair<size_t, double>> t;
        for (size_t v : rout)
        {
            t.push_back({ v, weights[i][v] });
        }

        auto routs = a.calculate(t).first;

        //auto routs = algorithms::clark_right::clark_right(dist_matrixes[i], dist_matrixes[i].size(), 60);
        // auto osman = new Osman(dist_matrixes[i], dist_matrixes[i].size(), 60);
        lenght += utils::length_routs(routs, dist_matrixes[i]);
    }
    cout << "lenght: " << lenght << endl;
}

// lenght: 146182
// time : 506.847 sec

// lenght: 186045
// time : 94.158 sec


// time: 259378
int main()
{
    unsigned int start_time = clock(); // начальное время
    read_file();
    unsigned int end_time = clock(); // конечное время
    unsigned int search_time = end_time - start_time; // искомое время
    cout << "time: " << search_time;
    return 0;
    srand((unsigned int)time(0));
    dist_mat = utils::fill_matrix(x, y, count_point);

    auto routs_dichotomous_division = balancedVRP::clustering::dichotomous_division_weight(dist_mat, weight, 600, 5);
    //auto routs_osman = osman::osman(routs_dichotomous_division, new osman::checker_change_smoll_local());
    auto TSP_opt_rout = routs_dichotomous_division;
    for (size_t i = 0; i < TSP_opt_rout.size(); ++i)
    {
        TSP::local_opt::TSP_2_opt(TSP_opt_rout[i], dist_mat);
        TSP::local_opt::TSP_3_opt(TSP_opt_rout[i], dist_mat);
        TSP::local_opt::TSP_2_opt(TSP_opt_rout[i], dist_mat);
        TSP::local_opt::TSP_3_opt(TSP_opt_rout[i], dist_mat);
    }

    cout << "dichotomous_division: lenght= " << utils::length_routs(TSP_opt_rout, dist_mat) << endl;
    for (const vector<size_t>& rout : TSP_opt_rout)
    {
        cout << "[" << 0 << ", ";
        for (auto vertex : rout)
            cout << vertex << ", ";
        cout << 0 << "]," << endl;
    }
    /*auto rout = TSP::Lin_Kernighan::Lin_Kernighan(x, y, 51);
    TSP::local_opt::TSP_2_opt(rout, dist_mat);
    TSP::local_opt::TSP_3_opt(rout, dist_mat);
    TSP::local_opt::TSP_2_opt(rout, dist_mat);
    TSP::local_opt::TSP_3_opt(rout, dist_mat);
    cout << "[" << 0 << ", ";
    for (auto vertex : rout)
        cout << vertex << ", ";
    cout << 0 << "]," << endl;
    std::cout << "Lin_Kernighan: lenght= " <<
        utils::length_rout(rout, dist_mat) << std::endl;

    auto cutting_routs = balancedVRP::cutting_rout(rout, dist_mat, 6);

    for (const vector<size_t>& rout : cutting_routs)
    {
        cout << "[" << 0 << ", ";
        for (auto vertex : rout)
            cout << vertex << ", ";
        cout << 0 << "]," << endl;
    }
    cout << "cutting_rout: lenght= " << utils::length_routs(cutting_routs, dist_mat) << endl;

    balancedVRP::clustering::radian_sort(x, y, 51);
    //auto routs_dichotomous_division = balancedVRP::clustering::dichotomous_division(dist_mat, 5);
    auto routs_dichotomous_division = balancedVRP::clustering::sweeping(x, y, 50, 5);
    //auto routs_osman = osman::osman(routs_dichotomous_division, new osman::checker_change_smoll_local());
    auto TSP_opt_rout = routs_dichotomous_division;
    for (size_t i = 0; i < TSP_opt_rout.size(); ++i)
    {
        TSP::local_opt::TSP_2_opt(TSP_opt_rout[i], dist_mat);
        TSP::local_opt::TSP_3_opt(TSP_opt_rout[i], dist_mat);
        TSP::local_opt::TSP_2_opt(TSP_opt_rout[i], dist_mat);
        TSP::local_opt::TSP_3_opt(TSP_opt_rout[i], dist_mat);
    }

    cout << "dichotomous_division: lenght= " << utils::length_routs(TSP_opt_rout, dist_mat) << endl;
    for (const vector<size_t>& rout : TSP_opt_rout)
    {
        cout << "[" << 0 << ", ";
        for (auto vertex : rout)
            cout << vertex << ", ";
        cout << 0 << "]," << endl;
    }

    auto routs_genetic_alg = genetic::genetic_alg();
    cout << "genetic_alg: lenght= " << length_routs(routs_genetic_alg) << endl;
    for (const vector<size_t>& rout : routs_genetic_alg)
    {
        cout << "[" << 0 << ", ";
        for (auto vertex : rout)
            cout << vertex << ", ";
        cout << 0 << "]," << endl;
    }
    
    auto routs_clark_right = algorithms::clark_right::clark_right(dist_mat, count_point, 5);
    cout << "clark_right: lenght= " << utils::length_routs(routs_clark_right, dist_mat) <<endl;
    for (const vector<size_t>& rout : routs_clark_right)
    {
        cout << "[" << 0 << ", ";
        for (auto vertex : rout)
            cout << vertex << ", ";
        cout << 0 << "]," << endl;
    }
    algorithms::Osman osman(dist_mat, count_point, 5);
    auto routs_osman = osman.start(routs_clark_right, 1);
    cout << "osman: lenght= " << utils::length_routs(routs_osman, dist_mat) << endl;
    for (const vector<size_t>& rout : routs_osman)
    {
        cout << "[" << 0 << ", ";
        for (auto vertex : rout)
            cout << vertex << ", ";
        cout << 0 << "]," << endl;
    }

    TSP_opt_rout = routs_osman;
    for (size_t i = 0; i < TSP_opt_rout.size(); ++i)
    {
        TSP::local_opt::TSP_2_opt(TSP_opt_rout[i], dist_mat);
        TSP::local_opt::TSP_3_opt(TSP_opt_rout[i], dist_mat);
        TSP::local_opt::TSP_2_opt(TSP_opt_rout[i], dist_mat);
        TSP::local_opt::TSP_3_opt(TSP_opt_rout[i], dist_mat);
    }

    cout << "TSP_opt_rout: lenght= " << utils::length_routs(TSP_opt_rout, dist_mat) << endl;
    for (const vector<size_t>& rout : TSP_opt_rout)
    {
        cout << "[" << 0 << ", ";
        for (auto vertex : rout)
            cout << vertex << ", ";
        cout << 0 << "]," << endl;
    }*/
}
