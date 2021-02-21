#include "headers/algorithms.hpp"


using namespace std;
using namespace algorithms;

const size_t count_point = 51;
const size_t need_routs = 5;
const double x[] = { 0, 28, -7, 78, 23, -99, -96, -39, 45, -56, -2, -77, -25, 43, -72,
-31, 78, -49, -76, 23, -33, -22, -58, -79, 24, -6, -92, 41, -42, 92, -72, -25, -3,
74, 44, 46, -16, -2, 36, -21, -69, 19, 10, 39, 40, -84, 46, 15, -86, -65, 9 };
const double y[] = { 0, 62, 80, 60, 83, 36, 62, -93, 95, -22, 42, 77, -34, -66, -43,
-69, 98, 94, -15, -34, -51, 3, -66, -69, 90, -53, 75, -28, 53, 46, -7, 44, 62,
-23, 48, -22, -90, -45, -50, 46, -18, 30, -28, 10, 43, 83, -84, 64, 45, 7, -63 };
const vector<double> weight{ 0, 49, 68, 34, 27, 51, 3, 8, 76, 18, 76, 79, 41, 68,
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

int main()
{
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
