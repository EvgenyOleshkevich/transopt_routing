#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <random>
#include "headers/balancedVRP.hpp"

using namespace std;

const size_t count_point = 51;
const size_t need_routs = 5;
const double x[51] = { 0, 28, -7, 78, 23, -99, -96, -39, 45, -56, -2, -77, -25, 43, -72,
-31, 78, -49, -76, 23, -33, -22, -58, -79, 24, -6, -92, 41, -42, 92, -72, -25, -3,
74, 44, 46, -16, -2, 36, -21, -69, 19, 10, 39, 40, -84, 46, 15, -86, -65, 9 };
const double y[51] = { 0, 62, 80, 60, 83, 36, 62, -93, 95, -22, 42, 77, -34, -66, -43,
-69, 98, 94, -15, -34, -51, 3, -66, -69, 90, -53, 75, -28, 53, 46, -7, 44, 62,
-23, 48, -22, -90, -45, -50, 46, -18, 30, -28, 10, 43, 83, -84, 64, 45, 7, -63 };


vector<vector<double>> dist_mat;

double sqr(double x)
{
    return x * x;
}

void fill_matrix()
{
    dist_mat = vector<vector<double>>(count_point, vector<double>(count_point));
    for (size_t i = 0; i < count_point; ++i)
        for (size_t j = 0; j < count_point; ++j)
            dist_mat[i][j] = std::sqrt(sqr(x[i] - x[j]) + sqr(y[i] - y[j]));
}

double length_rout(const vector<size_t>& rout)
{
    double lenght = dist_mat[0][rout[0]] + dist_mat[rout.back()][0];
    for (size_t i = 1; i < rout.size(); ++i)
        lenght += dist_mat[rout[i - 1]][rout[i]];
    return lenght;
}

double length_routs(const vector<vector<size_t>>& routs)
{
    double lenght = 0;
    for (const vector<size_t>& rout : routs)
        lenght += length_rout(rout);
    return lenght;
}

// smth
void c()
{
    vector<size_t> g1 = { 42 };
    vector<size_t> g2 = { 37, 25, 50, 46, 13, 38, 33, 35, 27, 19 };
    vector<size_t> g3 = { 41, 44, 34, 29, 3, 16, 8, 24, 4, 1, 47, 32, 2, 17, 11, 45, 26, 6, 5, 48, 28, 31, 39, 10 };
    vector<size_t> g4 = { 12, 20, 15, 36, 7, 22, 23, 14, 9, 40, 18, 30, 49, 21 };
    vector<size_t> g5 = { 43 };
    vector<vector<size_t>> o = { g1, g2, g3, g4, g5 };
    cout << length_routs(o);
}


namespace clark_right
{
    vector<vector<double>> get_save_matrix()
    {
        auto save_matrix = vector<vector<double>>(count_point, vector<double>(count_point));
        for (size_t i = 0; i < count_point; ++i)
            for (size_t j = 0; j < count_point; ++j)
                save_matrix[i][j] = dist_mat[i][0] + dist_mat[0][j] - dist_mat[i][j];
        return save_matrix;
    }

    vector<vector<double>> sort_save_matrix(const vector<vector<double>>& save_matrix)
    {
        vector<vector<double>> sorted_save_matrix;
        for (size_t i = 1; i < count_point; ++i)
            for (size_t j = 1; j < count_point; ++j)
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
        for (size_t i = 1; i < sorted_save_matrix.size() / 2; ++i)
            swap(sorted_save_matrix[i], sorted_save_matrix[sorted_save_matrix.size() - i]);
        //for (auto save :sorted_save_matrix)
        //    cout << save[0] << "; i=" << save[1] << "; j=" << save[2] << endl;
        return sorted_save_matrix;
    }

    vector<vector<size_t>> init_routs()
    {
        auto routs = vector<vector<size_t>>(count_point);
        for (size_t i = 1; i < count_point; ++i)
            routs[i].push_back(i);
        return routs;
    }


    vector<vector<size_t>> clark_right()
    {
        auto s = get_save_matrix();
        auto sotred_save_matrix = sort_save_matrix(s);
        auto routs_start = init_routs();
        auto routs_end = init_routs();
        size_t count_routs = count_point - 1;
        while (count_routs != need_routs)
        {
            size_t i = sotred_save_matrix[0][1];
            size_t j = sotred_save_matrix[0][2];
            sotred_save_matrix.erase(sotred_save_matrix.begin());
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


        vector<vector<size_t>> res_rout;
        for (const vector<size_t>& rout : routs_start)
        {
            if (rout.empty())
                continue;
            res_rout.push_back(rout);
        }
        return res_rout;
    }
}

namespace osman
{
    vector<vector<int>> get_move_table()
    {
        return vector<vector<int>>(need_routs, vector<int>(count_point, -1000000000));
    }


    vector<vector<size_t>> BSTM;
    vector<vector<pair<pair<size_t, size_t>, size_t>>> RECM;
    size_t ts;// длина списка запрета
    vector<vector<int>> move_table; // rout*vertex
    // вершины и способ обмена
    // 0 - обмен, 1 - левая в правую перед указанной, 2 - правая в левую перед указанной
    // 3 - левая в правую в конец, 4 - правая в левую в конец, 5 - ничего


    class checker_change
    {
    public:
        void virtual check_change(const vector<vector<size_t>>& routs, size_t p, size_t q) const = 0;
    };

    class checker_change_widht_local : public checker_change
    {
    public:
        void check_change(const vector<vector<size_t>>& routs, size_t p, size_t q) const override
        {
            if (p > q)
                swap(p, q);
            double base_lenght = length_rout(routs[p]);
            base_lenght += length_rout(routs[q]);

            double best_lenght = base_lenght * 1.2;

            auto rout1 = routs[p];
            auto rout2 = routs[q];
            pair<pair<size_t, size_t>, size_t> decision = { {0,0}, 5 };

            // обмен
            for (size_t i = 0; i < rout1.size(); ++i)
                for (size_t j = 0; j < rout2.size(); ++j)
                {
                    swap(rout1[i], rout2[j]);
                    double lenght = length_rout(rout1);
                    lenght += length_rout(rout2);
                    if (best_lenght > lenght)
                    {
                        best_lenght = lenght;
                        decision = { {i,j}, 0 };
                    }
                    swap(rout1[i], rout2[j]);
                }

            // 1 идет в 2
            if (rout1.size() > 1)
            {
                for (size_t i = 0; i < rout1.size(); ++i)
                    for (size_t j = 0; j < rout2.size(); ++j)
                    {
                        rout2.emplace(rout2.begin() + j, rout1[i]);
                        rout1.erase(rout1.begin() + i);
                        double lenght = length_rout(rout1);
                        lenght += length_rout(rout2);
                        if (best_lenght > lenght)
                        {
                            best_lenght = lenght;
                            decision = { {i,j}, 1 };
                        }
                        rout1.emplace(rout1.begin() + i, rout2[j]);
                        rout2.erase(rout2.begin() + j);
                    }

                // вставка в конец
                for (size_t i = 0; i < rout1.size(); ++i)
                {
                    rout2.push_back(rout1[i]);
                    rout1.erase(rout1.begin() + i);
                    double lenght = length_rout(rout1);
                    lenght += length_rout(rout2);
                    if (best_lenght > lenght)
                    {
                        best_lenght = lenght;
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
                    for (size_t j = 0; j < rout1.size(); ++j)
                    {
                        rout1.emplace(rout1.begin() + j, rout2[i]);
                        rout2.erase(rout2.begin() + i);
                        double lenght = length_rout(rout1);
                        lenght += length_rout(rout2);
                        if (best_lenght > lenght)
                        {
                            best_lenght = lenght;
                            decision = { {j,i}, 2 };
                        }
                        rout2.emplace(rout2.begin() + i, rout1[j]);
                        rout1.erase(rout1.begin() + j);
                    }

                // вставка в конец
                for (size_t i = 0; i < rout2.size(); ++i)
                {
                    rout1.push_back(rout2[i]);
                    rout2.erase(rout2.begin() + i);
                    double lenght = length_rout(rout1);
                    lenght += length_rout(rout2);
                    if (best_lenght > lenght)
                    {
                        best_lenght = lenght;
                        decision = { {0,i}, 4 };
                    }
                    rout2.emplace(rout2.begin() + i, rout1.back());
                    rout1.pop_back();
                }
            }

            BSTM[p][q] = best_lenght - base_lenght;
            BSTM[q][p] = 0;
            RECM[p][q] = decision;
            RECM[q][p] = { {0, 0}, 6 };
        }
    };

    class checker_change_smoll_local : public checker_change
    {
    public:
        void check_change(const vector<vector<size_t>>& routs, size_t p, size_t q) const override
        {
            if (p > q)
                swap(p, q);
            double base_lenght = length_rout(routs[p]);
            base_lenght += length_rout(routs[q]);

            double best_lenght = base_lenght * 1.2;

            auto rout1 = routs[p];
            auto rout2 = routs[q];
            pair<pair<size_t, size_t>, size_t> decision = { {0,0}, 5 };

            // обмен
            for (size_t i = 0; i < rout1.size(); ++i)
                for (size_t j = 0; j < rout2.size(); ++j)
                {
                    swap(rout1[i], rout2[j]);
                    double lenght = length_rout(rout1);
                    lenght += length_rout(rout2);
                    if (best_lenght > lenght)
                    {
                        best_lenght = lenght;
                        decision = { {i,j}, 0 };
                    }
                    swap(rout1[i], rout2[j]);
                }

            BSTM[p][q] = best_lenght - base_lenght;
            BSTM[q][p] = 0;
            RECM[p][q] = decision;
            RECM[q][p] = { {0, 0}, 6 };
        }
    };


    void check_change(const vector<vector<size_t>>& routs, size_t p, size_t q)
    {
        if (p > q)
            swap(p, q);
        double base_lenght = length_rout(routs[p]);
        base_lenght += length_rout(routs[q]);

        double best_lenght = base_lenght * 1.2;

        auto rout1 = routs[p];
        auto rout2 = routs[q];
        pair<pair<size_t, size_t>, size_t> decision = { {0,0}, 5 };

        // обмен
        for (size_t i = 0; i < rout1.size(); ++i)
            for (size_t j = 0; j < rout2.size(); ++j)
            {
                swap(rout1[i], rout2[j]);
                double lenght = length_rout(rout1);
                lenght += length_rout(rout2);
                if (best_lenght > lenght)
                {
                    best_lenght = lenght;
                    decision = { {i,j}, 0 };
                }
                swap(rout1[i], rout2[j]);
            }

        // 1 идет в 2
        if (rout1.size() > 1)
        {
            for (size_t i = 0; i < rout1.size(); ++i)
                for (size_t j = 0; j < rout2.size(); ++j)
                {
                    rout2.emplace(rout2.begin() + j, rout1[i]);
                    rout1.erase(rout1.begin() + i);
                    double lenght = length_rout(rout1);
                    lenght += length_rout(rout2);
                    if (best_lenght > lenght)
                    {
                        best_lenght = lenght;
                        decision = { {i,j}, 1 };
                    }
                    rout1.emplace(rout1.begin() + i, rout2[j]);
                    rout2.erase(rout2.begin() + j);
                }

            // вставка в конец
            for (size_t i = 0; i < rout1.size(); ++i)
            {
                rout2.push_back(rout1[i]);
                rout1.erase(rout1.begin() + i);
                double lenght = length_rout(rout1);
                lenght += length_rout(rout2);
                if (best_lenght > lenght)
                {
                    best_lenght = lenght;
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
                for (size_t j = 0; j < rout1.size(); ++j)
                {
                    rout1.emplace(rout1.begin() + j, rout2[i]);
                    rout2.erase(rout2.begin() + i);
                    double lenght = length_rout(rout1);
                    lenght += length_rout(rout2);
                    if (best_lenght > lenght)
                    {
                        best_lenght = lenght;
                        decision = { {j,i}, 2 };
                    }
                    rout2.emplace(rout2.begin() + i, rout1[j]);
                    rout1.erase(rout1.begin() + j);
                }

            // вставка в конец
            for (size_t i = 0; i < rout2.size(); ++i)
            {
                rout1.push_back(rout2[i]);
                rout2.erase(rout2.begin() + i);
                double lenght = length_rout(rout1);
                lenght += length_rout(rout2);
                if (best_lenght > lenght)
                {
                    best_lenght = lenght;
                    decision = { {0,i}, 4 };
                }
                rout2.emplace(rout2.begin() + i, rout1.back());
                rout1.pop_back();
            }
        }

        BSTM[p][q] = best_lenght - base_lenght;
        BSTM[q][p] = 0;
        RECM[p][q] = decision;
        RECM[q][p] = { {0, 0}, 6 };
    }

    void fill_BSTM_RECM(const vector<vector<size_t>>& routs, const checker_change* const checker)
    {
        BSTM = vector<vector<size_t>>(need_routs, vector < size_t>(need_routs, 0));
        RECM = vector<vector<pair<pair<size_t, size_t>, size_t>>>(need_routs,
            vector<pair<pair<size_t, size_t>, size_t>>(need_routs, { {0,0}, 6 }));
        for (size_t p = 0; p < need_routs; ++p)
            for (size_t q = p + 1; q < need_routs; ++q)
                checker->check_change(routs, p, q);
    }

    void recalculate_BSTM_RECM(const vector<vector<size_t>>& routs, size_t p, size_t q, const checker_change* const checker)
    {
        for (size_t i = 0; i < need_routs; ++i)
            if (i != p)
                checker->check_change(routs, p, i);

        for (size_t i = 0; i < need_routs; ++i)
            if (i != p && i != q)
                checker->check_change(routs, i, q);
    }

    pair<size_t, size_t> get_max_free_BSTM(size_t iter, const vector<vector<size_t>>& routs)
    {
        auto res = pair<size_t, size_t>(0, 0);
        double max = 0;
        for (size_t p = 0; p < need_routs; ++p)
            for (size_t q = p + 1; q < need_routs; ++q)
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


    vector<vector<size_t>> osman(vector<vector<size_t>> routs, const checker_change* const checker)
    {
        auto best_routs = routs;
        auto best_lenght = length_routs(routs);
        move_table = get_move_table();
        ts = max(7.0, 9.6 * log(count_point * need_routs) - 40);// 13
        size_t max_iter = 340 + 0.000353 * 5 * sqr(count_point * need_routs);
        fill_BSTM_RECM(routs, checker);

        for (size_t i = 0; i < max_iter; ++i)
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
                swap(routs[p][a], routs[q][b]);
                break;
            }
            case 1:
            {
                move_table[p][routs[p][a]] = i;
                routs[q].emplace(routs[q].begin() + b, routs[p][a]);
                routs[p].erase(routs[p].begin() + a);
                break;
            }
            case 2:
            {
                move_table[q][routs[q][b]] = i;
                routs[p].emplace(routs[p].begin() + a, routs[q][b]);
                routs[q].erase(routs[q].begin() + b);
                break;
            }
            case 3:
            {
                move_table[p][routs[p][a]] = i;
                routs[q].push_back(routs[p][a]);
                routs[p].erase(routs[p].begin() + a);
                break;
            }
            case 4:
            {
                move_table[q][routs[q][b]] = i;
                routs[p].push_back(routs[q][b]);
                routs[q].erase(routs[q].begin() + b);
                break;
            }
            default:
                break;
            }
            auto len = length_routs(routs);
            if (best_lenght > len)
            {
                best_lenght = len;
                best_routs = routs;
            }
            recalculate_BSTM_RECM(routs, p, q, checker);
        }
        return best_routs;
    }
}

namespace local_opt
{
    vector<size_t> TSP_2_opt(vector<size_t> rout)
    {
        auto best_len = length_rout(rout);
        auto best_rout = rout;
        for (int i = 0; i < rout.size() - 2; ++i)
            for (int j = i + 1; j < rout.size() - 1; ++j)
            {
                rout = vector<size_t>();
                for (int k = 0; k < i; ++k)
                    rout.push_back(best_rout[k]);
                for (int k = j; k >= i; --k)
                    rout.push_back(best_rout[k]);
                for (int k = j + 1; k < best_rout.size(); ++k)
                    rout.push_back(best_rout[k]);

                auto len = length_rout(rout);
                if (best_len > len)
                {
                    best_len = len;
                    best_rout = rout;
                }
            }
        return best_rout;
    }

    vector<size_t> TSP_3_opt(vector<size_t> rout)
    {
        auto best_len = length_rout(rout);
        auto best_rout = rout;
        // 3 ребра смежны
        for (size_t i = 0; i < rout.size() - 3; ++i)
        {
            swap(rout[i + 1], rout[i + 2]);

            auto len = length_rout(rout);
            if (best_len > len)
            {
                best_len = len;
                best_rout = rout;
            }
            else
                swap(rout[i + 1], rout[i + 2]);
        }

        // 2 ребра смежны, пусть вторые 2

        for (size_t i = 0; i < rout.size() - 4; ++i)
            for (size_t j = i + 2; j < rout.size() - 1; ++j)
            {
                vector<size_t> rout1;
                vector<size_t> rout2;

                for (size_t k = 0; k <= i; ++k)
                {
                    rout1.push_back(best_rout[k]);
                    rout2.push_back(best_rout[k]);
                }

                rout1.push_back(best_rout[j + 1]);
                for (size_t k = i + 1; k <= j; ++k)
                    rout1.push_back(best_rout[k]);
                for (size_t k = j + 2; k < best_rout.size(); ++k)
                    rout1.push_back(best_rout[k]);

                for (size_t k = j; k > i; --k)
                    rout2.push_back(best_rout[k]);
                for (size_t k = j + 1; k < best_rout.size(); ++k)
                    rout2.push_back(best_rout[k]);

                auto len = length_rout(rout1);
                if (best_len > len)
                {
                    best_len = len;
                    best_rout = rout1;
                }

                len = length_rout(rout2);
                if (best_len > len)
                {
                    best_len = len;
                    best_rout = rout2;
                }
            }

        // 2 ребра смежны, пусть первые 2

        for (size_t j = 0; j < rout.size() - 1; ++j)
            for (size_t i = j + 3; i < rout.size() - 4; ++i)
            {
                vector<size_t> rout1;
                vector<size_t> rout2;

                for (size_t k = 0; k <= j; ++k)
                {
                    rout1.push_back(best_rout[k]);
                    rout2.push_back(best_rout[k]);
                }


                for (size_t k = j + 2; k <= i; ++k)
                    rout1.push_back(best_rout[k]);
                rout1.push_back(best_rout[j + 1]);
                for (size_t k = i + 1; k < best_rout.size(); ++k)
                    rout1.push_back(best_rout[k]);

                for (size_t k = i; k > j; --k)
                    rout2.push_back(best_rout[k]);
                for (size_t k = i + 1; k < best_rout.size(); ++k)
                    rout2.push_back(best_rout[k]);

                auto len = length_rout(rout1);
                if (best_len > len)
                {
                    best_len = len;
                    best_rout = rout1;
                }

                len = length_rout(rout2);
                if (best_len > len)
                {
                    best_len = len;
                    best_rout = rout2;
                }
            }


        // 3 ребра не смежны

        for (size_t i = 0; i < rout.size() - 5; ++i)
            for (size_t j = i + 2; j < rout.size() - 3; ++j)
                for (size_t k = j + 2; i < rout.size() - 1; ++i)
                {
                    vector<size_t> rout1;
                    vector<size_t> rout2;

                    for (size_t l = 0; l <= i; ++l)
                    {
                        rout1.push_back(best_rout[l]);
                        rout2.push_back(best_rout[l]);
                    }

                    for (size_t l = k; l >= j + 1; --l)
                        rout1.push_back(best_rout[l]);
                    for (size_t l = i + 1; l <= j; ++l)
                        rout1.push_back(best_rout[l]);
                    for (size_t l = k + 1; l < best_rout.size(); ++l)
                        rout1.push_back(best_rout[l]);


                    for (size_t l = j; l >= i + 1; --l)
                        rout2.push_back(best_rout[l]);
                    for (size_t l = k; l >= j + 1; --l)
                        rout2.push_back(best_rout[l]);
                    for (size_t l = k + 1; l < best_rout.size(); ++l)
                        rout2.push_back(best_rout[l]);

                    auto len = length_rout(rout1);
                    if (best_len > len)
                    {
                        best_len = len;
                        best_rout = rout1;
                    }

                    len = length_rout(rout2);
                    if (best_len > len)
                    {
                        best_len = len;
                        best_rout = rout2;
                    }
                }
        return best_rout;
    }
}

namespace genetic
{
    double mutation_prob = 0.05;

    size_t rand_size_t(size_t max)
    {
        --max;
        return round(((double)rand()) * max / RAND_MAX);
    }

    struct gen
    {
        vector<size_t> rout;
        vector<size_t> sep;
    };


    double lenght_gen(const gen& gen)
    {
        size_t sep_in = 0;
        double lenght = dist_mat[0][gen.rout[0]] + dist_mat[gen.rout.back()][0];
        for (size_t i = 1; i < gen.rout.size(); ++i)
            if (sep_in < gen.sep.size() && i == gen.sep[sep_in])
            {
                lenght += dist_mat[0][gen.rout[i]] + dist_mat[gen.rout[i - 1]][0];;
                ++sep_in;
            }
            else
                lenght += dist_mat[gen.rout[i - 1]][gen.rout[i]];
        return lenght;
    }

    void check_repeat(const gen& gen)
    {
        auto h = vector<size_t>(51, false);
        for (auto t : gen.rout)
        {
            if (h[t])
                int ascac = 8;
            h[t] = true;
        }
    }

    gen get_child(const gen& parent1, const gen& parent2, size_t a, size_t b)
    {
        auto child = parent1;

        for (int i = a; i <= b; ++i)
        {
            auto miss_vertex = child.rout[i];
            child.rout[i] = parent2.rout[i];
            for (int j = 0; j < a; ++j)
                if (child.rout[j] == child.rout[i])
                    child.rout[j] = miss_vertex;
            for (int j = i + 1; j < parent1.rout.size(); ++j)
                if (child.rout[j] == child.rout[i])
                    child.rout[j] = miss_vertex;
        }

        auto i = rand_size_t(child.sep.size());
        auto j = rand_size_t(child.sep.size());

        child.sep[i] = rand_size_t(std::min(parent1.sep[i] + parent1.sep[j], count_point - 2)) + 1;
        sort(child.sep.begin(), child.sep.end());
        for (size_t k = 1; k < child.sep.size(); ++k)
            if (child.sep[k] == child.sep[k - 1])
            {
                child.sep = parent1.sep;
                break;
            }


        if ((double)rand() / RAND_MAX < mutation_prob)
        {
            auto i = rand_size_t(child.rout.size());
            auto j = rand_size_t(child.rout.size());
            swap(child.rout[i], child.rout[j]);
        }

        return child;
    }

    pair<gen, gen> crossover(const gen& parent1, const gen& parent2)
    {
        auto i = rand_size_t(parent1.rout.size());
        auto j = rand_size_t(parent1.rout.size());
        while (i == j)
            j = rand_size_t(parent1.rout.size());

        if (i > j)
            swap(i, j);
        auto child1 = get_child(parent1, parent2, i, j);

        i = rand_size_t(parent1.rout.size());
        j = rand_size_t(parent1.rout.size());
        while (i == j)
            j = rand_size_t(parent1.rout.size());

        if (i > j)
            swap(i, j);
        auto child2 = get_child(parent2, parent1, i, j);

        return { child1, child2 };
    }

    vector<size_t> generate_sep()
    {
        vector<size_t> seprators;
        for (size_t i = 0; i < need_routs - 1; ++i)
        {
            auto sep = rand_size_t(count_point - 2) + 1;
            bool is_break = false;
            for (size_t j = 0; j < seprators.size(); ++j)
                if (seprators[j] == sep)
                {
                    --i;
                    is_break = true;
                    continue;
                }
            if (is_break)
                continue;
            seprators.push_back(sep);
        }
        sort(seprators.begin(), seprators.end());
        return seprators;
    }

    vector<pair<double, gen>> generate_populatin(size_t size_populatin)
    {
        vector<pair<double, gen>> populatin;
        auto rng = std::default_random_engine{};

        for (size_t i = 0; i < size_populatin; ++i)
        {
            vector<size_t> vec(count_point - 1);
            for (size_t j = 1; j < count_point; ++j)
                vec[j - 1] = j;
            std::shuffle(vec.begin(), vec.end(), rng);
            auto sep = generate_sep();
            pair<double, gen> gen = { 0, { vec, sep } };
            gen.first = lenght_gen(gen.second);
            populatin.push_back(gen);
        }

        return populatin;
    }

    vector<vector<size_t>> routs_by_gen(const gen& gen)
    {
        size_t i = 0;
        auto routs = vector<vector<size_t>>(need_routs);
        for (size_t j = 0; j < gen.sep.size(); j++)
            for (; i < gen.sep[j]; ++i)
                routs[j].push_back(gen.rout[i]);

        for (; i < gen.rout.size(); ++i)
            routs[gen.sep.size()].push_back(gen.rout[i]);
        return routs;
    }


    vector<vector<size_t>> genetic_alg()
    {
        size_t size_populatin = 7000;
        size_t ages = 400;
        size_t size_update = 900;

        auto population = generate_populatin(size_populatin);

        for (size_t i = 0; i < ages; ++i)
        {
            for (size_t j = 0; j < size_update; ++j)
            {
                auto a = rand_size_t(size_populatin);
                auto b = rand_size_t(size_populatin);
                while (a == b)
                    b = rand_size_t(size_populatin);

                auto children = crossover(population[a].second, population[b].second);
                pair<double, gen> gen1 = { 0, children.first };
                gen1.first = lenght_gen(gen1.second);
                population.push_back(gen1);
                pair<double, gen> gen2 = { 0, children.second };
                gen2.first = lenght_gen(gen2.second);
                population.push_back(gen2);
            }

            sort(population.begin(), population.end(), [](const pair<double, gen>& a, const pair<double, gen>& b)
                {
                    return a.first < b.first;
                    //return a.second.sep.size() < b.second.sep.size();
                });
            for (size_t j = 0; j < size_update * 2; ++j)
                population.pop_back();
        }

        return routs_by_gen(population.front().second);
    }
}

int main()
{
    srand(time(0));
    fill_matrix();
    
    //auto routs_dichotomous_division = balancedVRP::clustering::dichotomous_division(dist_mat, 5);
    auto routs_dichotomous_division = balancedVRP::clustering::sweeping(x, y, 50, 5);
    //auto routs_osman = osman::osman(routs_dichotomous_division, new osman::checker_change_smoll_local());
    auto TSP_opt_rout = routs_dichotomous_division;
    for (size_t i = 0; i < TSP_opt_rout.size(); ++i)
    {
        TSP_opt_rout[i] = local_opt::TSP_2_opt(TSP_opt_rout[i]);
        TSP_opt_rout[i] = local_opt::TSP_3_opt(TSP_opt_rout[i]);
        TSP_opt_rout[i] = local_opt::TSP_2_opt(TSP_opt_rout[i]);
        TSP_opt_rout[i] = local_opt::TSP_3_opt(TSP_opt_rout[i]);
    }

    cout << "dichotomous_division: lenght= " << length_routs(TSP_opt_rout) << endl;
    for (const vector<size_t>& rout : TSP_opt_rout)
    {
        cout << "[" << 0 << ", ";
        for (auto vertex : rout)
            cout << vertex << ", ";
        cout << 0 << "]," << endl;
    }
    return 0;

    /*auto routs_genetic_alg = genetic::genetic_alg();
    cout << "genetic_alg: lenght= " << length_routs(routs_genetic_alg) << endl;
    for (const vector<size_t>& rout : routs_genetic_alg)
    {
        cout << "[" << 0 << ", ";
        for (auto vertex : rout)
            cout << vertex << ", ";
        cout << 0 << "]," << endl;
    }*/
    
    auto routs_clark_right = clark_right::clark_right();
    cout << "clark_right: lenght= " << length_routs(routs_clark_right) <<endl;
    for (const vector<size_t>& rout : routs_clark_right)
    {
        cout << "[" << 0 << ", ";
        for (auto vertex : rout)
            cout << vertex << ", ";
        cout << 0 << "]," << endl;
    }
    auto routs_osman = osman::osman(routs_clark_right, new osman::checker_change_smoll_local());
    cout << "osman: lenght= " << length_routs(routs_osman) << endl;
    for (const vector<size_t>& rout : routs_osman)
    {
        cout << "[" << 0 << ", ";
        for (auto vertex : rout)
            cout << vertex << ", ";
        cout << 0 << "]," << endl;
    }

    TSP_opt_rout = routs_osman;
    TSP_opt_rout[2] = local_opt::TSP_2_opt(TSP_opt_rout[2]);
    TSP_opt_rout[2] = local_opt::TSP_3_opt(TSP_opt_rout[2]);
    TSP_opt_rout[3] = local_opt::TSP_2_opt(TSP_opt_rout[3]);
    TSP_opt_rout[3] = local_opt::TSP_3_opt(TSP_opt_rout[3]);

    cout << "TSP_opt_rout: lenght= " << length_routs(TSP_opt_rout) << endl;
    for (const vector<size_t>& rout : TSP_opt_rout)
    {
        cout << "[" << 0 << ", ";
        for (auto vertex : rout)
            cout << vertex << ", ";
        cout << 0 << "]," << endl;
    }
}
