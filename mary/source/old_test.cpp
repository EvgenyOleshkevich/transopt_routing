#ifdef DEBUG_TEST_07_05_2021

void Ant_test_old()
{
    read_file();
    read_transports();
    dist_mat = utils::fill_matrix_with_end_point(x, y);

    Ant_algorithm ant(dist_mat, weights, transports);
    ant.run();
    auto mat = ant.res;
    cout << "length: " << ant.length() << endl;
    cout << "check: " << checkers::check_matrix(mat, dist_mat.size()) << endl;
    print::print_file(mat);
}

void osman_test_old()
{
    read_file();
    read_transports();
    auto P = utils::fill_matrix_and_sort(x, y);
    dist_mat = P.first;

    project::GreadyBase greedy(dist_mat, P.second, weights, transports);
    greedy.run();
    cout << "length greedy: " << greedy.length() << endl;
    cout << "check: " << checkers::check_matrix(greedy.res, dist_mat.size()) << endl;

    auto data = project::Osman::transform_data(greedy.res);

    project::Osman osman(dist_mat, weights, transports, data.first,
        data.second, data.first.size());

    osman.run();
    osman.add_zero_vertex();

    print::print(osman.res);
    cout << "length osman: " << osman.length() << endl;
    cout << "check: " << checkers::check_matrix(osman.res, dist_mat.size()) << endl;
}

void cluster_test1()
{
    dist_mat = utils::fill_matrix_with_end_point(x, y);

    auto clusters = clustering::dichotomous_division_weight(dist_mat, weights, clust_capacity, 7);

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


    auto vertex_x_clust = clustering::get_number_cluster_by_vertex(clusters);

    cout << "vertex_x_clust.size(): " << vertex_x_clust.size() << endl;
    //    sum volume : 4215.56
    //    frequence : 396
    //    acc_weight : 1225.62
}

void cluster_test2()
{
    dist_mat = utils::fill_matrix_with_end_point(x, y);
    auto clusters = clustering::dichotomous_division_weight(dist_mat, weights, clust_capacity, 7);

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


    auto clusters2 = project::clusters_with_multiple_duplicates(clusters, frequence);

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


#endif // DEBUG


#ifdef DEBUG_TEST_26_04_2021

const size_t count_point = 51;
const size_t need_routs = 5;
vector<double> x{ 0, 28, -7, 78, 23, -99, -96, -39, 45, -56, -2, -77, -25, 43, -72,
-31, 78, -49, -76, 23, -33, -22, -58, -79, 24, -6, -92, 41, -42, 92, -72, -25, -3,
74, 44, 46, -16, -2, 36, -21, -69, 19, 10, 39, 40, -84, 46, 15, -86, -65, 9 };
vector<double> y{ 0, 62, 80, 60, 83, 36, 62, -93, 95, -22, 42, 77, -34, -66, -43,
-69, 98, 94, -15, -34, -51, 3, -66, -69, 90, -53, 75, -28, 53, 46, -7, 44, 62,
-23, 48, -22, -90, -45, -50, 46, -18, 30, -28, 10, 43, 83, -84, 64, 45, 7, -63 };
vector<double> weights{ 0, 49, 68, 34, 27, 51, 3, 8, 76, 18, 76, 79, 41, 68,
79, 87, 14, 17, 39, 47, 30, 2, 56, 97, 37, 54, 39, 45, 88, 31, 65, 49, 54, 92, 8,
83, 31, 65, 80, 66, 40, 85, 10, 64, 40, 54, 86, 80, 63, 31, 24 };

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
    TSP::local_opt::opt_2(rout, dist_mat);
    TSP::local_opt::opt_3(rout, dist_mat);
    TSP::local_opt::opt_2(rout, dist_mat);
    TSP::local_opt::opt_3(rout, dist_mat);
    balancedVRP::VND_STS a(dist_mat, need_routs, 600, weights);

    std::vector<std::pair<size_t, double>> t;
    for (size_t v : rout)
    {
        t.push_back({ v, weights[v] });
    }

    auto routs = a.calculate(t).first;

    cout << "dynamic: length= " << utils::length_routs(routs, dist_mat) << endl;
    for (const vector<size_t>& rout : routs)
    {
        double w = 0;
        cout << "[" << 0 << ", ";
        for (auto vertex : rout)
        {
            w += weights[vertex];
            cout << vertex << ", ";
        }
        cout << 0 << "]," << endl;
    }
}


bool check(const vector<size_t> rout) {
    vector<size_t> rout_checker(rout.size() - 1, 0);
    bool check = rout[0] == 0 && rout[rout.size() - 1] == 0;
    for (size_t i = 1; i < rout.size() - 1; ++i) {
        ++rout_checker[rout[i]];
        check &= rout_checker[rout[i]] == 1;
        if (!check)
        {
            int c = i;
            int b = rout[i];
            int aa = rout_checker[rout[i]];
            int a = 8932;

        }
    }
    return check;
}

void test1()
{
    x.clear();
    y.clear();
    weights.clear();

    x.push_back(0);
    y.push_back(0);
    weights.push_back(0);
    std::ifstream in("file.csv");
    while (!in.eof())
    {
        string line;
        getline(in, line);
        double w, X, Y;
        in >> w >> X >> Y;
        x.push_back(X);
        y.push_back(Y);
        weights.push_back(w);
    }

    in.close();

    auto dist = utils::fill_matrix_with_end_point(x, y);

    auto clusters = balancedVRP::clustering::dichotomous_division(dist, 4);
    auto dist_matrixes = balancedVRP::clustering::get_dist_inner_cluster(dist, clusters);
    auto weights_by_clust = balancedVRP::clustering::get_weight_inner_cluster(weights, clusters);
    double length = 0;
    for (size_t i = 0; i < dist_matrixes.size(); ++i) {


        balancedVRP::VND_STS a(dist_matrixes[i], 60, 24, weights_by_clust[i]);

        vector<size_t> rout(clusters[i].size());
        for (size_t i = 0; i < rout.size(); ++i)
            rout[i] = i + 1;

        unsigned int start_time = clock(); // начальное время
        TSP::local_opt::opt_3_fast(rout, dist_matrixes[i]);
        unsigned int end_time = clock(); // конечное время
        unsigned int search_time = end_time - start_time; // искомое время
        cout << "time: " << search_time << endl;

        TSP::local_opt::opt_3(rout, dist_matrixes[i]);
        TSP::local_opt::opt_2(rout, dist_matrixes[i]);
        TSP::local_opt::opt_3(rout, dist_matrixes[i]);

        std::vector<std::pair<size_t, double>> t;
        for (size_t v : rout)
        {
            t.push_back({ v, weights_by_clust[i][v] });
        }

        auto routs = a.calculate(t).first;

        //auto routs = algorithms::clark_right::clark_right(dist_matrixes[i], dist_matrixes[i].size(), 60);
        // auto osman = new Osman(dist_matrixes[i], dist_matrixes[i].size(), 60);
        length += utils::length_routs(routs, dist_matrixes[i]);
    }
    cout << "length: " << length << endl;
}

void opt_2_test()
{
    /*
        length: 79203.3
        check: 1

        opt_2_symmetrical time: 1915
        opt_2_symmetrical check: 1
        length: 10144.1

        opt_2_fast2 time: 2266
        opt_2_fast2 check: 1
        length: 10144.1

        opt_2 time: 64725
        opt_2 check: 1
        length: 11501.9

        opt_2_fast time: 23676
        opt_2_fast check: 1
        length: 11501.9
    */
    x.clear();
    y.clear();
    weights.clear();

    x.push_back(0);
    y.push_back(0);
    weights.push_back(0);
    std::ifstream in("file.csv");
    while (!in.eof())
    {
        string line;
        getline(in, line);
        double w, X, Y;
        in >> w >> X >> Y;
        x.push_back(X);
        y.push_back(Y);
        weights.push_back(w);
    }

    in.close();

    auto dist = utils::fill_matrix_with_end_point(x, y);

    auto clusters = balancedVRP::clustering::dichotomous_division(dist, 4);
    auto dist_matrixes = balancedVRP::clustering::get_dist_inner_cluster(dist, clusters);
    auto weights_by_clust = balancedVRP::clustering::get_weight_inner_cluster(weights, clusters);
    double length = 0;


    vector<size_t> rout(clusters[0].size() + 2);

    rout[0] = 0;
    rout[rout.size() - 1] = 0;
    for (size_t i = 1; i < rout.size() - 1; ++i)
        rout[i] = i;
    cout << "length: " << utils::length_rout(rout, dist_matrixes[0]) << endl;
    cout << "check: " << check(rout) << endl;

    // opt_2_symmetrical
    {
        auto rout_test = rout;
        unsigned int start_time = clock(); // начальное время
        for (size_t i = 0; i < 10; ++i)
            TSP::local_opt::opt_2_symmetrical(rout_test, dist_matrixes[0]);
        unsigned int end_time = clock(); // конечное время
        unsigned int search_time = end_time - start_time; // искомое время
        cout << "opt_2_symmetrical time: " << search_time << endl;
        cout << "opt_2_symmetrical check: " << check(rout_test) << endl;
        cout << "length: " << utils::length_rout(rout_test, dist_matrixes[0]) << endl;
    }

    // opt_2_fast2
    {
        auto rout_test = rout;
        unsigned int start_time = clock(); // начальное время
        for (size_t i = 0; i < 10; ++i)
            TSP::local_opt::opt_2_fast2(rout_test, dist_matrixes[0]);
        unsigned int end_time = clock(); // конечное время
        unsigned int search_time = end_time - start_time; // искомое время
        cout << "opt_2_fast2 time: " << search_time << endl;
        cout << "opt_2_fast2 check: " << check(rout_test) << endl;
        cout << "length: " << utils::length_rout(rout_test, dist_matrixes[0]) << endl;
    }

    // opt_2
    {
        auto rout_test = rout;
        unsigned int start_time = clock(); // начальное время
        TSP::local_opt::opt_2(rout_test, dist_matrixes[0]);
        unsigned int end_time = clock(); // конечное время
        unsigned int search_time = end_time - start_time; // искомое время
        cout << "opt_2 time: " << search_time << endl;
        cout << "opt_2 check: " << check(rout_test) << endl;
        cout << "length: " << utils::length_rout(rout_test, dist_matrixes[0]) << endl;
    }

    // opt_2_fast
    {
        auto rout_test = rout;
        unsigned int start_time = clock(); // начальное время
        TSP::local_opt::opt_2_fast(rout_test, dist_matrixes[0]);
        unsigned int end_time = clock(); // конечное время
        unsigned int search_time = end_time - start_time; // искомое время
        cout << "opt_2_fast time: " << search_time << endl;
        cout << "opt_2_fast check: " << check(rout_test) << endl;
        cout << "length: " << utils::length_rout(rout_test, dist_matrixes[0]) << endl;
    }


}

void opt_2_test2()
{
    /*
            opt_2_best_parallel time: 21216 -> 7693
        opt_2_best_parallel check: 1
        length: 1.11833e+06

        opt_2_best time: 35940 -> 14946
        opt_2_best check: 1
        length: 1.11536e+06
    */
    x.clear();
    y.clear();
    weights.clear();

    x.push_back(0);
    y.push_back(0);
    weights.push_back(0);
    std::ifstream in("file.csv");
    while (!in.eof())
    {
        string line;
        getline(in, line);
        double w, X, Y;
        in >> w >> X >> Y;
        x.push_back(X);
        y.push_back(Y);
        weights.push_back(w);
    }

    in.close();

    auto dist = utils::fill_matrix_with_end_point(x, y);


    double length = 0;


    vector<size_t> rout(dist.size() + 1);

    rout[0] = 0;
    rout[rout.size() - 1] = 0;
    for (size_t i = 1; i < rout.size() - 1; ++i)
        rout[i] = i;
    cout << "length: " << utils::length_rout(rout, dist) << endl;
    cout << "check: " << check(rout) << endl;

    // opt_2_symmetrical
    //{
    //    auto rout_test = rout;
    //    unsigned int start_time = clock(); // начальное время
    //    for (size_t i = 0; i < 2; ++i)
    //        TSP::local_opt::opt_2_symmetrical(rout_test, dist);
    //    unsigned int end_time = clock(); // конечное время
    //    unsigned int search_time = end_time - start_time; // искомое время
    //    cout << endl << "opt_2_symmetrical time: " << search_time << endl;
    //    cout << "opt_2_symmetrical check: " << check(rout_test) << endl;
    //    cout << "length: " << utils::length_rout(rout_test, dist) << endl;
    //}

        // opt_2_best
    {
        auto rout_test = rout;
        unsigned int start_time = clock(); // начальное время

        TSP::local_opt::opt_2_best_parallel(rout_test, dist);
        unsigned int end_time = clock(); // конечное время
        unsigned int search_time = end_time - start_time; // искомое время
        cout << endl << "opt_2_best_parallel time: " << search_time << endl;
        cout << "opt_2_best_parallel check: " << check(rout_test) << endl;
        cout << "length: " << utils::length_rout(rout_test, dist) << endl;
    }

    // opt_2_best
    {
        auto rout_test = rout;
        unsigned int start_time = clock(); // начальное время
        TSP::local_opt::opt_2_best_step(rout_test, dist);
        unsigned int end_time = clock(); // конечное время
        unsigned int search_time = end_time - start_time; // искомое время
        cout << endl << "opt_2_best time: " << search_time << endl;
        cout << "opt_2_best check: " << check(rout_test) << endl;
        cout << "length: " << utils::length_rout(rout_test, dist) << endl;
    }


}


int old_main()
{
    /*auto t1 = std::thread(&out);
    auto t2 = std::thread(&out2);
    t1.join();
    t2.join();
    out();
    out2();*/
    unsigned int start_time = clock(); // начальное время
    opt_2_test2();
    unsigned int end_time = clock(); // конечное время
    unsigned int search_time = end_time - start_time; // искомое время
    cout << "time: " << search_time;
    return 0;
    srand((unsigned int)time(0));
    dist_mat = utils::fill_matrix(x, y, count_point);

    auto routs_dichotomous_division = balancedVRP::clustering::dichotomous_division_weight(dist_mat, weights, 600, 5);
    //auto routs_osman = osman::osman(routs_dichotomous_division, new osman::checker_change_smoll_local());
    auto TSP_opt_rout = routs_dichotomous_division;
    for (size_t i = 0; i < TSP_opt_rout.size(); ++i)
    {
        TSP::local_opt::opt_2(TSP_opt_rout[i], dist_mat);
        TSP::local_opt::opt_3(TSP_opt_rout[i], dist_mat);
        TSP::local_opt::opt_2(TSP_opt_rout[i], dist_mat);
        TSP::local_opt::opt_3(TSP_opt_rout[i], dist_mat);
    }

    cout << "dichotomous_division: length= " << utils::length_routs(TSP_opt_rout, dist_mat) << endl;
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
    std::cout << "Lin_Kernighan: length= " <<
        utils::length_rout(rout, dist_mat) << std::endl;

    auto cutting_routs = balancedVRP::cutting_rout(rout, dist_mat, 6);

    for (const vector<size_t>& rout : cutting_routs)
    {
        cout << "[" << 0 << ", ";
        for (auto vertex : rout)
            cout << vertex << ", ";
        cout << 0 << "]," << endl;
    }
    cout << "cutting_rout: length= " << utils::length_routs(cutting_routs, dist_mat) << endl;

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

    cout << "dichotomous_division: length= " << utils::length_routs(TSP_opt_rout, dist_mat) << endl;
    for (const vector<size_t>& rout : TSP_opt_rout)
    {
        cout << "[" << 0 << ", ";
        for (auto vertex : rout)
            cout << vertex << ", ";
        cout << 0 << "]," << endl;
    }

    auto routs_genetic_alg = genetic::genetic_alg();
    cout << "genetic_alg: length= " << length_routs(routs_genetic_alg) << endl;
    for (const vector<size_t>& rout : routs_genetic_alg)
    {
        cout << "[" << 0 << ", ";
        for (auto vertex : rout)
            cout << vertex << ", ";
        cout << 0 << "]," << endl;
    }

    auto routs_clark_right = algorithms::clark_right::clark_right(dist_mat, count_point, 5);
    cout << "clark_right: length= " << utils::length_routs(routs_clark_right, dist_mat) <<endl;
    for (const vector<size_t>& rout : routs_clark_right)
    {
        cout << "[" << 0 << ", ";
        for (auto vertex : rout)
            cout << vertex << ", ";
        cout << 0 << "]," << endl;
    }
    algorithms::Osman osman(dist_mat, count_point, 5);
    auto routs_osman = osman.start(routs_clark_right, 1);
    cout << "osman: length= " << utils::length_routs(routs_osman, dist_mat) << endl;
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

    cout << "TSP_opt_rout: length= " << utils::length_routs(TSP_opt_rout, dist_mat) << endl;
    for (const vector<size_t>& rout : TSP_opt_rout)
    {
        cout << "[" << 0 << ", ";
        for (auto vertex : rout)
            cout << vertex << ", ";
        cout << 0 << "]," << endl;
    }*/
}

#endif // DEBUG