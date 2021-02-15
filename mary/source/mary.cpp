#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <random>
using namespace std;

namespace mary
{
    /* enum class - this is a C++ 11 feature */
    enum class FieldType {
        X = 1,
        O = -1,
        Empty = 0
    };

    constexpr size_t GAME_SIZE = 3;
    using GAME_BOARD = vector<vector<FieldType>>;

    using Turn = pair<unsigned short, unsigned short>;

    class BotVasiliy;
    class Me;
    class Manager;

    class GameBoard {
    public:
        GameBoard() : desk(3, vector<FieldType>(GAME_SIZE)) {
            clean(desk);
        }

        const GAME_BOARD& getFieldsStates() const {
            return desk;
        }

        static char getSymbol(const FieldType& state) {
            if (state == FieldType::Empty) return '-';
            if (state == FieldType::X) return 'x';
            if (state == FieldType::O) return 'o';

            return '#';
        }

    private:
        GAME_BOARD desk;
        friend ostream& operator << (ostream& out, const GameBoard& gameBoard) {
            /*
            o|o|-
            x|x|-
            -|-|o
            */

            //GAME_BOARD* toPrint = gameBoard.getFieldsStates();

            for (int i = 0; i < GAME_SIZE; i++) {
                for (int j = 0; j < GAME_SIZE; j++) {
                    out << GameBoard::getSymbol(gameBoard.desk[i][j]);
                    if (j != (GAME_SIZE - 1)) out << "|";
                }

                out << std::endl;
            }

            return out;
        }
        friend BotVasiliy;
        friend Me;
        friend Manager;

        static void clean(GAME_BOARD& desk)
        {
            for (int i = 0; i < GAME_SIZE; i++)
                for (int j = 0; j < GAME_SIZE; j++)
                    desk[i][j] = FieldType::Empty;
        }
    };

    class Player {

    protected:
        FieldType me;
        FieldType oponent;
        string name;
        friend GameBoard;
    public:
        virtual FieldType getPlayerType() {
            return me;
        }

        virtual const string& getPlayerName()
        {
            return name;
        }
        virtual Turn makeTurn(GameBoard& board) = 0;

        Player(FieldType playerType, string name) : me(playerType), name(name)
        {
            if (FieldType::X == me)
                oponent = FieldType::O;
            else
                oponent = FieldType::X;
        }
    };

    class BotVasiliy : public Player
    {
    public:
        BotVasiliy(FieldType type, string name) : Player(type, name)
        {}

        Turn makeTurn(GameBoard& board)
        {
            Turn turn;
            if (checkWin(board, turn))
                return turn;
            if (checkDanger(board, turn))
                return turn;

            if (board.desk[1][1] == FieldType::Empty)
            {
                board.desk[1][1] = me;
                return { 1, 1 };
            }

            auto count_op = countOponent(board);
            if (count_op == 1)
            {
                singleOponent(board, turn);
                return turn;
            }
            if (me == FieldType::O && count_op == 2 && checkDoubleOponent(board, turn))
            {
                return turn;
            }
            another(board, turn);
            return turn;
        }

    private:
        bool checkWin(GameBoard& board, Turn& turn)
        {
            int count = 0;
            for (size_t i = 0; i < GAME_SIZE; ++i)
            {
                count = 0;
                for (size_t j = 0; j < GAME_SIZE; ++j)
                {
                    if (oponent == board.desk[i][j])
                        --count;
                    else if (me == board.desk[i][j])
                        ++count;
                    else
                        turn = { i, j };
                }
                if (count == 2)
                {
                    board.desk[turn.first][turn.second] = me;
                    return true;
                }
            }
            for (size_t i = 0; i < GAME_SIZE; ++i)
            {
                count = 0;
                for (size_t j = 0; j < GAME_SIZE; ++j)
                {
                    if (oponent == board.desk[j][i])
                        --count;
                    else if (me == board.desk[j][i])
                        ++count;
                    else
                        turn = { j, i };
                }
                if (count == 2)
                {
                    board.desk[turn.first][turn.second] = me;
                    return true;
                }
            }
            count = 0;
            for (size_t i = 0; i < GAME_SIZE; ++i)
            {
                if (oponent == board.desk[i][i])
                    --count;
                else if (me == board.desk[i][i])
                    ++count;
                else
                    turn = { i, i };
            }
            if (count == 2)
            {
                board.desk[turn.first][turn.second] = me;
                return true;
            }
            count = 0;
            for (size_t i = 0; i < GAME_SIZE; ++i)
            {
                if (oponent == board.desk[i][GAME_SIZE - i - 1])
                    --count;
                else if (me == board.desk[i][GAME_SIZE - i - 1])
                    ++count;
                else
                    turn = { i, GAME_SIZE - i - 1 };
            }
            if (count == 2)
            {
                board.desk[turn.first][turn.second] = me;
                return true;
            }
            return false;
        }

        bool checkDanger(GameBoard& board, Turn& turn)
        {
            int count = 0;
            for (size_t i = 0; i < GAME_SIZE; ++i)
            {
                count = 0;
                for (size_t j = 0; j < GAME_SIZE; ++j)
                {
                    if (oponent == board.desk[i][j])
                        ++count;
                    else if (me == board.desk[i][j])
                        --count;
                    else
                        turn = { i, j };
                }
                if (count == 2)
                {
                    board.desk[turn.first][turn.second] = me;
                    return true;
                }
            }
            for (size_t i = 0; i < GAME_SIZE; ++i)
            {
                count = 0;
                for (size_t j = 0; j < GAME_SIZE; ++j)
                {
                    if (oponent == board.desk[j][i])
                        ++count;
                    else if (me == board.desk[j][i])
                        --count;
                    else
                        turn = { j, i };
                }
                if (count == 2)
                {
                    board.desk[turn.first][turn.second] = me;
                    return true;
                }
            }
            count = 0;
            for (size_t i = 0; i < GAME_SIZE; ++i)
            {
                if (oponent == board.desk[i][i])
                    ++count;
                else if (me == board.desk[i][i])
                    --count;
                else
                    turn = { i, i };
            }
            if (count == 2)
            {
                board.desk[turn.first][turn.second] = me;
                return true;
            }
            count = 0;
            for (size_t i = 0; i < GAME_SIZE; ++i)
            {
                if (oponent == board.desk[i][GAME_SIZE - i - 1])
                    ++count;
                else if (me == board.desk[i][GAME_SIZE - i - 1])
                    --count;
                else
                    turn = { i, GAME_SIZE - i - 1 };
            }
            if (count == 2)
            {
                board.desk[turn.first][turn.second] = me;
                return true;
            }
            return false;
        }

        void singleOponent(GameBoard& board, Turn& turn)
        {
            if (board.desk[1][1] == oponent)
            {
                turn = { 0, 0 };
                board.desk[0][0] = me;
                return;
            }
            if (board.desk[0][0] == oponent)
            {
                turn = { 2, 2 };
                board.desk[2][2] = me;
                return;
            }
            if (board.desk[0][2] == oponent)
            {
                turn = { 2, 0 };
                board.desk[2][0] = me;
                return;
            }
            if (board.desk[2][0] == oponent)
            {
                turn = { 0, 2 };
                board.desk[0][2] = me;
                return;
            }
            if (board.desk[2][2] == oponent)
            {
                turn = { 0, 0 };
                board.desk[0][0] = me;
                return;
            }
            if (board.desk[0][1] == oponent)
            {
                turn = { 2, 0 };
                board.desk[2][0] = me;
                return;
            }
            if (board.desk[1][2] == oponent)
            {
                turn = { 2, 0 };
                board.desk[2][0] = me;
                return;
            }
            if (board.desk[1][0] == oponent)
            {
                turn = { 0, 2 };
                board.desk[0][2] = me;
                return;
            }
            if (board.desk[2][1] == oponent)
            {
                turn = { 0, 2 };
                board.desk[0][2] = me;
                return;
            }

        }

        void doubleOponent(GameBoard& board, Turn& turn)
        {
            if (board.desk[0][0] == FieldType::Empty)
            {
                turn = { 0, 0 };
                board.desk[0][0] = me;
                return;
            }
            if (board.desk[2][0] == FieldType::Empty)
            {
                turn = { 2, 0 };
                board.desk[2][0] = me;
                return;
            }
            if (board.desk[0][2] == FieldType::Empty)
            {
                turn = { 0, 2 };
                board.desk[0][2] = me;
                return;
            }
            if (board.desk[2][2] == FieldType::Empty)
            {
                turn = { 2, 2 };
                board.desk[2][2] = me;
                return;
            }
        }

        bool checkDoubleOponent(GameBoard& board, Turn& turn)
        {
            if (board.desk[0][1] == oponent && board.desk[1][2] == oponent)
            {
                turn = { 0, 2 };
                board.desk[0][2] = me;
                return true;
            }
            if (board.desk[1][2] == oponent && board.desk[2][1] == oponent)
            {
                turn = { 2, 2 };
                board.desk[2][2] = me;
                return true;
            }
            if (board.desk[2][1] == oponent && board.desk[1][0] == oponent)
            {
                turn = { 2, 0 };
                board.desk[2][0] = me;
                return true;
            }
            if (board.desk[1][0] == oponent && board.desk[0][1] == oponent)
            {
                turn = { 0, 0 };
                board.desk[0][0] = me;
                return true;
            }
            if (board.desk[1][1] == oponent)
            {
                if (board.desk[0][0] == oponent && board.desk[2][2] == me)
                {
                    turn = { 0, 2 };
                    board.desk[0][2] = me;
                    return true;
                }
                if (board.desk[2][2] == oponent && board.desk[0][0] == me)
                {
                    turn = { 2, 0 };
                    board.desk[2][0] = me;
                    return true;
                }
                if (board.desk[2][0] == oponent && board.desk[0][2] == me)
                {
                    turn = { 0, 0 };
                    board.desk[0][0] = me;
                    return true;
                }
                if (board.desk[0][2] == oponent && board.desk[2][0] == me)
                {
                    turn = { 0, 0 };
                    board.desk[0][0] = me;
                    return true;
                }
            }
            return false;
        }

        void another(GameBoard& board, Turn& turn)
        {
            vector<vector<size_t>> MVP(3, vector<size_t>(3, 0)); // most valuable point
            int count = 1;
            for (size_t i = 0; i < GAME_SIZE; ++i)
            {
                count = 1;
                for (size_t j = 0; j < GAME_SIZE; ++j)
                {
                    if (oponent == board.desk[i][j])
                        count = -10;
                    else if (me == board.desk[i][j])
                        ++count;
                }
                if (count > 0)
                {
                    for (size_t j = 0; j < GAME_SIZE; ++j)
                        if (board.desk[i][j] == FieldType::Empty)
                            MVP[i][j] += count;
                }
            }
            for (size_t i = 0; i < GAME_SIZE; ++i)
            {
                count = 1;
                for (size_t j = 0; j < GAME_SIZE; ++j)
                {
                    if (oponent == board.desk[j][i])
                        count = -10;
                    else if (me == board.desk[j][i])
                        ++count;
                }
                if (count > 0)
                {
                    for (size_t j = 0; j < GAME_SIZE; ++j)
                        if (board.desk[i][j] == FieldType::Empty)
                            MVP[j][i] += count;
                }
            }
            count = 1;
            for (size_t i = 0; i < GAME_SIZE; ++i)
            {
                if (oponent == board.desk[i][i])
                    --count;
                else if (me == board.desk[i][i])
                    ++count;
                else
                    turn = { i, i };
            }
            if (count > 0)
            {
                for (size_t i = 0; i < GAME_SIZE; ++i)
                    if (board.desk[i][i] == FieldType::Empty)
                        MVP[i][i] += count;
            }
            count = 0;
            for (size_t i = 0; i < GAME_SIZE; ++i)
            {
                if (oponent == board.desk[i][GAME_SIZE - i - 1])
                    --count;
                else if (me == board.desk[i][GAME_SIZE - i - 1])
                    ++count;
                else
                    turn = { i, GAME_SIZE - i - 1 };
            }
            if (count > 0)
            {
                for (size_t i = 0; i < GAME_SIZE; ++i)
                    if (board.desk[i][GAME_SIZE - i - 1] == FieldType::Empty)
                        MVP[i][GAME_SIZE - i - 1] += count;
            }
            for (size_t i = 0; i < GAME_SIZE; ++i)
                for (size_t j = 0; j < GAME_SIZE; ++j)
                    if (board.desk[i][j] == FieldType::Empty)
                    {
                        turn = { i, j };
                        i = GAME_SIZE; // for break outer cicle
                        break;
                    }
            for (size_t i = 0; i < GAME_SIZE; ++i)
                for (size_t j = 0; j < GAME_SIZE; ++j)
                    if (MVP[turn.first][turn.second] < MVP[i][j])
                        turn = { i, j };
            board.desk[turn.first][turn.second] = me;
        }

        size_t countOponent(GameBoard& board)
        {
            size_t count = 0;
            for (size_t i = 0; i < GAME_SIZE; ++i)
                for (size_t j = 0; j < GAME_SIZE; ++j)
                    if (oponent == board.desk[i][j])
                        ++count;
            return count;
        }
    };

    class Me : public Player
    {
    public:
        Me(FieldType type, string name) : Player(type, name)
        {}

        Turn makeTurn(GameBoard& board)
        {
            size_t x, y;
            do
            {
                system("cls");
                cout << board << "write your turn correct:  ";
                cin >> x >> y;
                --x;
                --y;

            } while (x > 2 || y > 2 || board.desk[x][y] != FieldType::Empty);
            board.desk[x][y] = me;
            return { x, y };
        }
    };

    class Manager
    {
    public:
        void Play()
        {
            char x;
            cout << "write 'x' if you want to play: ";
            cin >> x;
            while (x == 'x')
            {
                Start();
                system("cls");
                Player* winner = nullptr;
                while (!winner)
                {
                    first->makeTurn(board);
                    winner = is_not_finished();
                    if (!winner)
                        second->makeTurn(board);
                    else
                        break;
                    winner = is_not_finished();
                }
                system("cls");
                cout << board;
                cout << winner->getPlayerName() << " won" << endl;
                cout << "write 'x' if you want to play: ";
                if (winner != first && winner != second)
                    delete winner;
                delete first;
                delete second;
                cin >> x;
            }
        }

    private:
        GameBoard board;
        Player* first;
        Player* second;

        void Start()
        {
            board.clean(board.desk);
            cout << "write 'x' if you want to play for it or any char: ";
            char x;
            cin >> x;
            if (x == 'x')
            {
                first = new Me(FieldType::X, "Me");
                second = new BotVasiliy(FieldType::O, "BotVasiliy");
            }
            else
            {
                first = new BotVasiliy(FieldType::X, "BotVasiliy");
                second = new Me(FieldType::O, "Me");
            }
        }

        Player* is_not_finished()
        {
            int count = 0;
            for (size_t i = 0; i < GAME_SIZE; ++i)
            {
                count = 0;
                for (size_t j = 0; j < GAME_SIZE; ++j)
                {
                    if (first->getPlayerType() == board.desk[i][j])
                        --count;
                    else if (second->getPlayerType() == board.desk[i][j])
                        ++count;
                }
                if (count == -3)
                    return first;
                if (count == 3)
                    return second;
            }
            for (size_t i = 0; i < GAME_SIZE; ++i)
            {
                count = 0;
                for (size_t j = 0; j < GAME_SIZE; ++j)
                {
                    if (first->getPlayerType() == board.desk[j][i])
                        --count;
                    else if (second->getPlayerType() == board.desk[j][i])
                        ++count;
                }
                if (count == -3)
                    return first;
                if (count == 3)
                    return second;
            }
            count = 0;
            for (size_t i = 0; i < GAME_SIZE; ++i)
            {
                if (first->getPlayerType() == board.desk[i][i])
                    --count;
                else if (second->getPlayerType() == board.desk[i][i])
                    ++count;
            }
            if (count == -3)
                return first;
            if (count == 3)
                return second;
            count = 0;
            for (size_t i = 0; i < GAME_SIZE; ++i)
            {
                if (first->getPlayerType() == board.desk[i][GAME_SIZE - i - 1])
                    --count;
                else if (second->getPlayerType() == board.desk[i][GAME_SIZE - i - 1])
                    ++count;
            }
            if (count == -3)
                return first;
            if (count == 3)
                return second;

            count = 0;
            for (size_t i = 0; i < GAME_SIZE; ++i)
                for (size_t j = 0; j < GAME_SIZE; ++j)
                    if (FieldType::Empty == board.desk[i][j])
                        ++count;
            if (count == 0)
                return new Me(FieldType::Empty, "draw, friend");
            return nullptr;
        }
    };


    class My_Player
    {
    public:
        My_Player(size_t _in_sug, size_t _out_sug) :
            in_sug(_in_sug), out_sug(_out_sug), money(0) {}

        void offer_to(My_Player& player)
        {
            if (player.accept_offer(100 - out_sug))
                money += out_sug;
        }

        size_t Money() const
        {
            return money;
        }
        size_t In() const
        {
            return in_sug;
        }
        size_t Out() const
        {
            return out_sug;
        }

        void Reset_money()
        {
            money = 0;
        }

        const bool operator < (const My_Player& player) const
        {
            return money > player.money;
        }
        const bool operator > (const My_Player& player) const
        {
            return money < player.money;
        }
    private:
        bool accept_offer(size_t x)
        {
            if (x >= in_sug)
            {
                money += x;
                return true;
            }
            else
                return false;
        }
        size_t in_sug;
        size_t out_sug;
        size_t money;
    };

    /*int main()
    {
        using std::vector;
        vector<My_Player> vec;
        for (size_t i = 0; i < 60; ++i)
            for (size_t j = 0; j < 90; ++j)
                vec.push_back(My_Player(i + 1, j + 1));
        for (size_t i = 0; i < vec.size(); ++i)
            for (size_t j = 0; j < vec.size(); ++j)
                vec[i].offer_to(vec[j]);

        std::sort(vec.begin(), vec.end());

        int u = 0;
        for (size_t i = 0; i < vec.size(); ++i)
                vec[i].Reset_money();
        for (size_t i = 0; i < 4000; ++i)
            for (size_t j = 0; j < 4000; ++j)
                vec[i].offer_to(vec[j]);
        std::sort(vec.begin(), vec.end());
        u = 0;

        for (size_t i = 0; i < vec.size(); ++i)
            vec[i].Reset_money();
        for (size_t i = 0; i < 3000; ++i)
            for (size_t j = 0; j < 3000; ++j)
                vec[i].offer_to(vec[j]);
        std::sort(vec.begin(), vec.end());
        u = 0;

        for (size_t i = 0; i < vec.size(); ++i)
            vec[i].Reset_money();
        for (size_t i = 0; i < 2000; ++i)
            for (size_t j = 0; j < 2000; ++j)
                vec[i].offer_to(vec[j]);
        std::sort(vec.begin(), vec.end());
        u = 0;

        for (size_t i = 0; i < vec.size(); ++i)
            vec[i].Reset_money();
        for (size_t i = 0; i < 1000; ++i)
            for (size_t j = 0; j < 1000; ++j)
                vec[i].offer_to(vec[j]);
        std::sort(vec.begin(), vec.end());
        u = 0;

        for (size_t i = 0; i < vec.size(); ++i)
            vec[i].Reset_money();
        for (size_t i = 0; i < 500; ++i)
            for (size_t j = 0; j < 500; ++j)
                vec[i].offer_to(vec[j]);
        std::sort(vec.begin(), vec.end());
        u = 0;

        for (size_t i = 0; i < vec.size(); ++i)
            vec[i].Reset_money();
        for (size_t i = 0; i < 200; ++i)
            for (size_t j = 0; j < 200; ++j)
                vec[i].offer_to(vec[j]);
        std::sort(vec.begin(), vec.end());
        u = 0;

        for (size_t i = 0; i < vec.size(); ++i)
            vec[i].Reset_money();
        for (size_t i = 0; i < 100; ++i)
            for (size_t j = 0; j < 100; ++j)
                vec[i].offer_to(vec[j]);
        std::sort(vec.begin(), vec.end());
        u = 0;

        for (size_t i = 0; i < vec.size(); ++i)
            vec[i].Reset_money();
        for (size_t i = 0; i < 50; ++i)
            for (size_t j = 0; j < 50; ++j)
                vec[i].offer_to(vec[j]);
        std::sort(vec.begin(), vec.end());
        u = 0;



        for (size_t i = 0; i < 1000; ++i)
            if (vec[i].In() > u)
                u = vec[i].In();
        int y = 0;
        return 0;
        GameBoard game1;
        std::cout << "Hello from AWS Cloud9!" << std::endl;
        std::cout << game1 << std::endl;

        GAME_BOARD board1 = game1.getFieldsStates();
        std::cout << "At position [1][1]: " << GameBoard::getSymbol(board1[1][1]) << std::endl;
        Manager m;
        m.Play();
    }*/
}