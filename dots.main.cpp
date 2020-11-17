#include <bits/stdc++.h>

#define C0  "\033[0m"
#define C1  "\033[31m"
#define C2  "\033[34m"
#define C1B "\033[1;31m"
#define C2B "\033[1;34m"

#ifdef PROUT
#define PR(x) x
#else
#define PR(x)
#endif

/* Mask for the color bits */
#define C_MASK   3U

/* Bit set if the dot is present */
#define B_SET    (C_MASK + 1)

#define RED_DOT  (1 | B_SET)
#define BLUE_DOT (2 | B_SET)

/* How many bits in color + presence */
#define DOT_BITS 3

/* Mask for color + presence */
#define DOT_MASK ((1 << DOT_BITS) - 1)


#define CHUNK_SIZE 1024

#ifdef NN
#define N NN
#else
#define N 30
#endif


inline static uint64_t tsc() {
    uint64_t a, d;
    asm volatile("rdtsc" : "=a"(a), "=d"(d));
    return (d<<32) + a;
}


using namespace std;

string inFile, outFile, dotFile;

struct MCTS_Node {
    explicit MCTS_Node()
        : parent(0), next_neighb(0), next_level(0), wins(0), games(0), x(0), y(0), color(0)
    {}

    MCTS_Node *parent, *next_neighb, *next_level;
    double wins;
    uint64_t games;
    uint8_t x, y, color;

    MCTS_Node(MCTS_Node const&) = delete;
};

class MCTS_Tree {
public:
    explicit MCTS_Tree()
        : next_free(0), n(CHUNK_SIZE), n_chunks(0)
    {}

    MCTS_Node* new_node(MCTS_Node * p, int x = 0, int y = 0) {
        if (n++ == CHUNK_SIZE) {
            if (n_chunks > 100000) {
                throw new runtime_error("Number of MCTS_Tree chunks exceeded");
            }
            next_free = new MCTS_Node[CHUNK_SIZE];
            n = 1;
            ++n_chunks;
        }
        next_free->parent = p;
        if (p) {
            next_free->next_neighb = p->next_level;
            p->next_level = next_free;
            next_free->color = (p->color == 1) ? 2 : 1;
        } else {
            next_free->color = 2;
        }
        next_free->x = x;
        next_free->y = y;
        return next_free++;
    }

    void save(MCTS_Node *root, const string& fileName) {
        ofstream file;
        file.open(fileName, ios::binary | ios::trunc);
        saveNode(0, root, file);
        file.close();
    }

    void saveNode(int level, MCTS_Node *node, ofstream& file) {
        file.write((char*)&level, sizeof(level));
        file.write((char*)&node->wins, sizeof(node->wins)); 
        file.write((char*)&node->games, sizeof(node->games));
        file.write((char*)&node->x, 1);
        file.write((char*)&node->y, 1);
        if (node->next_level) {
            saveNode(level + 1, node->next_level, file);
        }
        if (node->next_neighb) {
            saveNode(level, node->next_neighb, file);
        }
    }

    void saveToDot(MCTS_Node *root, const string& fileName) {
        ofstream file;
        file.open(fileName, ios::trunc);
        file << "digraph {\n";
        saveNodeToDot(root, file);
        file << '}';
        file.close();
    }

    void saveNodeToDot(MCTS_Node *node, ofstream& file) {
        file << (uintptr_t)node << "[label=\"[" << (char)(node->x + 96) << "," << (char)(node->y + 96) << "] ";
        file << node->wins << '/' << node->games << "\",color=" << (node->color == 1 ? "red" : "blue") << "]\n";
        if (node->parent)
            file << (uintptr_t)node->parent << "->" << (uintptr_t)node << '\n'; 
        if (node->next_level) {
            saveNodeToDot(node->next_level, file);
        }
        if (node->next_neighb) {
            saveNodeToDot(node->next_neighb, file);
        }
    }

    MCTS_Node* restore(const string& fileName) {
        ifstream file;
        MCTS_Node *node, *root, *prevNode;
        int prevLevel, level;
        double wins;
        uint64_t games;
        uint8_t x, y;
        file.open(fileName, ios::binary);
        file.read((char*)&level, sizeof(level));
        file.read((char*)&wins, sizeof(wins)); 
        file.read((char*)&games, sizeof(games));
        file.read((char*)&x, 1);
        file.read((char*)&y, 1);
        root = new_node(0);
        root->wins = wins;
        root->games = games;
        prevNode = root;
        prevLevel = level;
        while (file) {
            file.read((char*)&level, sizeof(level));
            file.read((char*)&wins, sizeof(wins)); 
            file.read((char*)&games, sizeof(games));
            file.read((char*)&x, 1);
            file.read((char*)&y, 1);
            if (level > prevLevel) {
                node = new_node(prevNode, x, y);
            } else {
                while (level < prevLevel--) {
                    prevNode = prevNode->parent;
                }
                node = new_node(prevNode->parent, x, y);
            }
            node->wins = wins;
            node->games = games;
            prevNode = node;
            prevLevel = level;
        }
        file.close();
        return root;
    }

private:
    MCTS_Node* next_free;
    size_t n;
    int n_chunks;
};

struct Dot {
    explicit Dot(int inColor = 0)
        : color(inColor), next(0), curve(this), prevEmpty(0), nextEmpty(0)
    {}

    operator int () const {
        return color & C_MASK;
    }

    int color, x, y;
    Dot *next, *curve, *prevEmpty, *nextEmpty;

    Dot(Dot const&) = delete;
};

void concatCurves(Dot &a, Dot &b) {
    Dot* d = a.curve;
    Dot* a_next = d->next;
    d->next = b.curve;
    while (d->next) {
        d = d->next;
        d->curve = a.curve;
    }
    d->next = a_next;
}

class Field {
public:
    explicit Field() {
        resetField();
        PR(printField();)
    }

    void resetField() {
        nEmpty = N * N;
        memset(dots, 0, (N + 2) * (N + 2) * sizeof(Dot));
        memset(moves, 0, 2 * N * N * sizeof(int));
        moves_p = moves;
        dots[0][0].nextEmpty = &dots[1][1];
        dots[1][1].prevEmpty = &dots[0][0];
        for (auto i{1}; i <= N; ++i) {
            for (auto j{1}; j <= N - 1; ++j) {
                dots[i][j].nextEmpty = &dots[i][j + 1];
                dots[i][j + 1].prevEmpty = &dots[i][j];
                dots[i][j].x = i;
                dots[i][j].y = j;
                dots[i][j].curve = &dots[i][j];
            }
            dots[i][N].nextEmpty = &dots[i + 1][1];
            dots[i + 1][1].prevEmpty = &dots[i][N];
            dots[i][N].x = i;
            dots[i][N].y = N;
            dots[i][N].curve = &dots[i][N];
        }
        dots[N][N].nextEmpty = 0;
    }

    Dot* getFirstEmpty() {
        return &dots[0][0];
    }

    void bfs(int x, int y, int color, bool *inside, int num_bfs) {
        bool a{false}, b{false}, c{false}, d{false};
        if (bfs1(x - 1, y, color, inside, num_bfs)) a = true, b = true;
        if (bfs1(x + 1, y, color, inside, num_bfs)) c = true, d = true;
        if (bfs1(x, y - 1, color, inside, num_bfs)) a = true, c = true;
        if (bfs1(x, y + 1, color, inside, num_bfs)) b = true, d = true;
        if (a) bfs1(x - 1, y - 1, color, inside, num_bfs);
        if (b) bfs1(x - 1, y + 1, color, inside, num_bfs);
        if (c) bfs1(x + 1, y - 1, color, inside, num_bfs);
        if (d) bfs1(x + 1, y + 1, color, inside, num_bfs);
    }

    bool bfs1(int x, int y, int color, bool *inside, int num_bfs) {
        if (dots[x][y] != color) {
            if (x == 0 || y == 0 || x == N + 1 || y == N + 1) {
                *inside = false;
            } else if (!seen[x][y]) {
                seen[x][y] = num_bfs;
                bfs(x, y, color, inside, num_bfs);
                return true;
            }
        }
        return false;
    }

    void redrawCycle(int color) {
        memset(seen, 0, (N + 2) * (N + 2) * sizeof(int));
        memset(isInside, 0, N * N);
        int num_bfs{1};
        for (auto i{1}; i <= N; ++i) {
            for (auto j{1}; j <= N; ++j) {
                if (!seen[i][j]) {
                    seen[i][j] = num_bfs;
                    if (dots[i][j] != color) {
                        bool inside = true;
                        bfs(i, j, color, &inside, num_bfs);
                        if (inside) {
                            isInside[num_bfs] = 1;
                        }
                        ++num_bfs;
                    }
                }
            }
        }
        nEmpty = 0;
        dots[0][0].nextEmpty = 0;
        for (auto i{1}; i <= N; ++i) {
            for (auto j{1}; j <= N; ++j) {
                if (isInside[seen[i][j]]) {
                    dots[i][j].color &= ~C_MASK; // zero color bits
                    dots[i][j].color |= color;
                }
                if (!dots[i][j]) {
                    dots[i][j].nextEmpty = dots[0][0].nextEmpty;
                    if (dots[0][0].nextEmpty) {
                        dots[0][0].nextEmpty->prevEmpty = &dots[i][j];
                    }
                    dots[0][0].nextEmpty = &dots[i][j];
                    dots[i][j].prevEmpty = &dots[0][0];
                    ++nEmpty;
                }
            }
        }
        //printField();
    }

    bool checkForCycle(int x, int y) const {

        int c = dots[x][y];

        #ifdef NO_CHECK_FOR_CYCLE

            bool result = true;

        #else

            #define D(i,j,k,l) (dots[x + i][y + j].curve == dots[x + k][y + l].curve && dots[x + i][y + j] == c && dots[x + k][y + l] == c)

            bool result = (!D(-1,-1,-1,0) && !D(-1,-1,0,-1) && (D(-1,-1,-1,1) || D(-1,-1,0,1) || D(-1,-1,1,1) || D(-1,-1,1,0) || D(-1,-1,1,-1)))
                       || (!D(-1,0,0,-1) && !D(-1,0,0,1) && (D(-1,0,1,1) || D(-1,0,1,0) || D(-1,0,1,-1)))
                       || (!D(-1,1,-1,0) && !D(-1,1,0,1) && (D(-1,1,1,1) || D(-1,1,1,0) || D(-1,1,1,-1) || D(-1,1,0,-1)))
                       || (!D(0,-1,1,0) && (D(0,-1,0,1) || D(0,-1,1,1)))
                       || (!D(0,1,1,0)  && D(0,1,1,-1))
                       || (!D(1,-1,1,0) && D(1,-1,1,1));
            
            #undef D

        #endif // ndef NO_CHECK_FOR_CYCLE

        PR(cout << (result ? "There is a cycle (" : "There is no cycle yet (") << "\033[3" << c << "m" << (char)(x + 96) << (char)(y + 96) << C0 << ")\n";)

        return result;
    }

    void manageCurves(int x, int y) {
        Dot& d = dots[x][y];
        for (auto i{x - 1}; i < x + 2; ++i) {
            for (auto j{y - 1}; j < y + 2; ++j) {
                if ((dots[i][j].color >> DOT_BITS) == d && (i != x || j != y) && dots[i][j].curve != d.curve) {
                    concatCurves(dots[i][j], d);
                }
            }
        }
    }

    void addDot(int x, int y, int c) {
        *moves_p++ = x;
        *moves_p++ = y;
        if (x > N || y > N || dots[x][y]) {
            printField();
            cout << "IMPOSSIBLE MOVE (" << "\033[3" << c << "m" << (char)(x + 96) << (char)(y + 96) << C0 << ")\n";
            return;
        }
        if (dots[x][y].nextEmpty) {
            dots[x][y].nextEmpty->prevEmpty = dots[x][y].prevEmpty;
        }
        dots[x][y].prevEmpty->nextEmpty = dots[x][y].nextEmpty;
        --nEmpty;

        dots[x][y].color = c | B_SET | (c << DOT_BITS);
        bool cycle = checkForCycle(x, y);
        manageCurves(x, y);
        if (cycle)
            redrawCycle(dots[x][y]);
        PR(printField();)

    }

    void printField() const {
        cout << ' ';
        for (auto i{1}; i <= N; ++i) {
            cout << (char)(i + 96);
        }
        cout << '\n';
        for (auto i{1}; i <= N; ++i) {
            cout << (char)(i + 96);
            for (auto j{1}; j <= N; ++j) {
                cout << "\033[" << (dots[i][j] != (dots[i][j].color >> DOT_BITS) ? "1;3" : "3") <<
                    #ifdef COLOR_CURVES
                        (dots[i][j] ? (intptr_t)dots[i][j].curve % 7 + 1 : 0)
                    #else
                        (dots[i][j])
                    #endif
                    << 'm' << (dots[i][j].color & DOT_MASK) << C0;
            }
            cout << '\n';
        }
    }

    double gameSim(int color) {
        int blue{0}, red{0};
        while (nEmpty) {
            Dot *nextDot = &dots[0][0];
            int n = rand() % nEmpty + 1;
            while (n--) {
                nextDot = nextDot->nextEmpty;
            }
            addDot(nextDot->x, nextDot->y, color);
            (color == 1) ? color = 2 : color = 1;
        }
        for (auto i{1}; i <= N; ++i) {
            for (auto j{1}; j <= N; ++j) {
                red  += (dots[i][j].color & DOT_MASK) == RED_DOT;
                blue += (dots[i][j].color & DOT_MASK) == BLUE_DOT;
            }
        }

        //cout << "Result: " << red << " " << blue << '\n';
        return (red > blue) + ((double)(red == blue)) / 2;
    }

    void printMoves() const {
        for (int i = 0; i < 2 * N * N; i++) {
            cout << (char)(moves[i] + 96);
        }
        cout << '\n';
    }

private:
    Dot  dots[N + 2][N + 2];
    int  seen[N + 2][N + 2];
    int  moves[2 * N * N], *moves_p;
    char isInside[N * N];
    int  nEmpty;
};

void getGameState(MCTS_Node *node, Field *field) {
    if (!(node->parent)) {
        return;
    }
    getGameState(node->parent, field);
    field->addDot(node->x, node->y, node->color);
}




MCTS_Node* uct(MCTS_Node *cur) {
    MCTS_Node *cond = cur, *ans = cond;
    double C = 1.41421356237;
    while (cond->next_level) {
        cond = cond->next_level;
        ans = cond;
        double UCB1 = cond->wins / cond->games + C * pow((log(cur->games) / cond->games), 0.5);
        double MaxUCB = UCB1;
        while(cond->next_neighb) {
            cond = cond->next_neighb;
            UCB1 = cond->wins / cond->games + C * pow((log(cur->games) / cond->games), 0.5);
            if (UCB1 > MaxUCB) {
                MaxUCB = UCB1;
                ans = cond;
            }
        }   
        cond = ans;
    }
    return cond;
}

void updateTree(MCTS_Node *cur, double res) {
    while (cur) {
        ++cur->games;
        cur->wins += cur->color == 1 ? res : 1 - res;
        cur = cur->parent;
    }
}

void mcts(uint64_t steps, int numGames, bool printMoves, bool printEachMove, bool printField, int sims) {
    uint64_t t;
    double result = 0;
    MCTS_Tree tree;
    MCTS_Node *root;
    if (inFile.size()) {
        root = tree.restore(inFile);
    } else {
        root = tree.new_node(0);
    }
    if (dotFile.size()) {
        tree.saveToDot(root, "in-" + dotFile);
    }
    Field gameField;
    Field expField;
    Field simField;
    for (int g{0}; g < numGames; ++g) {
        gameField.resetField();
        MCTS_Node *curNode = root;
        while (true) {
            for (int s{0}; s < sims; ++s) {
                MCTS_Node *nextNodeToExp = uct(curNode);
                expField.resetField();
                getGameState(nextNodeToExp, &expField);
                Dot *nextDot = expField.getFirstEmpty();
                while (nextDot->nextEmpty) {
                    simField.resetField();
                    getGameState(nextNodeToExp, &simField);
                    nextDot = nextDot->nextEmpty;
                    simField.addDot(nextDot->x, nextDot->y, nextNodeToExp->color == 1 ? 2 : 1);
                    MCTS_Node *newNode = tree.new_node(nextNodeToExp, nextDot->x, nextDot->y);
                    t = tsc();
                    result = simField.gameSim(nextNodeToExp->color);
                    updateTree(newNode, result);
                    //cout << newNode->wins << "[" << tsc() - t << "] ";
                    if (!--steps) {
                        goto end;
                    }
                }
            }   
            MCTS_Node *nextMove = curNode->next_level, *ans = nextMove;
            if (!nextMove) {
                cerr << "Result: " << result << '\n';
                break;
            }
            double score = nextMove->wins / nextMove->games;
            double maxScore = score;
            while(nextMove->next_neighb) {
                nextMove = nextMove->next_neighb;
                score = nextMove->wins / nextMove->games;
                if (score > maxScore) {
                    maxScore = score;
                    ans = nextMove;
                }
            }
            gameField.addDot(ans->x, ans->y, ans->color);
            if (printEachMove)
                cout << (char)(ans->x + 96) << (char)(ans->y + 96) << "\n";
            if (printField)
                gameField.printField();
            curNode = ans;
        }
    }

end:
    if (outFile.size()) {
        tree.save(root, outFile);
    }
    if (dotFile.size()) {
        tree.saveToDot(root, dotFile);
    }
    if (printMoves) {
        simField.printMoves();
    }
}


void playGame(int c, int numGames, bool printField, bool printEachMove, int sims) {
    int color{c};
    char inp_x, inp_y;
    int x, y;
    uint64_t t;
    double result = 0;
    MCTS_Tree tree;
    MCTS_Node *root;
    if (inFile.size()) {
        root = tree.restore(inFile);
    } else {
        root = tree.new_node(0);
    }
    if (dotFile.size()) {
        tree.saveToDot(root, "in-" + dotFile);
    }


    Field gameField;
    Field expField;
    Field simField;
    //gameField.printField();
    //uint64_t t = tsc();
    for (int g{0}; g < numGames; ++g) {
        gameField.resetField();
        MCTS_Node *curNode = root;
        while (true) {
            for (int s{0}; s < sims; ++s) {
                MCTS_Node *nextNodeToExp = uct(curNode);
                expField.resetField();
                getGameState(nextNodeToExp, &expField);
                Dot *nextDot = expField.getFirstEmpty();
                while (nextDot->nextEmpty) {
                    simField.resetField();
                    getGameState(nextNodeToExp, &simField);
                    nextDot = nextDot->nextEmpty;
                    simField.addDot(nextDot->x, nextDot->y, nextNodeToExp->color == 1 ? 2 : 1);
                    MCTS_Node *newNode = tree.new_node(nextNodeToExp, nextDot->x, nextDot->y);
                    t = tsc();
                    result = simField.gameSim(nextNodeToExp->color);
                    updateTree(newNode, result);
                    //cout << newNode->wins << "[" << tsc() - t << "] ";
                    //if (!--steps) {
                      //  goto end;
                    //}
                }
            }   
            MCTS_Node *nextMove = curNode->next_level, *ans = nextMove;

            if (!nextMove) {
                cerr << "Result: " << result << '\n';
                break;
            }

            if (ans->color == color) {
                cin >> inp_x >> inp_y;
                x = inp_x - 96;
                y = inp_y - 96;
                if (!cin) {
                    break;
                }
                while(nextMove->next_neighb) {
                    nextMove = nextMove->next_neighb;
                    if (nextMove->x == x && nextMove->y == y) {
                        ans = nextMove;
                    }
                }
            } else {
                double score = nextMove->wins / nextMove->games;
                double maxScore = score;
                while(nextMove->next_neighb) {
                    nextMove = nextMove->next_neighb;
                    score = nextMove->wins / nextMove->games;
                    if (score > maxScore) {
                        maxScore = score;
                        ans = nextMove;
                    }
                }
            }
            gameField.addDot(ans->x, ans->y, ans->color);
            if (printEachMove && (ans->color != color))
                    cout << (char)(ans->x + 96) << (char)(ans->y + 96) << "\n";
            if (printField)
                gameField.printField();
            curNode = ans;
        }
    }


end:
    if (outFile.size()) {
        tree.save(root, outFile);
    }
    if (dotFile.size()) {
        tree.saveToDot(root, dotFile);
    }
    //if (printMoves) {
      //  simField.printMoves();
    //}
}

int main(int argc, char *argv[]) {
    bool play{false}, optPlay{false}, optPrintFinalMoves{false}, optPrintEachMove{false}, optPrintField{false}, optSteps{false}, optSims{false}, optSeed{false}, optNumGames{false}, optInpFile{false}, optOutFile{false}, optDotFile{false};
    uint64_t steps{-1UL};
    unsigned seed{1}, numGames{1}, sims{1};
    int color;
    for (auto i{1}; i < argc; i++) {
        for (auto j{0}; argv[i][j] && j < 10; j++) {
            if (optPlay) {
                color = atoi(argv[i] + j);
                if (!color) {
                    goto wrong_int;
                }
                optPlay = false;
                play = true;
                break;
            }
            if (optSteps) {
                steps = atol(argv[i] + j);
                if (!steps) {
                    goto wrong_int;
                }
                optSteps = false;
                break;
            }
            if (optSims) {
                sims = atoi(argv[i] + j);
                if (!sims) {
                    goto wrong_int;
                }
                optSteps = false;
                break;
            }
            if (optSeed) {
                seed = atoi(argv[i] + j);
                if (!seed) {
                    goto wrong_int;
                }
                optSeed = false;
                break;
            }
            if (optNumGames) {
                numGames = atol(argv[i] + j);
                if (!numGames) {
                    goto wrong_int;
                }
                optNumGames = false;
                break;
            }
            if (optInpFile) {
                inFile = argv[i] + j;
                optInpFile = false;
                break;
            }
            if (optOutFile) {
                outFile = argv[i] + j;
                optOutFile = false;
                break;
            }
            if (optDotFile) {
                dotFile = argv[i] + j;
                optDotFile = false;
                break;
            }

            switch (argv[i][j]) {
                case 'p':
                    optPlay = true;
                    break;
                case 'e':
                    optPrintFinalMoves = true;
                    break;
                case 'm':
                    optPrintEachMove = true;
                    break;
                case 'f':
                    optPrintField = true;
                    break;
                case 'r':
                    optSeed = true;
                    break;
                case 's':
                    optSteps = true;
                    break;
                case 'S':
                    optSims = true;
                    break;
                case 'g':
                    optNumGames = true;
                    break;
                case 'i':
                    optInpFile = true;
                    break;
                case 'o':
                    optOutFile = true;
                    break;
                case 'd':
                    optDotFile = true;
                    break;
                case '-':
                    if (j == 0)
                        break;
                    [[fallthrough]];
                default:
                    goto print_usage;
            }
        }
    }
    if (optSteps || optSeed) {
        goto wrong_int;
    }

    if (play) {
        playGame(color, numGames, optPrintField, optPrintEachMove, sims);
        return 0;
    }

    srand(seed);
    mcts(steps, numGames, optPrintFinalMoves, optPrintEachMove, optPrintField, sims);
    cout << '\n';

    return 0;

wrong_int:
    cout << "Expected integer n > 0 option value\n\n";

print_usage:
    cout << "Usage: " << argv[0] << " [-efm] [-p <n>] [-r <n>] [-s <n>]\n\n";
    cout << "Options:\n";
    cout << "  e:        print all moves at the end of the game\n";
    cout << "  f:        print field after each move\n";
    cout << "  m:        print each move\n";
    cout << "  p <n>:    interactive play with color n\n";
    cout << "  r <n>:    use n as a random seed\n";
    cout << "  s <n>:    do n simulations, then stop\n";
    cout << "  g <n>:    number of games\n";
    cout << "  S <n>:    do n mcts simulations, then stop\n";
    cout << "  i <file>: input mcts\n";
    cout << "  o <file>: output mcts\n";
    cout << "  d <file>: output dot\n";
    return -1;
}

