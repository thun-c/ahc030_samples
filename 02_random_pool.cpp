/*
ランダムに配置候補を生成
*/
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <unordered_map>
#include <deque>
#include <algorithm>
#include <bitset>
#include <chrono>
#include <tuple>
#include <iterator>
#include <iomanip>
#include <cassert>
#include <cstdint>
#include <functional>

// #define cerr   \
//     if (false) \
//     cerr

using namespace std;

const vector<pair<int, int>> DIJ = {{0, 1}, {1, 0}, {0, -1}, {-1, 0}};

class Xorshift
{
public:
    using result_type = uint64_t;

    explicit Xorshift(uint64_t seed) : x_(seed) { assert(seed); }

    // [0, stop)
    uint64_t randrange(uint64_t stop)
    {
        assert(stop > 0);
        next();
        return x_ % stop;
    }

    // [start, stop)
    uint64_t randrange(uint64_t start, uint64_t stop)
    {
        assert(start < stop);
        next();
        return start + x_ % (stop - start);
    }

    // [a, b]
    uint64_t randint(uint64_t a, uint64_t b)
    {
        assert(a <= b);
        return randrange(a, b + 1);
    }

    // [0.0, 1.0]
    double random()
    {
        next();
        return static_cast<double>(x_) * (1.0 / static_cast<double>(UINT64_MAX));
    }

    // [a, b] or [b, a]
    double uniform(double a, double b) { return a + (b - a) * random(); }

    uint64_t next()
    {
        x_ ^= x_ << 13;
        x_ ^= x_ >> 17;
        x_ ^= x_ << 5;
        return x_;
    }

    bool gen_bool(double p)
    {
        return random() < p;
    }

    // Add this operator to use with std::shuffle
    uint64_t operator()()
    {
        return next();
    }

    static constexpr uint64_t min() { return 0; }
    static constexpr uint64_t max() { return UINT64_MAX; }

private:
    uint64_t x_;
};
Xorshift rng(1);

// 油田の配置についての情報
struct OilLayout
{
    uint64_t hash;
    // 油田配置をこの配置xだと過程したとき
    // 今までのクエリ結果Rが得られる確率の対数.
    // 対数尤度ともいう。
    double ln_pR_if_x;
    // 今までのクエリの結果Rから計算した,この配置になる事後確率P(x|R)
    double px_if_R;
    vector<size_t> top_lefts; // 油田の左上の座標
    vector<uint8_t> volume;   // ある位置の油の埋蔵量
};

// 油田の形についての情報
struct OilShape
{
    size_t max_i, max_j;                      // 油田が収まる正方形の大きさ
    vector<size_t> coordinate_ids;            // 座標(i,j)の組を1変数で表したもの
    vector<pair<size_t, size_t>> coordinates; // 座標(i,j)の組
    // 座標(i,j)の組をマスクしたもの
    // 島の大きさN<=20なので、20*20のビットセットで表現できる
    bitset<20 * 20> mask;
};

struct Input
{
    size_t n;              // 島の大きさ 10<=N<=20
    size_t n2;             // n*n
    size_t m;              // 油田の数 2<=M<=10
    double eps;            // 占い結果に用いるエラーパラメータ 0.01<=eps<=0.2
    vector<OilShape> oils; // 油田の形に関する情報 .size()==M
    size_t total;          // 島全体の油田の埋蔵量の合計

    // 油田の左上の座標を受け取り、埋蔵量が1以上のマスの集合をbitsetで返す
    bitset<20 * 20> get_positives(const vector<size_t> &top_lefts) const
    {
        bitset<20 * 20> positives;
        for (size_t oil_id = 0; oil_id < m; ++oil_id)
        {
            positives |= oils[oil_id].mask << top_lefts[oil_id];
        }
        return positives;
    }

    // M個の油田の左上の座標を受け取り、その位置の油の埋蔵量を返す
    vector<uint8_t> get_volume(const vector<size_t> &top_lefts) const
    {
        vector<uint8_t> volume(n2, 0);
        for (size_t oil_id = 0; oil_id < top_lefts.size(); ++oil_id)
        {
            size_t pij = top_lefts[oil_id];
            for (size_t ij : oils[oil_id].coordinate_ids)
            {
                volume[ij + pij] += 1;
            }
        }
        return volume;
    }
};

// プログラムが始まってからの時間を計測する
double get_time()
{
    // startをstaticにすることで、
    // 2回目以降のget_time()の呼び出しでも
    // プログラムが始まってからの時間を計測できる
    static auto start = chrono::steady_clock::now();
    auto now = chrono::steady_clock::now();
    return chrono::duration<double>(now - start).count();
}

Input read_input()
{
    Input input;
    cin >> input.n >> input.m >> input.eps;
    input.n2 = input.n * input.n;
    input.oils.resize(input.m);

    for (size_t oil_id = 0; oil_id < input.m; ++oil_id)
    {
        size_t t_size;
        cin >> t_size;
        input.oils[oil_id].coordinates.resize(t_size);
        for (size_t i = 0; i < t_size; ++i)
        {
            size_t x, y;
            cin >> x >> y;
            input.oils[oil_id].coordinates[i] = make_pair(x, y);
        }
    }
    input.total = 0;
    for (size_t oil_id = 0; oil_id < input.m; ++oil_id)
    {
        input.total += input.oils[oil_id].coordinates.size();
    }

    // 同じ図形かの判定を行うためにポリオミノをソートしてく
    sort(input.oils.begin(), input.oils.end(), [](const OilShape &a, const OilShape &b)
         { return a.coordinates < b.coordinates; });

    for (size_t oil_id = 0; oil_id < input.m; ++oil_id)
    {
        auto &oil = input.oils[oil_id];
        oil.max_i = 0;
        oil.max_j = 0;
        for (const auto &[i, j] : oil.coordinates)
        {
            oil.max_i = max(oil.max_i, i);
            oil.max_j = max(oil.max_j, j);
        }
        oil.coordinate_ids.resize(oil.coordinates.size());
        for (size_t i = 0; i < oil.coordinates.size(); ++i)
        {
            oil.coordinate_ids[i] = oil.coordinates[i].first * input.n + oil.coordinates[i].second;
        }
        oil.mask.reset();
        for (size_t ij : oil.coordinate_ids)
        {
            oil.mask.set(ij, true);
        }
    }

    return input;
}

const double SMALL_VALUE = 1e-6; // すごく小さい値

// 油田1個分の状態
struct OilState
{
public:
    // top_left_query_volumes[top_left][q]は、
    // この油田の左上座標がtop_leftにあるとき、
    // q番目のクエリで占った座標集合の埋蔵量の合計を示している
    vector<vector<uint8_t>> top_left_query_volumes;
    // この油田の左上座標ごとのzobrist hash
    vector<uint64_t> hashes;
    OilState(const Input &input) : top_left_query_volumes(input.n2),
                                   hashes(input.n2)
    {
    }
};

class State
{
public:
    vector<OilState> oil_states;   // 油田の状態についてのリスト: M個
    vector<size_t> top_lefts;      // 油田の左上の座標. oil_satesに含めたかったが、OilLayoutにコピーして使うので別で持つ: M個
    vector<uint8_t> volumes;       // ある位置の油の埋蔵量: N*N個
    vector<uint8_t> query_volumes; // q番目のクエリで占った座標集合の埋蔵量の合計
    uint64_t hash;
    const Input &input;

    // 全ての油田が0,0にあるときの状態を初期状態とする
    State(const Input &input) : input(input), oil_states(input.m, OilState(input))
    {
        for (size_t oil_id = 0; oil_id < input.m; ++oil_id)
        {
            auto &oil_state = oil_states[oil_id];
            if (oil_id > 0 && input.oils[oil_id - 1].coordinate_ids == input.oils[oil_id].coordinate_ids)
            {
                oil_state.hashes = oil_states[oil_id - 1].hashes;
            }
            else
            {
                for (size_t ij = 0; ij < input.n2; ++ij)
                {
                    oil_state.hashes[ij] = rng.next();
                }
            }
        }
        // zobrist hashを初期化
        hash = 0;
        for (auto &oil_state : oil_states)
        {
            hash ^= oil_state.hashes[0];
        }
        top_lefts.resize(input.m, 0);
    }

    // 油田oil_idをnew_top_leftに移動する
    // 移動自体はtop_leftsの更新だが、
    // 移動に伴い、hash,query_volumes,volumesを更新する
    void move_to(size_t oil_id, size_t new_top_left)
    {
        auto &oil_state = oil_states[oil_id];
        hash ^= oil_state.hashes[top_lefts[oil_id]] ^ oil_state.hashes[new_top_left];
        for (size_t q = 0; q < query_volumes.size(); ++q)
        {
            query_volumes[q] += oil_state.top_left_query_volumes[new_top_left][q] - oil_state.top_left_query_volumes[top_lefts[oil_id]][q];
        }
        if (!volumes.empty())
        {
            for (size_t ij : input.oils[oil_id].coordinate_ids)
            {
                volumes[ij + top_lefts[oil_id]] -= 1;
                volumes[ij + new_top_left] += 1;
            }
        }
        top_lefts[oil_id] = new_top_left;
    }

    // 占いクエリを投げた後、今回占った座標集合に含まれる埋蔵量を記録する
    void add_query(const vector<size_t> &query_coordinates)
    {
        vector<bool> in_query(input.n2, false);
        for (size_t ij : query_coordinates)
        {
            in_query[ij] = true;
        }
        for (size_t oil_id = 0; oil_id < input.m; ++oil_id)
        {
            auto &oil = input.oils[oil_id];
            auto &oil_state = oil_states[oil_id];
            for (size_t di = 0; di < input.n - oil.max_i; ++di)
            {
                for (size_t dj = 0; dj < input.n - oil.max_j; ++dj)
                {
                    size_t top_left = di * input.n + dj;
                    uint8_t c = 0;
                    for (size_t ij : input.oils[oil_id].coordinate_ids)
                    {
                        if (in_query[top_left + ij])
                        {
                            c += 1;
                        }
                    }
                    oil_state.top_left_query_volumes[top_left].push_back(c);
                }
            }
        }
        vector<uint8_t> volume = input.get_volume(top_lefts);
        uint8_t c = 0;
        for (size_t ij : query_coordinates)
        {
            c += volume[ij];
        }
        query_volumes.push_back(c);
    }
};

class Sim
{
public:
    size_t n, n2, m, total;
    double eps;
    // 過去の占いの(油田配置、占い結果)の集合
    vector<pair<vector<size_t>, size_t>> queries;
    // 既に油田配置を答えるクエリを投げて失敗した油田配置の集合
    vector<vector<size_t>> failed;
    // クエリサイズk、埋蔵量総量Sに対して、
    // pr_if_xがもつrの値の下限
    vector<vector<size_t>> pr_if_x_lb;
    // 真の配置xを過程したときに占い結果がrになる確率(尤度とも呼ぶ)
    // 真の配置xを過程したときはクエリサイズk、埋蔵量総量Sも固定されるため、
    // クエリサイズk、埋蔵量総量Sの時に占い結果がrになる確率を記録しておけばいい
    // 小さすぎる確率は無視するため、配列はr=lb以上のものだけ格納する。
    // pr_if_x[k][S][r-lb] = (prob, log(prob))
    vector<vector<vector<pair<double, double>>>> pr_if_x;
    // クエリごとにあり得る埋蔵量総量Sごとに
    // 埋蔵量Sのときにそのクエリで得られた結果になる確率P(r|S)の対数を記録
    vector<vector<double>> ln_pr_if_s_query;
    // 残りクエリ回数. 2*N*N回までクエリを投げられる
    size_t rem;

    Sim(const Input &input) : n(input.n), n2(input.n2), m(input.m), total(input.total), eps(input.eps), rem(input.n * input.n * 2)
    {
        // クエリサイズk、埋蔵量総量S、クエリの結果rに対する尤度は事前に計算しておく
        pr_if_x_lb.resize(n * n + 1, vector<size_t>(total + 1));
        pr_if_x.resize(n * n + 1, vector<vector<pair<double, double>>>(total + 1));
        for (size_t k = 1; k <= n * n; ++k)
        {
            for (size_t S = 0; S <= total; ++S)
            {
                // muとsigmaは問題文の計算式をそのまま使うだけで求まる
                double mu = ((double)k - (double)S) * eps + (double)S * (1.0 - eps);
                double sigma = sqrt(k * eps * (1.0 - eps));
                // 尤度が小さすぎるrがどこかわからないため、
                // 尤度が最も高いmuから順に尤度を求め、
                // 尤度が小さすぎたタイミングで止める
                for (int r = static_cast<int>(round(mu)); r >= 0; --r)
                {
                    double prob = likelihood(mu, sigma, r);
                    if (prob < SMALL_VALUE)
                    {
                        pr_if_x_lb[k][S] = r + 1;
                        break;
                    }
                    pr_if_x[k][S].emplace_back(prob, log(prob));
                }
                // rの値について降順になっているため、昇順になおす
                reverse(pr_if_x[k][S].begin(), pr_if_x[k][S].end());
                // muを基準に対称なので、muより大きいrについても同様に計算する
                for (int r = static_cast<int>(round(mu)) + 1;; ++r)
                {
                    double prob = likelihood(mu, sigma, r);
                    if (prob < SMALL_VALUE)
                    {
                        break;
                    }
                    pr_if_x[k][S].emplace_back(prob, log(prob));
                }
            }
        }
    }

    // 油田配置をあてるクエリを投げる
    // 失敗した場合は失敗した座標郡をfailedに記録する
    bool ans(const vector<size_t> &T)
    {
        if (rem == 0)
        {
            cerr << "!log giveup " << endl;
            exit(0);
        }
        --rem;
        cout << "a " << T.size() << " ";
        for (size_t ij : T)
        {
            cout << ij / n << " " << ij % n << " ";
        }
        cout << endl;
        size_t ret;
        cin >> ret;
        if (ret == 1)
        {
            return true;
        }
        failed.push_back(T);
        return false;
    }

    // 指定したマスの集合の埋蔵量を占う
    size_t query(const vector<size_t> &query_coordinates)
    {
        if (rem == 0)
        {
            cerr << "!log giveup " << endl;
            exit(0);
        }
        --rem;
        cout << "q " << query_coordinates.size() << " ";
        for (size_t ij : query_coordinates)
        {
            cout << ij / n << " " << ij % n << " ";
        }
        cout << endl;
        size_t ret;
        cin >> ret;
        // クエリの結果を記録する
        queries.emplace_back(query_coordinates, ret);
        // 結果retが得られた.
        // 指定したマス集合の真の埋蔵量総量がわからないため、
        // あり得る埋蔵量総量全てについて、
        // 埋蔵量Sのときに結果がretになる確率P(ret|S)を求める
        vector<double> ln_pr_if_s(total + 1);
        size_t k = query_coordinates.size();
        for (size_t S = 0; S <= total; ++S)
        {
            double kS = (double)k - (double)S;
            double kSeps = kS * eps;
            double meps = 1.0 - eps;
            double mu = kSeps + S * meps;
            double sigma = sqrt(k * eps * meps);
            ln_pr_if_s[S] = log(likelihood(mu, sigma, ret));
        }

        ln_pr_if_s_query.push_back(ln_pr_if_s);
        return ret;
    }

    // 指定したマスの埋蔵量を取得する
    // 1マスなら正確な値がわかる
    size_t mine(size_t i, size_t j)
    {
        if (rem == 0)
        {
            cerr << "!log giveup" << endl;
            exit(0);
        }
        --rem;
        cout << "q 1 " << i << " " << j << endl;
        size_t ret;
        cin >> ret;
        return ret;
    }

    // 平均mean, 標準偏差std_devに従う正規分布において、
    // resが観測される確率を求める
    double likelihood(double mean, double std_dev, size_t res)
    {
        // 占い結果はμとσ^2に従う正規分布からサンプルされたxそのものではなく、
        // max(0,round(x))である。
        // res=rounnd(x)について考えると、
        // 占い結果がresになる確率はres-0.5<=x<res+0.5の範囲の確率分布の面積と同等である。
        // なぜなら、round(res-0.5)=res, round(res+0.5)=res+1となり、
        // res-0.5<=x<res+0.5の範囲全てで同じ値をとるからである。
        // また、max(0,round(x))について,
        // x<0の場合は常に0であるから、
        // x=-1でもx=0.5でもx=0でもres=0である。
        // よって、res=0の場合は
        // -∞<=x<0.5の確率分布の面積を求める。
        // 区間[a,b)の確率分布の面積は
        // 累積分布関数の差cdf(b)-cdf(a)で求められる。
        // 区間[-∞,b)の確率分布の面積は
        // 累積分布関数cdf(b)で求められる。
        double b = (double)res + 0.5;

        if (res == 0)
        {
            return normal_cdf(mean, std_dev, b);
        }
        else
        {
            double a = (double)res - 0.5;
            return normal_cdf(mean, std_dev, b) - normal_cdf(mean, std_dev, a);
        }
    }

    // 平均mean, 標準偏差std_devに従う正規分布において、a以上b以下の確率を求める
    double probability_in_range(double mean, double std_dev, double a, double b)
    {
        double p_a = normal_cdf(mean, std_dev, a);
        double p_b = normal_cdf(mean, std_dev, b);
        return p_b - p_a;
    }

    // 累積分布関数
    // 平均mean, 標準偏差std_devに従う正規分布において、x以下の確率を求める
    double normal_cdf(double mean, double std_dev, double x)
    {
        return 0.5 * (1.0 + erf((x - mean) / (std_dev * sqrt(2.0))));
    }

    /// 時間切れしたときはBFSで掘る。
    void giveup()
    {
        cerr << "!log giveup" << endl;
        deque<pair<size_t, size_t>> que;
        que.emplace_back(n / 2, n / 2);
        vector<size_t> list;
        size_t rem = total;
        vector<vector<bool>> used(n, vector<bool>(n, false));

        while (!que.empty())
        {
            auto [i, j] = que.front();
            que.pop_front();
            if (used[i][j])
                continue;
            used[i][j] = true;

            size_t ret = mine(i, j);
            if (ret > 0)
            {
                list.push_back(i * n + j);
                rem -= ret;
                if (rem == 0)
                    break;
            }

            for (const auto &[di, dj] : DIJ)
            {
                size_t i2 = i + di;
                size_t j2 = j + dj;
                if (i2 >= 0 && i2 < n && j2 >= 0 && j2 < n)
                {
                    if (ret == 0)
                    {
                        que.emplace_back(i2, j2);
                    }
                    else
                    {
                        que.emplace_front(i2, j2);
                    }
                }
            }
        }
        ans(list);
    }

    // volumesとfailed_coordinatesが異なるかどうかを返す
    bool is_different(const vector<uint8_t> &volumes, const vector<size_t> &failed_coordinates) const
    {
        for (const auto &ij : failed_coordinates)
        {
            if (volumes[ij] == 0)
            {
                return true;
            }
        }
        return false;
    }

    // 油田配置がtop_leftsにある場合、
    // q番目のクエリで占った油田集合の埋蔵量合計
    uint8_t get_query_volume(const vector<OilState> &oil_states, size_t q, const vector<size_t> &top_lefts) const
    {
        uint8_t S = 0;
        for (size_t oil_id = 0; oil_id < top_lefts.size(); ++oil_id)
        {
            auto &oil_state_p = oil_states[oil_id];
            auto ij = top_lefts[oil_id];
            size_t p_volume = oil_state_p.top_left_query_volumes[ij][q];
            if (p_volume > 0)
            {
                S += p_volume;
            }
        }
        return S;
    }

    // 油田配置がこの状態になる確率を求める
    // vs: 各座標の埋蔵量
    // top_lefts: 油田の左上座標郡
    double get_ln_pR_if_x(const vector<OilState> &oil_states, const vector<uint8_t> &volumes, const vector<size_t> &top_lefts) const
    {
        // 既に失敗した油田配置の集合に含まれる油田配置の場合、対数尤度を非常に小さい値にする
        for (const auto &failed_coordinates : failed)
        {
            bool skip = is_different(volumes, failed_coordinates);

            if (!skip)
            {
                size_t size = 0;
                for (size_t ij = 0; ij < n2; ++ij)
                {
                    if (volumes[ij] > 0)
                    {
                        ++size;
                    }
                }
                if (size == failed_coordinates.size())
                {
                    return -1e20;
                }
            }
        }

        // 今までのクエリ結果Rから、
        // 配置x=oil_statesにおける尤度P(R|x)を求める
        // 公式は以下の通り
        // P(R|x) = ΠP(ret|x) for ret in R
        // P(ret|x)は非常に小さい値でdoubleに収まらない可能性があるため、対数尤度を計算する
        // 対数に変形すると
        // log(P(R|x)) = Σlog(P(ret|x)) for ret in R
        // となる
        double ln_pR_if_x = 0.0;
        for (size_t q = 0; q < queries.size(); ++q)
        {
            // この油田配置において、q番目のクエリで占った油田の埋蔵量の合計を求める
            // q番目のクエリを打った時のlog(P(ret|S))は記録済みであり、
            // 配置xにおけるSを求めることで、
            // log(P(ret|x)) = log(P(ret|S))を求めることができる
            uint8_t S = get_query_volume(oil_states, q, top_lefts);
            ln_pR_if_x += ln_pr_if_s_query[q][S];
        }
        return ln_pR_if_x;
    }

    double ln_prob_state(const State &state) const
    {
        for (const auto &failed_coordinates : failed)
        {
            bool skip = is_different(state.volumes, failed_coordinates);
            if (!skip)
            {
                size_t size = 0;
                for (size_t ij = 0; ij < n2; ++ij)
                {
                    if (state.volumes[ij] > 0)
                    {
                        ++size;
                    }
                }
                if (size == failed_coordinates.size())
                {
                    return -1e20;
                }
            }
        }
        double prob = 0.0;
        for (size_t q = 0; q < ln_pr_if_s_query.size(); ++q)
        {
            prob += ln_pr_if_s_query[q][state.query_volumes[q]];
        }
        return prob;
    }
};

struct Query
{
    vector<bool> in_query;  // ある位置の油の埋蔵量がクエリされているか : N*N個
    vector<uint8_t> volume; // 油田の埋蔵量のリスト : M個
    size_t coordinate_size; // クエリに含めるマスの数
    vector<OilLayout> pool; // 油田の状態についてのリスト

    Query(const Input &input, vector<OilLayout> &pool) : in_query(input.n2, false), volume(pool.size(), 0), coordinate_size(0), pool(pool)
    {
    }
    // 指定したマスをクエリに含めるかどうかを反転する
    void flip(size_t ij)
    {
        if (in_query[ij])
        {
            in_query[ij] = false;
            for (size_t x = 0; x < pool.size(); ++x)
            {
                volume[x] -= pool[x].volume[ij];
            }
            --coordinate_size;
        }
        else
        {
            in_query[ij] = true;
            for (size_t x = 0; x < pool.size(); ++x)
            {
                volume[x] += pool[x].volume[ij];
            }
            ++coordinate_size;
        }
    }
    // 占うマス集合を取得する
    vector<size_t> get_coordinates() const
    {
        vector<size_t> result;
        for (size_t ij = 0; ij < in_query.size(); ++ij)
        {
            if (in_query[ij])
            {
                result.push_back(ij);
            }
        }
        return result;
    }
    // 占いの良さを評価する
    // add_k: クエリに含める座標の数. 0ならすでにインスタンスに含めているものだけ使う
    // add_v: クエリに含める座標の埋蔵量合計. 0ならすでにインスタンスに含めているものだけ使う
    double eval(Sim &sim, size_t add_k, uint8_t add_v) const
    {
        size_t k = coordinate_size + add_k;
        // ln_pr[r] = クエリ結果がrとなる確率の対数
        // 最後にlogをとるまで、普通の確率として計算する
        auto ln_pr = vector<double>();
        for (size_t x = 0; x < pool.size(); ++x)
        {
            size_t v = volume[x] + add_v;
            size_t lb = sim.pr_if_x_lb[k][v];
            if (ln_pr.size() < lb + sim.pr_if_x[k][v].size())
            {
                ln_pr.resize(lb + sim.pr_if_x[k][v].size());
            }
            double px = pool[x].px_if_R;
            // 公式 p(r)=Σp(r|x)p(x)
            // この公式を使って、p(r)を求める
            for (int pi = 0; pi < sim.pr_if_x[k][v].size(); ++pi)
            {
                const auto &[pr_if_x, ln_pr_if_x] = sim.pr_if_x[k][v][pi];
                ln_pr[lb + pi] += pr_if_x * px;
            }
        }
        for (size_t x = 0; x < ln_pr.size(); ++x)
        {
            ln_pr[x] = log(ln_pr[x]);
        }
        // I(X;R) = ΣΣp(x,r)log(p(x,r)/(p(x)p(r))
        //        = ΣΣp(r|x)p(x)log(p(r|x)p(x)/p(x)p(r))
        //        = ΣΣp(r|x)p(x)log(p(r|x)/p(r))
        //        = ΣΣp(r|x)p(x)(log(p(r|x))-log(p(r)))
        // この式を使って、相互情報量を求める
        double info = 0.0;
        for (size_t x = 0; x < pool.size(); ++x)
        {
            double px = pool[x].px_if_R;
            size_t v = volume[x] + add_v;
            size_t lb = sim.pr_if_x_lb[k][v];
            for (int pi = 0; pi < sim.pr_if_x[k][v].size(); ++pi)
            {
                const auto &[pr_if_x, ln_pr_if_x] = sim.pr_if_x[k][v][pi];
                const auto &ln_prr = ln_pr[lb + pi];
                info += pr_if_x * px * (ln_pr_if_x - ln_prr);
            }
        }
        // コストは1/sqrt(k)なので、
        // sqrt(k)倍することは、コストで割ることと同じ
        // コストあたりの相互情報慮を求める
        return info * sqrt((double)k);
    }
};

// プールの確率を正規化する
void normalize_pool(vector<OilLayout> &pool)
{
    double total = 0;
    for (const auto &layout : pool)
    {
        total += layout.px_if_R;
    }
    for (auto &layout : pool)
    {
        layout.px_if_R /= total;
    }
}

// プールの油田配置の全座標の埋蔵量を計算する
void set_volume(vector<OilLayout> &pool, const Input &input)
{
    for (auto &layout : pool)
    {
        layout.volume = input.get_volume(layout.top_lefts);
    }
}

// プールを小さくして探索範囲を狭める
void concat_pool(vector<OilLayout> &pool, const Input &input, size_t ITER)
{
    double tmp1 = ITER * 0.01;
    double tmp2 = min(3.0 - get_time(), 1.0);
    size_t size = min(pool.size(), max((size_t)(tmp1 * tmp2), (size_t)2));
    // 尤度が低すぎる要素を削除
    while (size > 2 && pool[0].px_if_R * 1e-4 > pool[size - 1].px_if_R)
    {
        --size;
    }
    if (size > 0)
    {
        pool.resize(size);
    }
}

// 占いクエリを取得する
vector<size_t> getDivinationQuery(
    const Input &input, vector<OilLayout> &pool,
    Sim &sim)
{
    const auto size = pool.size();
    // 全ての油田配置が同じ埋蔵量を持つかどうかを座標ごとに調べる
    vector<bool> same(input.n2, true);
    for (size_t x = 1; x < size; ++x)
    {
        const auto &layout = pool[x];
        for (size_t ij = 0; ij < input.n2; ++ij)
        {
            same[ij] = same[ij] && layout.volume[ij] == pool[0].volume[ij];
        }
    }

    // クエリを作成
    Query query(input, pool);
    // 全ての配置が同じ埋蔵量だと想定される座標以外で1点だけ占うクエリを評価する
    vector<pair<double, size_t>> query_coordinate_evals;
    for (size_t ij = 0; ij < input.n2; ++ij)
    {
        if (!same[ij])
        {
            // 評価前後にflipを挟むことで、
            // この座標1個だけでクエリを作成するときの情報量を計算する
            query.flip(ij);                    // クエリにこの座標を含める
            double ev = query.eval(sim, 0, 0); // クエリにこの座標を含めたときの情報量を計算
            query.flip(ij);                    // クエリにこの座標を含めない
            query_coordinate_evals.emplace_back(ev, ij);
        }
    }
    sort(query_coordinate_evals.begin(), query_coordinate_evals.end(), greater<pair<double, size_t>>());

    // 全ての配置が同じ埋蔵量だと想定される座標をクエリに含める
    vector<size_t> no_info_coordinates;
    for (size_t ij = 0; ij < input.n2; ++ij)
    {
        if (same[ij])
        {
            no_info_coordinates.push_back(ij);
        }
    }

    // クエリに含める座標をランダムにソート
    std::vector<std::pair<size_t, size_t>> evaluation_values;
    for (size_t ij : no_info_coordinates)
    {
        size_t evaluation_value = pool[0].volume[ij] * 1000 + rng.randrange(1000);
        evaluation_values.emplace_back(ij, evaluation_value);
    }

    std::sort(evaluation_values.begin(), evaluation_values.end(),
              [](const std::pair<size_t, size_t> &a, const std::pair<size_t, size_t> &b)
              {
                  return a.second > b.second;
              });
    no_info_coordinates = std::vector<size_t>(evaluation_values.size());
    std::transform(evaluation_values.begin(), evaluation_values.end(), no_info_coordinates.begin(),
                   [](const std::pair<size_t, size_t> &pair)
                   {
                       return pair.first;
                   });

    double crt = -1e100;
    size_t add_k = 0; // 情報量0のマスを追加する数
    size_t add_v = 0; // 情報量0のマスを追加することで増える埋蔵量
    const auto &best_layout = pool[0];
    for (int _ = 0; _ < 3; ++_)
    {
        bool change = false;
        // クエリに含める座標を1つずつ増やすか減らすかして、
        // 相互情報量/コストを最大化する山登り
        for (const auto &[_, ij] : query_coordinate_evals)
        {
            query.flip(ij);
            double eval = query.eval(sim, add_k, add_v);

            if (crt < eval)
            {
                crt = eval;
                change = true;
            }
            else
            {
                query.flip(ij);
            }
        }
        {
            // 情報量0のマスの追加
            while (add_k < no_info_coordinates.size())
            {
                double eval = query.eval(sim, add_k + 1, add_v + best_layout.volume[no_info_coordinates[add_k]]);
                if (crt < eval)
                {
                    crt = eval;
                    add_v += pool[0].volume[no_info_coordinates[add_k]];
                    add_k += 1;
                    change = true;
                }
                else
                {
                    break;
                }
            }
            while (add_k > 0)
            {
                double eval = query.eval(sim, add_k - 1, add_v - best_layout.volume[no_info_coordinates[add_k - 1]]);
                if (crt < eval)
                {
                    crt = eval;
                    add_v -= pool[0].volume[no_info_coordinates[add_k - 1]];
                    add_k -= 1;
                    change = true;
                }
                else
                {
                    break;
                }
            }
        }

        if (!change)
        {
            // どの座標を1点変えても相互情報量/コストが変わらない場合は終了
            break;
        }
    }
    vector<size_t> query_coordinates = query.get_coordinates();
    query_coordinates.insert(query_coordinates.end(), no_info_coordinates.begin(), no_info_coordinates.begin() + add_k);
    return query_coordinates;
}

void sort_pool(vector<OilLayout> &pool)
{
    sort(pool.begin(), pool.end(), [](const OilLayout &a, const OilLayout &b)
         { return b.ln_pR_if_x < a.ln_pR_if_x; });
}

int main()
{
    Input input = read_input();

    Sim sim(input);
    State state(input);
    vector<OilLayout> pool;
    size_t ITER = 4'000'000 / (2 * input.n2);
    for (size_t t = 0;; ++t)
    {
        if (sim.rem == 0)
        {
            cerr << "!There is no more query" << endl;
            break;
        }
        if (get_time() > 2.9)
        {
            sim.giveup();
            break;
        }

        // プールに存在する配置と対数尤度を記録しておく
        unordered_map<uint64_t, double> hash_ln_lilelihood;
        for (const auto &layout : pool)
        {
            hash_ln_lilelihood[layout.hash] = layout.ln_pR_if_x;
        }
        // ランダムに候補を生成する
        for (size_t _ = 0; _ < ITER; ++_)
        {
            for (size_t oil_id = 0; oil_id < input.m; ++oil_id)
            {
                const auto &oil = input.oils[oil_id];
                size_t i = rng.randrange(input.n - oil.max_i);
                size_t j = rng.randrange(input.n - oil.max_j);
                state.move_to(oil_id, i * input.n + j);
            }
            // まだプールに追加していない配置なら追加する
            if (!hash_ln_lilelihood.count(state.hash))
            {
                hash_ln_lilelihood[state.hash] = 0.0;
                pool.emplace_back(
                    OilLayout{state.hash, 0.0, 0.0, state.top_lefts, state.volumes});
            }
        }

        // 各配置の対数尤度を更新する
        // t=0のときはpoolになにも入っていないので何もしない
        for (auto &layout : pool)
        {
            if (layout.volume.empty() && !sim.failed.empty())
            {
                layout.volume = input.get_volume(layout.top_lefts);
            }
            layout.ln_pR_if_x = sim.get_ln_pR_if_x(state.oil_states, layout.volume, layout.top_lefts);
        }
        // 同じ尤度の配置を散らすためにシャッフル
        shuffle(pool.begin(), pool.end(), rng);
        // 対数尤度が高い順に配置候補をソート
        sort_pool(pool);

        // この時点でpoolはソート済み
        auto max_prob = pool[0].ln_pR_if_x;
        // 対数尤度から尤度に戻す
        // p(R|x) = exp(log(p(R|x)))
        for (auto &layout : pool)
        {
            layout.px_if_R = exp(layout.ln_pR_if_x - max_prob);
        }

        // 尤度から事後確率に変換する
        // p(x|R) = p(R|x)/Σp(R|x)
        normalize_pool(pool);

        // 尤度が低すぎる配置を削除
        while (pool.size() > 1 && pool[pool.size() - 1].px_if_R < 1e-9)
        {
            pool.pop_back();
        }
        const auto &best_layout = pool[0];
        auto best_bits = input.get_positives(best_layout.top_lefts);

        double best_pool_prob = best_layout.px_if_R;
        concat_pool(pool, input, ITER);
        // poolを切り取ったことで合計100%じゃなくなったので再度正規化
        normalize_pool(pool);
        // 油田配置の全座標の埋蔵量を計算
        set_volume(pool, input);

        if (best_pool_prob > 0.9)
        {
            // 自信があるとき、推測をう
            vector<size_t> T_vec;
            for (size_t ij = 0; ij < input.n2; ++ij)
            {
                if (best_bits[ij])
                {
                    T_vec.push_back(ij);
                }
            }
            if (sim.ans(T_vec))
            {
                break;
            }
            else if (sim.failed.size() == 1)
            {
                state.volumes = input.get_volume(state.top_lefts);
            }
        }
        else
        {
            auto query_coordinates = getDivinationQuery(input, pool, sim);
            // 占いの評価値の方が高い場合、占いを行う
            sim.query(query_coordinates);
            state.add_query(query_coordinates);
        }
    }

    cerr << "!Time = " << get_time() << endl;
    cerr << "!log miss " << sim.failed.size() << endl;
    cerr << "!main end" << endl;
}