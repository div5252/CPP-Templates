#include <bits/stdc++.h>

using namespace std;
#define ll long long
#define vll vector<ll>

#define LSOne(S) (S & (-S))

class FenwickTree {
private:
  vll ft;

public:
  FenwickTree() {}
  // initialization: n + 1 zeroes, ignore index 0
  FenwickTree(ll n) { ft.assign(n + 1, 0); }

  ll rsq(ll b) {                                     // returns RSQ(1, b)
    ll sum = 0; for (; b; b -= LSOne(b)) sum += ft[b];
    return sum; }

  ll rsq(ll a, ll b) {                              // returns RSQ(a, b)
    return rsq(b) - (a == 1 ? 0 : rsq(a - 1)); }

  // adjusts value of the k-th element by v (v can be +ve/inc or -ve/dec)
  void adjust(ll k, ll v) {                    // note: n = ft.size() - 1
    for (; k < (int)ft.size(); k += LSOne(k)) ft[k] += v; }
};

int main() {
    int f[] = { 2,4,5,5,6,6,6,7,7,8,9 }; // m = 11 scores
    FenwickTree ft(10); // declare a Fenwick Tree for range [1..10]
    // insert these scores manually one by one into an empty Fenwick Tree
    for (int i = 0; i < 11; i++) ft.adjust(f[i], 1); // this is O(k log n)
    printf("%d\n", ft.rsq(1, 1)); // 0 => ft[1] = 0
    printf("%d\n", ft.rsq(1, 2)); // 1 => ft[2] = 1
    printf("%d\n", ft.rsq(1, 6)); // 7 => ft[6] + ft[4] = 5 + 2 = 7
    printf("%d\n", ft.rsq(1, 10)); // 11 => ft[10] + ft[8] = 1 + 10 = 11
    printf("%d\n", ft.rsq(3, 6)); // 6 => rsq(1, 6) - rsq(1, 2) = 7 - 1
    ft.adjust(5, 2); // update demo
    printf("%d\n", ft.rsq(1, 10)); // now 13
} // return 0;
