#include <bits/stdc++.h>

using namespace std;
#define ll long long
#define vll vector<ll>

class SegmentTree { // the segment tree is stored like a heap array
private:
    ll n;
    vll st, A, lazy;
    ll left (ll p) { return p << 1; } // same as binary heap operations
    ll right(ll p) { return (p << 1) + 1; }

    void build(ll p, ll L, ll R) { // O(n)
        if (L == R) // as L == R, either one is fine
            st[p] = A[L]; // store the index
        else { // recursively compute the values
            build(left(p) , L , (L + R) / 2);
            build(right(p), (L + R) / 2 + 1, R );
            ll p1 = st[left(p)], p2 = st[right(p)];
            st[p] = max(p1,p2);
        }
    }
    void push(int p) {
        st[p*2] += lazy[p];
        lazy[p*2] += lazy[p];
        st[p*2+1] += lazy[p];
        lazy[p*2+1] += lazy[p];
        lazy[p] = 0;
}
    ll rmq(ll p, ll L, ll R, ll i, ll j) { // O(log n)
        if (i > j) return LLONG_MIN;
        if(i==L && j==R) return st[p];
        push(p);
        ll m=(L+R)/2;
        ll p1 = rmq(left(p) , L , (L+R) / 2, i, min(j,m));
        ll p2 = rmq(right(p), (L+R) / 2 + 1, R , max(i,m+1), j);
        return max(p1,p2); // as in build routine
    }
    void update(ll p, ll L, ll R, ll l, ll r, ll addend){
        if (l > r)
        return;
        if (l == L && r == R) {
            st[p] += addend;
            lazy[p] += addend;
        } else {
            push(p);
            ll m = (L + R) / 2;
            update(left(p), L, m, l, min(r, m), addend);
            update(right(p), m+1, R, max(l, m+1), r, addend);
            st[p] = max(st[left(p)], st[right(p)]);
        }
    }

public:
    SegmentTree(const vll &_A) {
        A = _A; n = (int)A.size(); // copy content for local usage
        st.assign(4 * n, 0); // create large enough vector of zeroes
        lazy.assign(4*n,0);
        build(1, 0, n - 1); // recursive build
    }
    ll rmq(ll i, ll j) { return rmq(1, 0, n - 1, i, j); } // overloading
    void update(ll l, ll r, ll addend) {return update(1,0,n-1,l,r,addend); }
};


