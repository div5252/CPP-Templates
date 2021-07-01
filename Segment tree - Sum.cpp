#include <bits/stdc++.h>

using namespace std;
#define ll long long
#define vll vector<ll>

class SegmentTree { // the segment tree is stored like a heap array
private:
    vll st, A;
    ll n;
    ll left (ll p) { return p << 1; } // same as binary heap operations
    ll right(ll p) { return (p << 1) + 1; }

    void build(ll p, ll L, ll R) { // O(n)
        if (L == R)
            st[p] = A[L];
        else {
            build(left(p) , L , (L + R) / 2);
            build(right(p), ((L + R) / 2) + 1, R );
            ll p1 = st[left(p)], p2 = st[right(p)];
            st[p] = p1+p2;
        }
    }

    ll rsq(ll p, ll L, ll R, ll i, ll j) { // O(log n)
        if (i > j) return 0;
        if(i==L && j==R) return st[p];
        ll m=(L+R)/2;
        return rsq(left(p),L,m,i,min(j,m))+rsq(right(p),m+1,R,max(i,m+1),j);
    }
    void update(ll p, ll L, ll R, ll pos, ll val){
        if (L == R)
            st[p] = val;
        else{
            ll m=(L+R)/2;
            if(pos<=m) update(left(p),L,m,pos,val);
            else update(right(p),m+1,R,pos,val);
            ll p1 = st[left(p)], p2 = st[right(p)];
            st[p] = p1+p2;
        }
    }

public:
    SegmentTree(const vll &_A) {
        A = _A; n = (int)A.size(); // copy content for local usage
        st.assign(4 * n, 0); // create large enough vector of zeroes
        build(1, 0, n - 1); // recursive build
    }
    ll rsq(ll i, ll j) { return rsq(1, 0, n - 1, i, j); } // overloading
    void update(ll pos, ll val) {return update(1,0,n-1,pos,val); }
};
