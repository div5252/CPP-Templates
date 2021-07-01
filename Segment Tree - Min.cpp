#include <bits/stdc++.h>

using namespace std;
#define ll long long
#define vll vector<ll>

class SegmentTree { // the segment tree is stored like a heap array
private:
    vll st, A; // recall that vi is: typedef vector<int> vi;
    ll n;
    ll left (ll p) { return p << 1; } // same as binary heap operations
    ll right(ll p) { return (p << 1) + 1; }

    void build(ll p, ll L, ll R) { // O(n)
        if (L == R) // as L == R, either one is fine
            st[p] = L; // store the index
        else { // recursively compute the values
            build(left(p) , L , (L + R) / 2);
            build(right(p), (L + R) / 2 + 1, R );
            ll p1 = st[left(p)], p2 = st[right(p)];
            st[p] = (A[p1] <= A[p2]) ? p1 : p2;
        }
    }

    ll rmq(ll p, ll L, ll R, ll i, ll j) { // O(log n)
        if (i > j) return 0;
        if(i==L && j==R) return st[p];
        ll m=(L+R)/2;
        ll p1 = rmq(left(p) , L , (L+R) / 2, i, min(j,m));
        ll p2 = rmq(right(p), (L+R) / 2 + 1, R , max(i,m+1), j);
        return (A[p1] <= A[p2]) ? p1 : p2; // as in build routine
    }
    void update(ll p, ll L, ll R, ll pos, ll val){
        if (L == R)
            A[L] = val;
        else{
            ll m=(L+R)/2;
            if(pos<=m) update(left(p),L,m,pos,val);
            else update(right(p),m+1,R,pos,val);
            ll p1 = st[left(p)], p2 = st[right(p)];
            st[p] = (A[p1] <= A[p2]) ? p1 : p2;
        }
    }

public:
    SegmentTree(const vll &_A) {
        A = _A; n = (int)A.size(); // copy content for local usage
        st.assign(4 * n, 0); // create large enough vector of zeroes
        build(1, 0, n - 1); // recursive build
    }
    ll rmq(ll i, ll j) { return rmq(1, 0, n - 1, i, j); } // overloading
    void update(ll pos, ll val) {return update(1,0,n-1,pos,val); }
};

int main() {
    int arr[] = { 18, 17, 13, 19, 15, 11, 20 }; // the original array
    vll A(arr, arr + 7);
    SegmentTree st(A);
    printf("RMQ(1, 3) = %d\n", st.rmq(1, 3)); // answer = index 2
    printf("RMQ(4, 6) = %d\n", st.rmq(4, 6)); // answer = index 5
    st.update(1,10);
    cout<<st.rmq(0,6);
}

