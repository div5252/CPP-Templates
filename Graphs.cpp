//DFS
void dfs(int u) {
    dfs_num[u] = 1;
    for (int j = 0; j < (int)AdjList[u].size(); j++) {
        ll v = AdjList[u][j];
        if (dfs_num[v] == 0)
            dfs(v);
    }
}


//BFS
// inside int main()---no recursion
vll d(V, INF);
d[s] = 0;
queue<ll> q;
q.push(s);
while (!q.empty())
{
    ll u = q.front(); q.pop(); // queue: layer by layer!
    for (int j = 0; j < (int)AdjList[u].size(); j++)
    {
        pll v = AdjList[u][j]; // for each neighbor of u
        if (d[v.first] == INF)
        {
            d[v.first] = d[u] + 1;
            q.push(v.first);
        }
    }
}


//Connected component-Undirected
// inside int main()---this is the DFS solution
numCC = 0;
dfs_num.assign(V, UNVISITED); // sets all vertices’ state to UNVISITED
for (int i = 0; i < V; i++) // for each vertex i in [0..V-1]
    if (dfs_num[i] == UNVISITED) // if vertex i is not visited yet
        printf("CC %d:", ++numCC), dfs(i), printf("\n"); // 3 lines here!


//Topological sort
vi ts; // global vector to store the toposort in reverse order
void dfs2(int u) {
    dfs_num[u] = VISITED;
    for (int j = 0; j < (int)AdjList[u].size(); j++) {
        ii v = AdjList[u][j];
        if (dfs_num[v.first] == UNVISITED)
            dfs2(v.first);
    }
    ts.push_back(u);
}


//Kahn Algo-TS
enqueue vertices with zero incoming degree into a (priority) queue Q;
while (Q is not empty) {
    vertex u = Q.dequeue(); put vertex u into a topological sort list;
    remove this vertex u and all outgoing edges from this vertex;
    if such removal causes vertex v to have zero incoming degree
        Q.enqueue(v); }


// To check graph is bipartite
// inside int main()
queue<int> q; q.push(s);
vi color(V, INF); color[s] = 0;
bool isBipartite = true;
while (!q.empty() && isBipartite)
{
    int u = q.front(); q.pop();
    for (int j = 0; j < (int)AdjList[u].size(); j++) {
        ll v = AdjList[u][j];
        if (color[v.first] == INF) {
            color[v.first] = 1 - color[u];
            q.push(v.first); }
        else if (color[v.first] == color[u]) {
            isBipartite = false; break; 
        } 
    } 
}


//Graph edges property
void graphCheck(int u) { // DFS for checking graph edge properties
    dfs_num[u] = EXPLORED; // color u as EXPLORED instead of VISITED
    for (int j = 0; j < (int)AdjList[u].size(); j++) {
        ii v = AdjList[u][j];
        if (dfs_num[v.first] == UNVISITED) { // Tree Edge, EXPLORED->UNVISITED
            dfs_parent[v.first] = u; // parent of this children is me
            graphCheck(v.first);
        }
        else if (dfs_num[v.first] == EXPLORED) { // EXPLORED->EXPLORED
            if (v.first == dfs_parent[u]) // to differentiate these two cases
                printf(" Two ways (%d, %d)-(%d, %d)\n", u, v.first, v.first, u);
            else // the most frequent application: check if the graph is cyclic
                printf(" Back Edge (%d, %d) (Cycle)\n", u, v.first);
        }
        else if (dfs_num[v.first] == VISITED) // EXPLORED->VISITED
            printf(" Forward/Cross Edge (%d, %d)\n", u, v.first);
    }
    dfs_num[u] = VISITED; // after recursion, color u as VISITED (DONE)
}
// inside int main()
    dfs_num.assign(V, UNVISITED);
    dfs_parent.assign(V, 0); // new vector


//Articulation Point and Bridges
void articulationPointAndBridge(int u) {
    dfs_low[u] = dfs_num[u] = dfsNumberCounter++; // dfs_low[u] <= dfs_num[u]
    for (int j = 0; j < (int)AdjList[u].size(); j++) {
        ii v = AdjList[u][j];
        if (dfs_num[v.first] == UNVISITED) { // a tree edge
            dfs_parent[v.first] = u;
            if (u == dfsRoot) rootChildren++; // special case if u is a root
            articulationPointAndBridge(v.first);
            if (dfs_low[v.first] >= dfs_num[u] && u!=dfsRoot) // for articulation point
                articulation_vertex[u] = true; // store this information first
            if (dfs_low[v.first] > dfs_num[u]) // for bridge
                printf(" Edge (%d, %d) is a bridge\n", u, v.first);
            dfs_low[u] = min(dfs_low[u], dfs_low[v.first]); // update dfs_low[u]
        }
        else if (v.first != dfs_parent[u]) // a back edge and not direct cycle
            dfs_low[u] = min(dfs_low[u], dfs_num[v.first]); // update dfs_low[u]
} }
// inside int main()
    dfsNumberCounter = 0; dfs_num.assign(V, UNVISITED); dfs_low.assign(V, 0);
    dfs_parent.assign(V, 0); articulation_vertex.assign(V, 0);
    printf("Bridges:\n");
    for (int i = 0; i < V; i++)
        if (dfs_num[i] == UNVISITED) {
            dfsRoot = i; rootChildren = 0; articulationPointAndBridge(i);
            articulation_vertex[dfsRoot] = (rootChildren > 1); } // special case
    printf("Articulation Points:\n");
    for (int i = 0; i < V; i++)
        if (articulation_vertex[i])
            printf(" Vertex %d\n", i);


//SCC
vi dfs_num, dfs_low, S, visited; // global variables
void tarjanSCC(int u) {
    dfs_low[u] = dfs_num[u] = dfsNumberCounter++; // dfs_low[u] <= dfs_num[u]
    S.push_back(u); // stores u in a vector based on order of visitation
    visited[u] = 1;
    for (int j = 0; j < (int)AdjList[u].size(); j++) {
        ii v = AdjList[u][j];
        if (dfs_num[v.first] == UNVISITED)
            tarjanSCC(v.first);
        if (visited[v.first])   // condition for update
            dfs_low[u] = min(dfs_low[u], dfs_low[v.first]); 
    }

    if (dfs_low[u] == dfs_num[u]) {     // if this is a root (start) of an SCC
        printf("SCC %d:", ++numSCC);    // this part is done after recursion
        while (1) {
            int v = S.back(); S.pop_back(); visited[v] = 0;
            printf(" %d", v);
            if (u == v) break; 
        }
        printf("\n");
    } 
}
// inside int main()
    dfs_num.assign(V, UNVISITED); dfs_low.assign(V, 0); visited.assign(V, 0);
    dfsNumberCounter = numSCC = 0;
    for (int i = 0; i < V; i++)
        if (dfs_num[i] == UNVISITED)
            tarjanSCC(i);


//Kruskal-MST
// inside int main()
vector< pair<int, ii> > EdgeList; // (weight, two vertices) of the edge
for (int i = 0; i < E; i++) {
    scanf("%d %d %d", &u, &v, &w); // read the triple: (u, v, w)
    EdgeList.push_back(make_pair(w, ii(u, v))); // (w, u, v)
} 
sort(EdgeList.begin(), EdgeList.end());
int mst_cost = 0;
UnionFind UF(V); // all V are disjoint sets initially
for (int i = 0; i < E; i++) { // for each edge, O(E)
    pair<int, ii> front = EdgeList[i];
    if (!UF.isSameSet(front.second.first, front.second.second)) { // check
        mst_cost += front.first; // add the weight of e to MST
        UF.unionSet(front.second.first, front.second.second); // link them
    } 
}


//Prim-MST
vi taken; // global boolean flag to avoid cycle
priority_queue<ii> pq; // priority queue to help choose shorter edges
void process(int vtx) {
    taken[vtx] = 1;
    for (int j = 0; j < (int)AdjList[vtx].size(); j++) {
        ii v = AdjList[vtx][j];
        if (!taken[v.first]) pq.push(ii(-v.second, -v.first));
} }
// inside int main()---assume the graph is stored in AdjList, pq is empty
    taken.assign(V, 0); // no vertex is taken at the beginning
    process(0); // take vertex 0 and process all edges incident to vertex 0
    mst_cost = 0;
    while (!pq.empty()) { // repeat until V vertices (E=V-1 edges) are taken
        ii front = pq.top(); pq.pop();
        u = -front.second, w = -front.first; // negate the id and weight again
        if (!taken[u]) // we have not connected this vertex yet
            mst_cost += w, process(u); // take u, process all edges incident to u
    } // each edge is in pq only once!
    printf("MST cost = %d\n", mst_cost);


//For Minimax using mst
void dfs(ll u, ll v_end, ll max_val)
{
    if (u == v_end)
        max_w = max_val;
    visited[u] = 1;
    for (int i = 0; i < AdjList[u].size(); i++)
    {
        auto v = AdjList[u][i];
        if (!visited[v.first])
            dfs(v.first, v_end, max(max_val, v.second));
    }
}


//SSSP on Unweighted graph
vi dist(V, INF); dist[s] = 0; // distance from source s to s is 0
queue<int> q; q.push(s);
vi p; // addition: the predecessor/parent vector
while (!q.empty()) {
    int u = q.front(); q.pop();
    for (int j = 0; j < (int)AdjList[u].size(); j++) {
        ii v = AdjList[u][j];
        if (dist[v.first] == INF) {
            dist[v.first] = dist[u] + 1;
            p[v.first] = u; // addition: the parent of vertex v.first is u
            q.push(v.first);
} } }
printPath(t), printf("\n"); // addition: call printPath from vertex t


//Dijkstra
vi dist(V, INF); 
dist[s] = 0;
priority_queue< pair<ll,ll>, vector<pair<ll,ll> >, greater<pair<ll,ll> > > pq; 
pq.push(make_pair(0, s));
while (!pq.empty())
{
    pair<ll,ll> front = pq.top(); 
    pq.pop();
    ll d = front.first, u = front.second;
    if (d > dist[u]) continue;
    for (int j = 0; j < (int)AdjList[u].size(); j++)
    {
        pair<ll,ll> v = AdjList[u][j];
        if (dist[u] + v.second < dist[v.first])
        {
            dist[v.first] = dist[u] + v.second;
            pq.push(make_pair(dist[v.first], v.first));
        }
    }
}

//Dijkstra K shortest
vll ans;
ll cnt=0;
vll dist(n, 0);
priority_queue< pair<ll,ll>, vector<pair<ll,ll> >, greater<pair<ll,ll> > > pq;
pq.push(make_pair(0, 0));
while (!pq.empty())
{
    pair<ll,ll> front = pq.top(); pq.pop();
    ll d = front.first, u = front.second;
    if (dist[u]>k) continue;
    if(u==n-1)
    {
        cnt++;
        ans.PB(d);
        if(cnt==k) break;
    }
    dist[u]++;
    for (int j = 0; j < (int)AdjList[u].size(); j++)
    {
        pair<ll,ll> v = AdjList[u][j];
        if(dist[v.first]>k) continue;
        pq.push(make_pair(d+v.second, v.first));
    }
}


//Bellman Ford
vi dist(V, INF); dist[s] = 0;
for (int i = 0; i < V - 1; i++) // relax all E edges V-1 times
    for (int u = 0; u < V; u++) // these two loops = O(E), overall O(VE)
        for (int j = 0; j < (int)AdjList[u].size(); j++) {
            pair<ll,ll> v = AdjList[u][j]; // record SP spanning here if needed
            dist[v.first] = min(dist[v.first], dist[u] + v.second); // relax
        }
bool hasNegativeCycle = false;
for (int u = 0; u < V; u++) {// one more pass to check
    for (int j = 0; j < (int)AdjList[u].size(); j++) {
        pair<ll,ll> v = AdjList[u][j];
        if (dist[v.first] > dist[u] + v.second) // if this is still possible
            hasNegativeCycle = true; // then negative cycle exists!
    }
}
printf("Negative Cycle Exist? %s\n", hasNegativeCycle ? "Yes" : "No");


//Floyd Warshall
// inside int main()
// AdjMat is a 32-bit signed integer array
for (int k = 0; k < V; k++) // remember that loop order is k->i->j
    for (int i = 0; i < V; i++)
        for (int j = 0; j < V; j++)
            AdjMat[i][j] = min(AdjMat[i][j], AdjMat[i][k] + AdjMat[k][j]);

//Printing shortest path
// let p be a 2D parent matrix, where p[i][j] is the last vertex before j
// on a shortest path from i to j, i.e. i -> ... -> p[i][j] -> j
for (int i = 0; i < V; i++)
    for (int j = 0; j < V; j++)
        p[i][j] = i; // initialize the parent matrix
for (int k = 0; k < V; k++)
    for (int i = 0; i < V; i++)
        for (int j = 0; j < V; j++) // this time, we need to use if statement
            if (AdjMat[i][k] + AdjMat[k][j] < AdjMat[i][j]) {
                AdjMat[i][j] = AdjMat[i][k] + AdjMat[k][j];
                p[i][j] = p[k][j]; // update the parent matrix
            }
// when we need to print the shortest paths, we can call the method below:
void printPath(int i, int j) {
    if (i != j) printPath(i, p[i][j]);
    printf(" %d", j);
}


//Edmonds Karp's
int res[MAX_V][MAX_V], mf, f, s, t; // global variables
vi p; // p stores the BFS spanning tree from s
void augment(int v, int minEdge) { // traverse BFS spanning tree from s->t
    if (v == s) { f = minEdge; return; } // record minEdge in a global var f
    else if (p[v] != -1) { augment(p[v], min(minEdge, res[p[v]][v]));
                           res[p[v]][v] -= f; res[v][p[v]] += f; } }
// inside int main()
    mf = 0;
    while (1) {
        f = 0;
        vi dist(MAX_V, INF); dist[s] = 0; queue<int> q; q.push(s);
        p.assign(MAX_V, -1);
        while (!q.empty()) {
            int u = q.front(); q.pop();
            if (u == t) break;
            for (int v = 0; v < MAX_V; v++) // note: this part is slow
                if (res[u][v] > 0 && dist[v] == INF)
                    dist[v] = dist[u] + 1, q.push(v), p[v] = u; // 3 lines in 1!
        }
        augment(t, INF); // find the min edge weight ‘f’ in this path, if any
        if (f == 0) break; // we cannot send any more flow (‘f’ = 0), terminate
        mf += f; // we can still send a flow, increase the max flow!
    }
    printf("%d\n", mf); // this is the max flow value


//Euler Tour
list<int> cyc; // we need list for fast insertion in the middle
void EulerTour(list<int>::iterator i, int u) {
    for (int j = 0; j < (int)AdjList[u].size(); j++) {
        ii v = AdjList[u][j];
        if (v.second) { // if this edge can still be used/not removed
            v.second = 0; // make the weight of this edge to be 0 (‘removed’)
            for (int k = 0; k < (int)AdjList[v.first].size(); k++) {
                ii uu = AdjList[v.first][k]; // remove bi-directional edge
                if (uu.first == u && uu.second) {
                    uu.second = 0;
                    break;
            } }
            EulerTour(cyc.insert(i, u), v.first);
} } }
// inside int main()
    cyc.clear();
    EulerTour(cyc.begin(), A); // cyc contains an Euler tour starting at A
    for (list<int>::iterator it = cyc.begin(); it != cyc.end(); it++)
        printf("%d\n", *it);


//Augmenting path for MCBM
vi match, vis; // global variables
int Aug(int l) { // return 1 if an augmenting path is found
    if (vis[l]) return 0; // return 0 otherwise
    vis[l] = 1;
    for (int j = 0; j < (int)AdjList[l].size(); j++) {
        int r = AdjList[l][j]; // edge weight not needed -> vector<vi> AdjList
        if (match[r] == -1 || Aug(match[r])) {
            match[r] = l; return 1; // found 1 matching
    } }
    return 0; // no matching
}
// inside int main()
    // build unweighted bipartite graph with directed edge left->right set
    int MCBM = 0;
    match.assign(V, -1); // V is the number of vertices in bipartite graph
    for (int l = 0; l < n; l++) { // n = size of the left set
        vis.assign(n, 0); // reset before each recursion
        MCBM += Aug(l);
    }
    printf("Found %d matchings\n", MCBM);
