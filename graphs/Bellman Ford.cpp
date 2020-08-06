//Bellman Ford
vi dist(V, INF); dist[s] = 0;
for (int i = 0; i < V - 1; i++) // relax all E edges V-1 times
    for (int u = 0; u < V; u++) // these two loops = O(E), overall O(VE)
        for (int j = 0; j < (int)AdjList[u].size(); j++) {
            pair<ll,ll> v = AdjList[u][j]; // record SP spanning here if needed
            dist[v.first] = min(dist[v.first], dist[u] + v.second); // relax
        }
bool hasNegativeCycle = false;
for (int u = 0; u < V; u++) // one more pass to check
    for (int j = 0; j < (int)AdjList[u].size(); j++) {
        pair<ll,ll> v = AdjList[u][j];
        if (dist[v.first] > dist[u] + v.second) // if this is still possible
            hasNegativeCycle = true; // then negative cycle exists!
    }
printf("Negative Cycle Exist? %s\n", hasNegativeCycle ? "Yes" : "No");
