//Dijkstra
vi dist(V, INF); dist[s] = 0;
priority_queue< pair<ll,ll>, vector<pair<ll,ll> >, greater<pair<ll,ll> > > pq; pq.push(make_pair(0, s));
while (!pq.empty())
{
    pair<ll,ll> front = pq.top(); pq.pop(); // greedy: get shortest unvisited vertex
    int d = front.first, u = front.second;
    if (d > dist[u]) continue;
    for (int j = 0; j < (int)AdjList[u].size(); j++)
    {
        pair<ll,ll> v = AdjList[u][j];
        if (dist[u] + v.second < dist[v.first])
        {
            dist[v.first] = dist[u] + v.second; // relax operation
            pq.push(make_pair(dist[v.first], v.first));
        }
    }
}
