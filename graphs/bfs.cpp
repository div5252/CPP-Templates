// inside int main()---no recursion
vll d(V, INF); d[s] = 0;
queue<int> q; q.push(s);
while (!q.empty())
{
    int u = q.front(); q.pop(); // queue: layer by layer!
    for (int j = 0; j < (int)AdjList[u].size(); j++)
    {
        ii v = AdjList[u][j]; // for each neighbor of u
        if (d[v.first] == INF)
        {
            d[v.first] = d[u] + 1;
            q.push(v.first);
        }
    }
}
