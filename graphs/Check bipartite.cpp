// To check graph is bipartite
// inside int main()
queue<int> q; q.push(s);
vi color(V, INF); color[s] = 0;
bool isBipartite = true;
while (!q.empty() && isBipartite)
{
    int u = q.front(); q.pop();
    for (int j = 0; j < (int)AdjList[u].size(); j++) {
        ii v = AdjList[u][j];
        if (color[v.first] == INF) {
            color[v.first] = 1 - color[u];
            q.push(v.first); }
        else if (color[v.first] == color[u]) {
            isBipartite = false; break; } } }
