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
