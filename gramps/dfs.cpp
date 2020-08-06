void dfs(int u) {
    dfs_num[u] = 1;
    for (int j = 0; j < (int)AdjList[u].size(); j++) {
        ll v = AdjList[u][j];
        if (dfs_num[v] == 0)
            dfs(v);
    }
}
