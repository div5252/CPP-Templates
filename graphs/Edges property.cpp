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
