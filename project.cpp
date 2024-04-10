#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <queue>

using namespace std;

class Vertex {
public:
    string name;
    double x, y;

    Vertex(string name, double x, double y) : name(name), x(x), y(y) {}
};

class Edge {
public:
    string start, end;
    double distance;

    Edge(string start, string end, double distance) : start(start), end(end), distance(distance) {}
};

class Mapping{
private:
    unordered_map<string, int> name_to_index;
    unordered_map<int, string> index_to_name;
    int next_index = 0;

public:
    int addString(string name){
        if(name_to_index.find(name) == name_to_index.end()){
            name_to_index[name] = next_index;
            index_to_name[next_index] = name;
            ++next_index;
        }

        return name_to_index[name];
    }

    string getString(int index){
        return index_to_name.at(index);
    }

    int getIndex(string name){
        return name_to_index.at(name);
    }
};

class Graph_bfs {
private:
    unordered_map<string, vector<string>> region;

public:
    void addEdges(const string& u, const string& v) {
        region[u].push_back(v);
        region[v].push_back(u);
    }

    const vector<string>& neighborsOfVertex(const string& u) const {
        static const vector<string> empty;
        auto it = region.find(u);
        return (it != region.end()) ? it->second : empty;
    }

    size_t getVertexesSize() const {
        return region.size();
    }
};

class BFS{
public:
    static unordered_map<string, unordered_set<string>> bfs(const Graph_bfs& region, const vector<Vertex>& vertices){
        unordered_map<string, unordered_set<string>> reachableVertex;

        for(const auto& vertex : vertices){
            const string& start = vertex.name;
            unordered_set<string> visited;
            queue<string> q;
            q.push(start);

            while(!q.empty()){
                const string& currentVertex = q.front();
                q.pop();
                visited.insert(currentVertex);
                for(const auto& neighbor : region.neighborsOfVertex(currentVertex)){
                    if(visited.find(neighbor) == visited.end()){
                        q.push(neighbor);
                        reachableVertex[start].insert(neighbor);
                        visited.insert(neighbor);
                    }
                }
            }
        }
        return reachableVertex;
    }
};

class Graph {
private:
    vector<Vertex> vertices;
    map<string, string> parent;
    map<string, int> rank;

    double calculateDistance(Vertex v1, Vertex v2) {
        return sqrt(pow(v1.x - v2.x, 2) + pow(v1.y - v2.y, 2));
    }

    string find(string i) {
        if (parent[i] == i)
            return i;
        return find(parent[i]);
    }

    void Union(string x, string y) {
        string xroot = find(x);
        string yroot = find(y);

        if (rank[xroot] < rank[yroot])
            parent[xroot] = yroot;
        else if (rank[xroot] > rank[yroot])
            parent[yroot] = xroot;
        else {
            parent[yroot] = xroot;
            rank[xroot]++;
        }
    }

public:
    Graph() {}

    void addVertex(string name, double x, double y) {
        vertices.push_back(Vertex(name, x, y));
        parent[name] = name;
        rank[name] = 0;
    }

    void addEdge(string start, string end) {
        Union(start, end);
    }

    vector<Edge> getAllEdges() {
        vector<Edge> edges;
        for (int i = 0; i < vertices.size(); ++i) {
            for (int j = i + 1; j < vertices.size(); ++j) {
                string start = vertices[i].name;
                string end = vertices[j].name;
                if (find(start) != find(end)) {
                    double distance = calculateDistance(vertices[i], vertices[j]);
                    edges.push_back(Edge(start, end, distance));
                }
            }
        }
        return edges;
    }

    string findRoot(string vertex) {
        return find(vertex);
    }

    const vector<Vertex>& getVertices() const {
        return vertices;
    }
};

class SpanningTreeBuilder {
private:
    Graph& graph;

public:
    SpanningTreeBuilder(Graph& graph) : graph(graph) {}

    vector<Edge> createSpanningTree() {
        vector<Edge> edges = graph.getAllEdges();

        sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
            return a.distance < b.distance;
        });

        vector<Edge> spanningTree;
        for (const Edge& edge : edges) {
            if (graph.findRoot(edge.start) != graph.findRoot(edge.end)) {
                spanningTree.push_back(edge);
                graph.addEdge(edge.start, edge.end);
            }
        }

        return spanningTree;
    }
};

int main() {
    int numVertices, numConnected;
    
    cout << "Enter the number of settlements: ";
    cin >> numVertices;

    Graph graph;
    Graph_bfs graph_bfs;

    for (int i = 0; i < numVertices; ++i) {
        string name;
        double x, y;
        cout << "Enter the name and the coordinates of settlement " << i + 1 << ": ";
        cin >> name >> x >> y;
        graph.addVertex(name, x, y);
    }

    cout << "Enter the number of connections: ";
    cin >> numConnected;
    for (int i = 0; i < numConnected; ++i) {
        string start, end;
        cout << "Enter the name of the settlements in connection " << i + 1 << ": ";
        cin >> start >> end;
        graph.addEdge(start, end);
        graph_bfs.addEdges(start, end);
    }
    
    std::unordered_map<string, std::unordered_set<string>> reachableVertices = BFS::bfs(graph_bfs, graph.getVertices());
    std::cout << "Reachable Vertices:\n";
    for (const auto& entry : reachableVertices) {
        std::cout << "From vertex " << entry.first << ": ";
        for (const string& v : entry.second) {
            std::cout << v << " ";
        }
        std::cout << std::endl;
    }

    SpanningTreeBuilder builder(graph);
    vector<Edge> spanningTree = builder.createSpanningTree();
    
    cout << "Number of added edges: " << spanningTree.size() << endl;

    cout << "New added edges is: " << endl;
    for (const Edge& edge : spanningTree) {
        cout << edge.start << " " << edge.end << endl;
    }

    return 0;
}

