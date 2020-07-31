//Daniel Bai 07/31/2020
//Dijkstra's shortest path algorithm
//Implemented using connectivity matrix
//For undirected graph (edge matrix symmetric)
//For ease of manual checking, the weight type was set to int
//Graph is randomly generated using 4 parameters: number of nodes, edge density, min and max of edge weights
//The algorithm stops as soon as the destination is found, thus only have nodes between source and destination nodes
//It can be modified to continue running to generate the entire shortest path tree from a given source node

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <climits>
#include <queue>
#include <tuple>  //key methods: make_tuple, get<i>
#include <cstddef>
#include <chrono>
using namespace std;

ostream& operator<< (ostream& out, vector<int> vec){
  for(int i=0; i<vec.size(); i++)
    cout << vec[i] << ", ";
  cout << endl;
  return out;
}


//function returning a random number between 0 and 1
double prob()
{ //must do seed outside of this func and only once, else it's not random!
  return static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
}

//return index of Min element of an array
int min_index(int* cost, int size)
{
  int min=INT_MAX;
  int minIdx;
  for(int i=0; i<size; i++){
    if(cost[i]<min)
      min = cost[i];
      minIdx = i;
  }
  return minIdx;
}

class graph{
public:
 graph(int V, double density, int dLow, int dHigh);
 ~graph();
 int get_num_node(){return V;}
 double get_density(){return density;}
 int get_dLow(){return dLow;}
 int get_dHigh(){return dHigh;}
 int get_num_edge();
 bool adjacent(int x, int y){if(edge[x][y]<get_dLow()) return true; else return false;}  //true if there is an edge between nodes x and y, false otherwise.
 int get_edge_value(int x, int y){return edge[x][y];}
 void set_edge_value(int x, int y, int val){edge[x][y] = edge[y][x] = val;}
 void add_edge(int x, int y, int weight);
 void delete_edge(int x, int y)
 {
   if(edge[x][y] != -1) edge[x][y] = edge[y][x] = -1;
   else cout << "No eddge exists between node " << x << " and " << y << endl;
 }  //remove an edge between nodes x and y
 vector<int> neighbors(int x);  //list of neighbors of node x
 void showAdjMatrix();
 int get_node_value(int x){return keys[x];}  //get key value of a node
 void set_node_value(int x, int a){keys[x] = a;}  //set key value of a node

private:
 int V;  //num nodes, numbered 0 to V-1
 int E;  //num edges
 double density;  //edge density, range 0.2 and 0.4 for this assignment
 int dLow, dHigh;  //edge distance range 1 to 10 for this assignment
 int** edge;  //edge weight matrix, 2D array size VxV, thus 2 stars
 int* keys;  //key values for each node, array of size V, assuming nodes are labeled as integers starting 0
};

//graph constructor
graph::graph(int V=0, double density=0.0, int dLow=1, int dHigh=10):
   V(V), density(density), dLow(dLow),dHigh(dHigh)
{
 //edge density [0.0,1.0], 0.2 and 0.4 for this assignment. Edge distance range: low, high [1, 10] for this assignment)

 const int size = V;
 if(size == 0) return;  //no need to populate edge matrix

 keys = new int [size]; //init keys
 for(int i=0; i<size; i++){
   keys[i] = i;
 }

 edge = new int* [size];  //heap allocation
 for(int i=0; i<size; i++){
   edge[i] = new int [size];
 }
 //microseconds for rand seed. time() returns second, no good - if running many runs they all end up with same seconds thus same seed!
 uint64_t us = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::
                  now().time_since_epoch()).count();
 srand(us);  //random seeed

 for(int i=0; i<size; i++){  //V: num Nodes of the graph, graph private member
   for(int j=0; j<size; j++){
      if(i==j)
       edge[i][j] = 0;
      else{
     //randomly initialize edges (0: no connection, else: random value between dLow and dHigh
        if(prob() < density){
          int tmp = static_cast<int>(dLow+prob()*(dHigh-dLow));
          edge[i][j] = edge[j][i] = tmp;
        }
        else
          edge[i][j] = edge [j][i] = 100000;  //-1 denotes no connection
       //cout << edge[i][j] << ",";
       //cout << get_edge_value(i,j) << ", ";
     }
   }
 }
}

int graph::get_num_edge(){
 const int size = V;
 int numEdge = 0;
 for(int i=0; i<size; i++){
   for(int j=0; j<size; j++){
     if(edge[i][j]>=get_dLow() && edge[i][j]<=get_dHigh()) ++numEdge;
   }
 }
 return numEdge/2;
}

void graph::add_edge(int x, int y, int weight){
 if(x==y){
   cout << "Can't add edge to self!" << endl;
   return;
 }
 if(edge[x][y] == -1)
   edge[x][y] = edge[y][x] = weight;
 else
   cout<<"Edge already exists between nodes " << x << " and " << y << "!" << endl;
 } //add an edge between nodes x and y

vector<int> graph::neighbors(int x){  //return all neighbors of node x in a vector
 vector<int> neighborList;
 const int size = V;
 for(int i=0; i<size; i++){
   if(x != i){  //excluding self
     if(edge[x][i]>=get_dLow() && edge[x][i]<=get_dHigh())
       neighborList.push_back(i);
   }
 }
 return neighborList;
}

void graph::showAdjMatrix(){
  const int size = V;
  for(int i=0; i<size; i++){
    for(int j=0; j<size; j++){
      cout << get_edge_value(i,j) << ", ";
    }
    cout << endl;
  }
}


graph::~graph(){  //destructor
 const int size = V;
// for(int i=0; i<size; i++){
//   delete [] edge[i];  //delete only takes pointer
// }
//delete [] edge;
}


//Below is class PriorityQueue
typedef pair<int, int> cost_node;  //cost_node is a type name of a pair, 1st is cost, 2nd is node ID
//typedef priority_queue<cost_node, vector<cost_node>, greater<cost_node>> cost_node_PQ;

class PriorityQueue{  //augmented from STL priority_queue, the built-in methods top(), push(), pop() still usable
  public:
    PriorityQueue(){}  //constructor
    int size();
    priority_queue<cost_node, vector<cost_node>, greater<cost_node>>& getPQ(){return pq;}  //return reference so inplace mutation can be done to pq
    bool contains(int node);
    void updatePriority(cost_node pair);  //update the priority queue if a node's cost is changed
    void push(cost_node pair){pq.push(pair);}
    void pop(){pq.pop();}
    bool empty(){return pq.empty();}
    cost_node top(){return pq.top();}
    void show();

  private:
    priority_queue<cost_node, vector<cost_node>, greater<cost_node>> pq;
};

int PriorityQueue::size(){  //return size of pq while keeping pq intact
  int s=0;
  cost_node tmpPair;
  vector<cost_node> tmpPairArr;

  while(!pq.empty()){
    int fst = (pq.top()).first;
    int snd = (pq.top()).second;
    tmpPairArr.push_back(make_pair(fst, snd));
    pq.pop();
  }
  //recover pq, or it will mutate pq
  for(int i=0; i<tmpPairArr.size(); i++)
    pq.push(tmpPairArr[i]);

  return tmpPairArr.size();
}

//return true if the pair with the specified node is in the pq
bool PriorityQueue::contains(int node){
  cost_node tmpPair;
  vector<cost_node> tmpPairArr;
  bool flag = false;

  while(!pq.empty()){
    int fst = (pq.top()).first;
    int snd = (pq.top()).second;
    tmpPairArr.push_back(make_pair(fst, snd));
    pq.pop();
    if(snd == node){
     flag=true;
     break;
    }
  }
  //recover pq, or it will mutate pq
  for(int i=0; i<tmpPairArr.size(); i++)
    pq.push(tmpPairArr[i]);

  return flag;
}

void PriorityQueue::updatePriority(cost_node pair){  //update the priority queue with the fed pair (cost, node), i.e. update cost for the node
  cost_node tmpPair;
  vector<cost_node> tmpPairArr;

  while(!pq.empty()){
    int fst = (pq.top()).first;
    int snd = (pq.top()).second;
    if(snd == pair.second){
      fst = pair.first;
      snd = pair.second;
     }
    tmpPairArr.push_back(make_pair(fst, snd));
    pq.pop();
    if(snd == pair.second)
      break;

  }
  //recover pq, or it will mutate pq
  for(int i=0; i<tmpPairArr.size(); i++)
    pq.push(tmpPairArr[i]);
}

void PriorityQueue::show(){
  cost_node tmpPair;
  vector<cost_node> tmpPairArr;

  if(pq.empty()){
    cout << "The priority queue is empty!" << endl;
    return;
  }
  cout << "PriorityQueue:" << endl;
  while(!pq.empty()){
    int fst = (pq.top()).first;
    int snd = (pq.top()).second;
    cout << "(" << fst << ", " << snd << ")" << endl;
    pq.pop();
    tmpPairArr.push_back(make_pair(fst, snd));
  }
  //recover pq, or it will mutate pq
  for(int i=0; i<tmpPairArr.size(); i++)
    pq.push(tmpPairArr[i]);
}

class ShortestPath{  //find shortest path for a given undirected graph g
  public:
    ShortestPath(graph& g, int u, int w):g(g), u(u), w(w){}
    graph& getGraph(){return g;}
    int getSrcNode(){return u;}
    int getDestNode(){return w;}
    vector<int> path(int u, int w);
    vector<int> vertices();  //call graph g by reference (when calling, just use g, not &g)
    void showPath(const vector<int>& path);
    void setTotCost(int cost){totCost = cost;}
    int getTotCost(){return totCost;}
    void showCost();
    int getNumSteps(){return numSteps;}
    void setNumSteps(int nSteps){numSteps = nSteps;}

  private:
    graph g;
    int u, w;
    int totCost; //cost from u to w shortest path
    int numSteps; //total steps from u to w shortest path

};

//returns the list of vertices of a graph g
vector<int> ShortestPath::vertices(){
  vector<int> nodes;
  for(int i=0; i<getGraph().get_num_node(); i++){
    nodes.push_back(i);
  }
  return nodes;
}

//find shortest path between source node u and destination node w and returns the sequence of vertices representing shorest path u-v1-v2-…-vn-w.
vector<int> ShortestPath::path(int srcNode, int destNode){  //u: source
  const int numV = getGraph().get_num_node();  //num vertices
  PriorityQueue openSet ;  //for nodes being considered but not finalized yet. See PriorityQueue class
  vector<cost_node> closeSet;
  int closeSetArr[numV];  //flag for nodes in closeSet
  int costList[numV];  //cost of each node to source
  int prev[numV];  //prev of each node along shortest path
  for(int i=0; i<numV; i++){  //init
    costList[i] = 100000;
    closeSetArr[i] = false;
    prev[i] = -1;  //initial no parents. When building the path later, need to ensure idx of prev != -1
    if(i == srcNode)
      costList[i]=0;
  }

//populate openSet with source node
  int currNode = srcNode;
  int tmpCost = costList[currNode];
  cost_node currPair = make_pair(tmpCost, currNode);
  openSet.push(currPair);

  while(!openSet.empty()){
    //openSet.show();
    cost_node topPair = openSet.top();
    tmpCost = topPair.first;
    currNode = topPair.second;
    closeSet.push_back(topPair);

    closeSetArr[currNode] = true;  //flag true if a node is put in closeSet
    openSet.pop();
    if(currNode == destNode) break;  //reached w, done

    vector<int> nnList = getGraph().neighbors(currNode);  //get all nearest neighbors of the node that was just placed in closeSet
    for(int i=0; i<nnList.size(); i++){  //find shortest path from currNode to its NN
      int newNode = nnList[i];
      int parent = currNode;
      if(closeSetArr[newNode]) continue;  //if the node is already in closeSet skip it
      int currCost = getGraph().get_edge_value(currNode, newNode);
      cost_node newPair;

       if(tmpCost + currCost < costList[newNode]){
        costList[newNode] = tmpCost + currCost;  //relaxation
        prev[newNode] = parent;
        newPair = make_pair(costList[newNode], newNode);

        if(!(openSet.contains(newNode)))
          openSet.push(newPair);
        else
          openSet.updatePriority(newPair);
       }
    }
  }

  vector<int> returnedPath;  //for function return
  if(closeSet[closeSet.size()-1].second != destNode){  //if w was found it should be the last of closeSet
    cout << "No path found between node " << srcNode << " and node " << destNode << "!" << endl;
    vector<int> emptyVec;
    return emptyVec;  //no path found
  }
  else{
    vector<int> path;
    currNode = destNode;
    path.push_back(currNode);
    while(prev[currNode]!=-1 && prev[currNode] <= numV){  //-1 indicates no parents
      path.push_back(prev[currNode]);
      currNode = prev[currNode];
    }
    //reverse order of path to src to dest
    for(int i=path.size()-1; i>=0; i--)
      returnedPath.push_back(path[i]);
    //save total cost
    setTotCost(costList[destNode]);
    setNumSteps(path.size()-1);
  }
  return returnedPath;
}

void ShortestPath::showPath(const vector<int>& path){
  cout << "The shortest path from node " << getSrcNode() << " to node " << getDestNode() << " is: " << endl;
  for(int i=0; i<path.size(); i++){
    cout << path[i];
    if(i!=path.size()-1)
      cout << "->";
  }
}

void ShortestPath::showCost(){
  int cost;
  cost = getTotCost();
  cout << "\nTotal cost from node " << getSrcNode() << " to node " << getDestNode() << ": " << cost << endl;

}

int main()
{
  int nNode=50, minW=1, maxW=10;
  double den=0.4;
  int srcNode=0, destNode=10;
  //cout << "Enter 4 numbers separated by space: NumNode (0-1000), density(0-1), min edge weight (e.g. 1), max edge weight (<1000)... " << endl;
  //cin >> nNode >> den >> minW >> maxW;
  //cout << "Enter 2 numbers, source node and destination node: (0-NumNode-1)" << endl;
  //cin >> srcNode >> destNode;

 //g1.showAdjMatrix();

 //Run Dijkstra shortest path
  //cout << "Num nodes, density, dLow (min weight), dHigh (max weight): " << endl;
  //cout << g1.get_num_node()<< ", "<< g1.get_density()<<", "<< g1.get_dLow()<<", "<< g1.get_dHigh()<<endl;
  //cout << "Number of edges: " << g1.get_num_edge() << endl;

  double avgSteps=0.0, avgCost = 0.0;
  int numRuns = 100;
  for(int i=0; i<numRuns; i++){

    graph g1(nNode, den, minW, maxW);
    const int size = g1.get_num_node();
    ShortestPath sp1(g1, srcNode, destNode);
    vector<int> path = sp1.path(srcNode, destNode);
    sp1.showPath(path);
    sp1.showCost();
    avgCost += sp1.getTotCost();
    avgSteps += sp1.getNumSteps();
  }

  avgCost = avgCost / numRuns;
  avgSteps = avgSteps / numRuns;
  cout << "Average cost is: " << avgCost << endl;
  cout << "Average steps is: " << avgSteps << endl;

  //cout << "Nodes of the graph: " << endl;
  //cout<< sp1.vertices();

/* BELOW ARE TEST CODES FOR DEBUGGING
 cout << "Adjacency matrix: " << endl;
 const int size = g1.get_num_node();
 cout << "Num edges: " << g1.get_num_edge() << endl;
 g1.showAdjMatrix();

 PriorityQueue pq1;
 pq1.push(make_pair(1,1));
 pq1.push(make_pair(2,2));
 pq1.push(make_pair(3,3));
 pq1.push(make_pair(4,5));
 cout << "pq1 size is: " << pq1.size() << endl;
 int node = 6;
 cout << "pq1 contains node " << node << "? " << pq1.contains(node) << endl;
 pq1.show();

 //add an edge
 g1.add_edge(3,4,6);
 cout << "Num edges: " << g1.get_num_edge() << endl;
 g1.showAdjMatrix();

 //delete an edge
 g1.delete_edge(3,4);
 cout << "Num edges: " << g1.get_num_edge() << endl;
 g1.showAdjMatrix();

//print neighbor lists
 for(int i=0; i<size; i++){
   vector<int> neighborList;
   neighborList = g1.neighbors(i);
   cout << "Neighbors of node " << i << ": ";
   for(int j=0; j<neighborList.size(); j++){
     cout << neighborList[j] << ", ";
   }
   cout << endl;
 }

//print edge values
 cout << "Edge values: ";
 for(int i=0; i<size; i++){
   for(int j=0; j<size; j++){
     cout << "(" << i << ", " << j << "): "<< g1.get_edge_value(i,j) << endl;
   }
 }
*/

}
