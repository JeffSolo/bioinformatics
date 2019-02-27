type NodeType = string | number;
type EdgeType = (NodeType)[];
type AdjacencyList = {
  [key: string]: EdgeType;
};

export default class Graph {
  public adjacencyList: AdjacencyList = {};

  public addNode(node: NodeType, edges: EdgeType = []): void {
    // TODO - handle duplicate values, rather than just overwriting
    this.adjacencyList[node] = edges;
  }

  public removeNode(node: NodeType): void {
    delete this.adjacencyList[node]; //tslint:disable-line no-dynamic-delete
  }

  public addEdge(node: NodeType, edges: EdgeType): void {
    const unique = new Set([...this.adjacencyList[node], ...edges]); // easier to read on separate line
    this.adjacencyList[node] = Array.from(unique);
  }

  public removeEdge(node: NodeType, edges: EdgeType): void {
    // TODO - maybe check that edge actually exists in node?
    this.adjacencyList[node] = this.adjacencyList[node].filter((edge) => {
     return edges.indexOf(edge) === -1;
    });
  }
}
