/*
 * A very basic graph class
 */
type NodeType = string | number;
type EdgeType = (NodeType)[];
type AdjacencyList = {
  [key: string]: EdgeType;
};

export default class Graph {
  public adjacencyList: AdjacencyList = {};

  constructor(adjacencyList: AdjacencyList = {}) {
      this.adjacencyList = adjacencyList;
  }

  public addNode(node: NodeType, edges: EdgeType = []): void {
    // TODO - handle duplicate values, rather than just overwriting
    this.adjacencyList[node] = [];
    this.addEdge(node, edges);
  }

  public removeNode(node: NodeType): void {
    // TODO - remove any edges pointing to node
    delete this.adjacencyList[node]; //tslint:disable-line no-dynamic-delete
  }

  public addEdge(node: NodeType, edges: EdgeType): void {
    // check that the edges point to existing nodes
    edges.forEach((e) => {
      if (!(e in this.adjacencyList)) {
        throw new Error(`Edge does not point to an existing node '${e}'`);
      }
    });

    const unique = new Set([...this.adjacencyList[node], ...edges]);
    this.adjacencyList[node] = Array.from(unique);
  }

  public removeEdge(node: NodeType, edges: EdgeType): void {
    // TODO - maybe check that edge actually exists in node?
    this.adjacencyList[node] = this.adjacencyList[node].filter((edge) => {
     return edges.indexOf(edge) === -1;
    });
  }
}
