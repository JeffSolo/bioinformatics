import assert from 'assert';
import Graph from '../bioinformatics/Graph';

const  initialValues = {
  Node1: [],
  Node2: ['Node1'],
  3: [],
  4: [3, 'Node2'],
  5: [3, 4]
};

describe('Graph class', () => {
  it('should initialize graph', () => {
    const graph = new Graph(initialValues);
    const blankGraph = new Graph();

    assert.deepStrictEqual(graph.adjacencyList, initialValues);
    assert.deepStrictEqual(blankGraph.adjacencyList, {});
  });

  it('should add nodes', () => {
    const expected = {
      Node1: [],
      Node2: ['Node1'],
      3: [],
      4: [3, 'Node2']
    };

    const graph = new Graph();
    graph.addNode('Node1');
    graph.addNode('Node2', ['Node1']);
    graph.addNode(3);
    graph.addNode(4, [3, 'Node2']);

    assert.deepStrictEqual(graph.adjacencyList, expected);
  });

  it('should remove nodes', () => {
    const expectedFirst =  {
      Node1: [],
      Node2: ['Node1'],
      4: [3, 'Node2'],
      5: [3, 4]
    };
    const expectedSecond =  {
      Node1: []
    };
    const graph = new Graph(initialValues);

    graph.removeNode(3);

    assert.deepStrictEqual(graph.adjacencyList, expectedFirst);

    graph.removeNode(5);
    graph.removeNode('Node2');
    graph.removeNode(4);

    assert.deepStrictEqual(graph.adjacencyList, expectedSecond);
  });

  it('should add edges');

  it('should add nodes and edges');

  it('should throw error for edges pointing to non-existent nodes');

  it('should remove edges');

  // TODO - decide how to handle duplications
  it('should handle duplicate nodes');
  it('should handle duplicate edges');

  // TODO - decide how to handle nonexistent removals - probably just console.warn()?
  it('should handle removal of nonexistent node');
  it('should handle removal of nonexistent edges');
});
