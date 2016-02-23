Viziru Luciana
Readme - Sudoku

Each process first reads from an imput file the line with its own adjacency list.
It will then write the its line in the topology matrix and save the number of
neighbour nodes in the numberOfNeighbours variable.

1. runProbeEcho method - establishing the topology
The node with rank 0 is chosen as the root of the spanning tree. It will send
a probe message to itself in order to start the probe-echo stage. Then the
behaviour of all the nodes is identical.
After this stage, each node will have determined:
- topology - a topology of the graph with no cycles
- numberOfChildren- the number of children in each of its subtrees (no link -> 0 children)
- routingTable - anything not in the subtree is routed to the parent node
In the end, the root will diffuse the complete topology to all the nodes.

2. Complete sudoku
The root reads the initial matrix from an input file. The problem dimension and
the initial matrix is diffused to the entire tree.
In order to allocate a square to each node in the subtree, each node will send
the number of squares already assigned. All the nodes will count from left to right
and up to down, starting from the index 0. If n nodes, have been assigned, a node
will complete the square numbered as n.
Each node will compute all the possible solutions for its square in the
initialResults matrix using the getInitialResults function.
To avoid unnecessary memory allocation, the valid solutions are counted and then
the needed space is allocated and populated - countInitialSolutions() and
saveInitialSolutions().
To combine multiple squares, from each son, the number of valid solutions are
received, then the number of possible valid solution are allocated and only the
the valid solutions are saved - countSolutions() and mergeSolutions().
The new, merged matrixes are in the newResults matrix. After all the data is
received from a child, newResults will replace the previous initialResults.
I chose the overhead of the extra computing in order to avoid using very much
unnecessary memory.

Each process will print its data in the output*.txt (* = process rank).
The data includes the routing table, the root also prints the complete topology.