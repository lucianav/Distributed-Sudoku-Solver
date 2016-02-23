//Viziru Luciana
#include <stdio.h>

//method to read and save problem matrix
int* readProblem(char* dataFilename, int* dim_probl);

//method to read and save node topology
void readTopology(char* topologyFilename, int rank, int size,
									int* topology, int* numberOfNeighbours);

//method for probe-echo algorithm
void runProbeEcho(int size, int rank, int* parent, int* topology,
							MPI_Status* status, int numberOfNeighbours,
									int* numberOfChildren, int* routingTable);

//method to diffuse information in tree starting from root
void diffuseInfo(int rank, int size, int* parent, int* information,
						int infosize, MPI_Status* status, int* routingTable);

//method to allocate a square to each node in tree
void allocSubmatrixes(int rank, int size, int* parent, int allocated,
	int* mysquare, MPI_Status* status, int* numberOfChildren, int* routingTable);

// method to check if a matrix is valiid
int checkSolution(int* solution, int matrixdim);

//method to get valid results after completing node's square
int** getInitialResults(int rank, int dim_probl, int* solutionCount,
										int* problemMatrix, int mysquare);

//method for counting how many valid initial results there are
void countInitialSolutions(int* checkMatrix, int dim_probl, int* freeValues,
	int* freePositions, int count, int* solutionNumber, int start, int mysquare);

//methid to save the valid initial results
void saveInitialSolutions(int** initialResults, int dim_probl, int* freeValues,
							int *freePositions,	int count, int* solutionNumber,
									int solutionCount, int start, int mysquare);

//method to count valid solutions with info received from child
void countSolutions(int** initialResults, int solutionCount, int* count,
							int matrixdim, int* tempMatrix, int* checkMatrix);

//method to save valid solutions with info received from child
void mergeSolutions(int** initialResults, int solutionCount, int matrixdim,
	int* tempMatrix, int** newResults, int* solutionNumber, int newSolutionCount);
