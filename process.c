//Viziru Luciana
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "process.h"

#define probe 1
#define echo 2
#define diffuse 3
#define countSol 4
#define sendSol 5

int main(int argc, char* argv[]) {

	//MPI needed variables
	MPI_Status status;
	MPI_Request request;
	int rank, size;
	int parent; // parent node for probe-echo stage

	int dim_probl, i, j, k, numberOfNeighbours;

	// files needed by program
	char *topologyFilename = (char *) calloc(32, sizeof(char));
	strcpy(topologyFilename, argv[1]);
	char *dataFilename = (char *) calloc(32, sizeof(char));
	strcpy(dataFilename, argv[2]);
	char *outputFilename = (char *) calloc(32, sizeof(char));
	strcpy(outputFilename, argv[3]);

	// init
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	//alloc arrays for info
	int* topology = calloc(size * size, sizeof(int));
	int* numberOfChildren = calloc(size, sizeof(int));
	int* routingTable = calloc(size, sizeof(int));

	//read neigbours from file
	readTopology(topologyFilename, rank, size, topology, &numberOfNeighbours);

	//processing probe-echo stage
	runProbeEcho(size, rank, &parent, topology, &status,
							numberOfNeighbours, numberOfChildren, routingTable);

	//now the root has the complete topology
	//diffuse topology to all nodes
	diffuseInfo(rank, size, &parent, topology, size*size, &status, routingTable);

	//containers for info
	int* problemMatrix, *tempMatrix, *checkMatrix;

	//allocate and read matrix in root
	if (rank == 0) {
		problemMatrix = (int*)readProblem(dataFilename, &dim_probl);
	}
	//diffuse problem dimension to all nodes
	diffuseInfo(rank, size, &parent, &dim_probl, 1, &status, routingTable);

	//problem matrix dimension
	int matrixdim = dim_probl * dim_probl;

	//allocate matrix in other nodes
	if (rank != 0) {
		problemMatrix = calloc(matrixdim * matrixdim, sizeof(int));
	}

	//diffuse problem matrix to all nodes
	diffuseInfo(rank, size, &parent, problemMatrix, matrixdim * matrixdim,
														&status, routingTable);

	//allocate one square to each node
	int allocated = 0, mysquare, solutionCount;
	allocSubmatrixes(rank, size, &parent, allocated, &mysquare,
									&status, numberOfChildren, routingTable);

	//get all possible options for mysquare 
	int** initialResults = getInitialResults(rank, dim_probl, &solutionCount,
													problemMatrix, mysquare);

	//buffer for MPI_Recv
	tempMatrix = malloc(matrixdim * matrixdim * sizeof(int));
	//matrix needed for solution checking
	checkMatrix = malloc(matrixdim * matrixdim * sizeof(int));

	//get partial solutions from children
	int newSolutionCount = 0, nr;
	int** newResults = NULL;
	//for every node
	for (i = 0; i < size; i++) {
		//if it is child
		if (numberOfChildren[i] > 0) {

			//receive number of solutions
			MPI_Recv(&nr, 1, MPI_INT, i, countSol, MPI_COMM_WORLD, &status);

			//receive solutions and merge
			for (j = 0; j < nr; j++) {

				//receive one solution from child
				MPI_Recv(tempMatrix, matrixdim * matrixdim, MPI_INT, i, sendSol,
														MPI_COMM_WORLD, &status);
				int count = 0;
				//count how many combined valid solutions there are 
				countSolutions(initialResults, solutionCount, &count, matrixdim,
				tempMatrix, checkMatrix);

				//if there are any
				if (count > 0) {

					//allocate memmory for them in newResults
					newResults = realloc(newResults,
									(newSolutionCount + count) * sizeof(int*));
					for (k = 0; k < count; k++) {
						newResults[k + newSolutionCount] =
									malloc(matrixdim * matrixdim * sizeof(int));
					}
					
					int solutionNumber = newSolutionCount;
					newSolutionCount += count;

					//save the new solutions to newResults
					mergeSolutions(initialResults, solutionCount, matrixdim,
					tempMatrix, newResults, &solutionNumber, newSolutionCount);
				}
			}

			//if there were valid solutions with merge info
			//and there should be
			if (newSolutionCount > 0) {

				//free initialResults
				for (j = 0; j < solutionCount; j++) {
					free(initialResults[j]);
				}
				free(initialResults);

				//save the new results in initialResults
				initialResults = newResults;
				solutionCount = newSolutionCount;

				newSolutionCount = 0;
				newResults = NULL;
			}
		}
	}

	//if node is different from root
	if (rank != parent) {

		//send number of solutions to parent
		MPI_Send(&solutionCount, 1, MPI_INT, parent, countSol, MPI_COMM_WORLD);

		//send solutions to parent
		for (k = 0; k < solutionCount; k++) {
			MPI_Send(initialResults[k], matrixdim * matrixdim, MPI_INT, parent,
														sendSol, MPI_COMM_WORLD);
		}
	}
	else {
		//the root writes output to file
		FILE *outputFile = fopen(outputFilename, "w");
		fprintf(outputFile, "%d\n", dim_probl);
		for (i = 0; i < matrixdim * matrixdim; i++) {
			fprintf(outputFile, "%d ", initialResults[0][i]);
			if ((i + 1) % matrixdim == 0) {
				fprintf(outputFile, "\n");
			}
		}
		fclose(outputFile);
	}

	//print in my file my routing table
	char name[32];
	sprintf(name, "output%d.txt", rank);
	FILE *myout = fopen(name, "w");

	fprintf(myout, "%d\n", rank);
	for (i = 0; i < size; i++) {
		fprintf(myout, "destinatie %d next hop %d\n", i, routingTable[i]);
	}

	//process 0 prints topology also
	if (rank == 0) {
		for (i = 0; i < size; i++) {
			for (j = 0; j < size; j++) {
				fprintf(myout, "%d\t", topology[i * size + j]);
			}
			fprintf(myout, "\n");
		}
	}

	fclose(myout);

	//free memory allocated dynamically
	free(topologyFilename);
	free(dataFilename);
	free(outputFilename);

	free(topology);
	free(routingTable);
	free(numberOfChildren);
	free(problemMatrix);
	if (initialResults != NULL) {
		for (i = 0; i < solutionCount; i++) {
			free(initialResults[i]);
		}
		free(initialResults);
	}

	free(tempMatrix);
	free(checkMatrix);


	// finalize
	MPI_Finalize();

	return 0;
}

int* readProblem(char* dataFilename, int* dim_probl) {

	FILE *dataFile = fopen(dataFilename, "r");	//open file

	fscanf(dataFile, "%d", dim_probl);

	//alloc matrix
	int i, matrixdim = (*dim_probl) * (*dim_probl);
	int* problemMatrix = calloc(matrixdim * matrixdim, sizeof(int));

	//read data from file
	for (i = 0; i < matrixdim * matrixdim; i++) {
		fscanf(dataFile, "%d", (problemMatrix + i));
	}

	fclose(dataFile);	//close file
	return problemMatrix;
}

void readTopology(char* topologyFilename, int rank, int size,
									int* topology, int* numberOfNeighbours)
{
	FILE *topologyFile = fopen(topologyFilename, "r");	//open file

	char buff[256];
	int nr;
	(*numberOfNeighbours) = 0;
	//read line from file
	while (fgets(buff, sizeof(buff), topologyFile)) {
		sscanf(buff, "%d", &nr);
		//if the line has my information
    	if (nr == rank) {
    		//save all to topology
    		int neighbour;
    		char delim[5] = " \n\t-";
    		char *token = strtok(buff, delim);
    		token = strtok(NULL, delim);
    		while(token != NULL) {
    			sscanf(token, "%d", &neighbour);
    			topology[rank * size + neighbour] = 1;
    			(*numberOfNeighbours) ++;
    			token = strtok(NULL, delim);
    		}
    		//stop reading from file
    		break;
    	}
    }
    fclose(topologyFile);	//close file
}

void runProbeEcho(int size, int rank, int* parent, int* topology,
							MPI_Status* status, int numberOfNeighbours,
									int* numberOfChildren, int* routingTable)
{

	int* emptyTopology = calloc(size * size, sizeof(int));
	int* auxTopology = calloc(size * size, sizeof(int));
	int i, j;

	//start probe-echo stage from 0 node
	if (rank == 0) {
		//send to myself
		MPI_Send(emptyTopology, size * size, MPI_INT, 0, probe, MPI_COMM_WORLD);
	}

	//receive first probe
	MPI_Recv(auxTopology, size * size, MPI_INT, MPI_ANY_SOURCE,
										MPI_ANY_TAG, MPI_COMM_WORLD, status);
	if ((*status).MPI_TAG == probe) 	//this sould be happening
		*parent = (*status).MPI_SOURCE;

	//set all routes to parent at the beginning
	for (i = 0; i < size; i++) {
		routingTable[i] = *parent;
	}
	//next hop to self invalid
	routingTable[rank] = -1;

	//send probes to all neighbours except parent
	for (i = 0; i < size; i++) {
		if (topology[rank * size + i] == 1 && i != (*parent)) {
			MPI_Send(emptyTopology, size * size, MPI_INT, i,
													probe, MPI_COMM_WORLD);
		}
	}

	//number of echoes to receive
	int nrEchoes = numberOfNeighbours;
	//all connections except the first probe parent to intermediary nodes
	if (rank != 0) {
		nrEchoes --;
	}
	//wait for all echoes and send dummy echoes to new probes
	while (nrEchoes > 0) {
		MPI_Recv(auxTopology, size * size, MPI_INT, MPI_ANY_SOURCE,
										MPI_ANY_TAG, MPI_COMM_WORLD, status);
		int msgSource = (*status).MPI_SOURCE;
		//dummy echo
		if ((*status).MPI_TAG == probe) {
			//send echo
			MPI_Send(emptyTopology, size * size, MPI_INT, msgSource,
														echo, MPI_COMM_WORLD);

			//break link in topology
			topology[rank * size + msgSource] = 0;
		}
		//save new info
		if ((*status).MPI_TAG == echo) {
			int sum, checkEmpty = 1, countChildren = 0;
			//every node in graph
			for (i = 0; i < size; i++) {
				sum = 0;
				for (j = 0; j < size; j++) {
					if (auxTopology[i * size + j] == 1) {
						topology[i * size + j] = 1;
						checkEmpty = 0;
						sum ++;
					}
				}
				//info about node found
				if (sum > 0) {
					//count as child in subtree
					countChildren ++;
					//can be accessed thorugh the source of the message
					routingTable[i] = msgSource;
				}
			}
			if (checkEmpty == 1) {
				//break link in topology
				topology[rank * size + msgSource] = 0;
			}
			//set number of nodes in subtree starting form msgSource
			numberOfChildren[msgSource] = countChildren;
			nrEchoes --;
		}
	}

	//send info to parent after all echoes are received
	if (rank != 0) {
		MPI_Send(topology, size * size, MPI_INT, *parent, echo, MPI_COMM_WORLD);
	}

	//wait for all to send partial topologies
	MPI_Barrier(MPI_COMM_WORLD);

	free(emptyTopology);
	free(auxTopology);
}

void diffuseInfo(int rank, int size, int* parent, int* information,
						int infosize, MPI_Status* status, int* routingTable)
{
	int i;
	if (rank == 0) {
		// start diffusion from root
		for (i = 0; i < size; i++) {
			if (routingTable[i] == i) //direct child of root
				//send information
				MPI_Send(information, infosize, MPI_INT, i, 
												diffuse, MPI_COMM_WORLD);
		}
	}
	else {
		//first receive from parent
		MPI_Recv(information, infosize, MPI_INT, *parent, diffuse,
													MPI_COMM_WORLD, status);
		//send to children
		for (i = 0; i < size; i++) {
			if (routingTable[i] == i && i != *parent) //direct child
				MPI_Send(information, infosize, MPI_INT, i,
													diffuse, MPI_COMM_WORLD);
		}
	}

	//wait for all to receive information
	MPI_Barrier(MPI_COMM_WORLD);
}

void allocSubmatrixes(int rank, int size, int* parent, int allocated,
	int* mysquare, MPI_Status* status, int* numberOfChildren, int* routingTable)
{
	int i;
	//if i am root, just send
	if (rank == 0) {
		(*mysquare) = allocated;	//choose square
		allocated ++;	//mark as allocated
		for (i = 0; i < size; i++) {
			if (routingTable[i] == i) {//direct child of root
				//send number of already taken squares
				MPI_Send(&allocated, 1, MPI_INT, i, diffuse, MPI_COMM_WORLD);
				//add number of squares in branch
				allocated += numberOfChildren[i];	
			}
		}
	}
	else {
		//first receive from parent
		MPI_Recv(&allocated, 1, MPI_INT, *parent, diffuse, MPI_COMM_WORLD, status);
		(*mysquare) = allocated;	//choose square
		allocated ++;	//count as taken
		//send to children
		for (i = 0; i < size; i++) {
			if (routingTable[i] == i && i != *parent) {//direct child
				//send to child  current number of taken square
				MPI_Send(&allocated, 1, MPI_INT, i, diffuse, MPI_COMM_WORLD);
				//count squares in branch
				allocated += numberOfChildren[i];
			}
		}
	}
}

int checkSolution(int* solution, int matrixdim) {
	int i, j, k;
	//check row and column conditions

	//for every position
	for (i = 0; i < matrixdim; i++) {
		for (j = 0; j < matrixdim; j++) {
			//check columns that follow
			for (k = j + 1; k < matrixdim; k++) {
				if (solution[i * matrixdim + j] != 0 &&
					solution[i * matrixdim + j] == solution[i * matrixdim + k]){
					return 0;
				}
			}
			//check lines that follow
			for (k = i + 1; k < matrixdim; k++) {
				if (solution[i * matrixdim + j] != 0 &&
					solution[i * matrixdim + j] == solution[k * matrixdim + j]) {
					return 0;
				}
			}
		}
	}
	//no problem found, solution is valid
	return 1;
}

int** getInitialResults(int rank, int dim_probl, int* solutionCount,
										int* problemMatrix, int mysquare)
{
	int i, j, l, c, k, matrixdim = dim_probl * dim_probl;
	int possibleValues[matrixdim + 1];
	//suppose all values are possible
	for (i = 0; i <= matrixdim; i++) {
		possibleValues[i] = 1;
	}

	//count empty cells in count variable
	int count = 0;
	l = (mysquare / dim_probl) * dim_probl;
	c = (mysquare % dim_probl) * dim_probl;
	//count number of values needed to fill square
	for (i = l; i < l + dim_probl; i++) {
		for (j = c; j < c + dim_probl; j++) {
			if (problemMatrix[i * matrixdim + j] == 0) {
				count ++;
			}
			else {
				//if value exists in square, is no longer valid
				possibleValues[problemMatrix[i * matrixdim + j]] = 0;

			}
		}
	}

	//containers for free positions in suqares and possible values for them
	int freePositions[count], freeValues[count];
	k = 0;
	for (i = 1; i <= matrixdim; i++) {
		//value i is missing from square
		if (possibleValues[i] == 1) {
			//save it to container
			freeValues[k ++] = i;
		}
	}
	k = 0;
	for (i = l; i < l + dim_probl; i++) {
		for (j = c; j < c + dim_probl; j++) {
			//position i * matrixdim + j is missing from square
			if (problemMatrix[i * matrixdim + j] == 0) {
				//save it to container
				freePositions[k ++] = (i - l) * dim_probl + (j - c);
			}
		}
	}

	//generate all solutions by permutations of the possible values
	int solutionNumber = 0;

	//count number of possible initial solutions in solutionNumber 
	countInitialSolutions(problemMatrix, dim_probl, freeValues,
	freePositions, count, &solutionNumber, 0, mysquare);

	(*solutionCount) = solutionNumber;

	//allocate memory for the itinial results matrix
	//and set them to problem matrix
	int** initialResults = malloc(solutionNumber * sizeof(int*));
	for (k = 0; k < solutionNumber; k++) {
		initialResults[k] = malloc(matrixdim * matrixdim * sizeof(int));
		for (i = 0; i < matrixdim * matrixdim; i++) {
			initialResults[k][i] = problemMatrix[i];
		}
	}

	//start filling initialResult form position 0
	solutionNumber = 0;

	//save solutions to initialResults
	saveInitialSolutions(initialResults, dim_probl, freeValues,
	  		freePositions, count, &solutionNumber, *solutionCount, 0, mysquare);

	//return 3d matrix of all correct possibilities
	return initialResults;
}

void swap(int* a, int* b) {
	int aux = *a;
	*a = *b;
	*b = aux;
}

void countInitialSolutions(int* checkMatrix, int dim_probl, int* freeValues,
	int* freePositions, int count, int* solutionNumber, int start, int mysquare)
{
	int i, j, k, matrixdim = dim_probl * dim_probl;
	int l = (mysquare / dim_probl) * dim_probl;
	int c = (mysquare % dim_probl) * dim_probl;

	//all postions were matched
	if (start >= count) {
		k = 0;
		//write solution to checkMatrix
		for (i = 0; i < dim_probl; i++) {
			for (j = 0; j < dim_probl; j++) {
				if (i * dim_probl + j == freePositions[k]) {
					checkMatrix[(l + i) * matrixdim + (c + j)] = freeValues[k++];
				}
			}
		}
		//check if solution is valid
		if (checkSolution(checkMatrix, matrixdim) == 1) {
			//if so, count it
			(*solutionNumber) ++;
		}
	}
	else {
		//set value i to start position
		for (i = start; i < count; i++) {

			swap((freeValues + start), (freeValues + i));

			//generate all following permutations
			countInitialSolutions(checkMatrix, dim_probl, freeValues,
				freePositions, count, solutionNumber, start + 1, mysquare);

			//swap back to initial state
			swap((freeValues + start), (freeValues + i));
		}
	}
}

void saveInitialSolutions(int** initialResults, int dim_probl, int* freeValues,
							int *freePositions,	int count, int* solutionNumber,
									int solutionCount, int start, int mysquare)
{
	//if all solutions were already found, return
	if (*solutionNumber >= solutionCount) {
		return;
	}

	int i, j, k, matrixdim = dim_probl * dim_probl;
	int l = (mysquare / dim_probl) * dim_probl;
	int c = (mysquare % dim_probl) * dim_probl;

	//all postions were matched
	if (start >= count) {
		k = 0;
		//write current matrix to initialResults
		for (i = 0; i < dim_probl; i++) {
			for (j = 0; j < dim_probl; j++) {
				if (i * dim_probl + j == freePositions[k]) {
					initialResults[*solutionNumber]
							[(l + i) * matrixdim + (c + j)] = freeValues[k++];
				}
			}
		}
		//if current matrix is a valid solution, count al leave it
		//otherwise, it will be overwritten
		if (checkSolution(initialResults[*solutionNumber], matrixdim) == 1) {
			(*solutionNumber) ++;
		}
	}
	else {
		//set value i to start position
		for (i = start; i < count; i++) {

			swap((freeValues + start), (freeValues + i));

			//generate all following permutations
 			saveInitialSolutions(initialResults, dim_probl, freeValues,
 					freePositions, count, solutionNumber, solutionCount,
 														start + 1, mysquare);

			//swap back to initial state
			swap((freeValues + start), (freeValues + i));
		}
	}
}

void countSolutions(int** initialResults, int solutionCount, int* count,
							int matrixdim, int* tempMatrix, int* checkMatrix)
{

	int i, j, k;
	//for every solution in initialResults
	for (k = 0; k < solutionCount; k++) {
		//combine initialResults[k] and matrix from child in checkMatrix
		for (i = 0; i < matrixdim; i++) {
			for (j = 0; j < matrixdim; j++) {
				checkMatrix[i * matrixdim + j] =
										initialResults[k][i * matrixdim + j];
				if (tempMatrix[i * matrixdim + j] != 0 &&
								initialResults[k][i * matrixdim + j] == 0) {
					checkMatrix[i * matrixdim + j] =
												tempMatrix[i * matrixdim + j];
				}
			}
		}

		//if solution is valid
		if (checkSolution(checkMatrix, matrixdim) == 1) {
			//count it
			(*count) ++;
		}
	}
}

void mergeSolutions(int** initialResults, int solutionCount, int matrixdim,
	int* tempMatrix, int** newResults, int* solutionNumber, int newSolutionCount)
{

	int i, j, k;
	//for every solution in initialResults
	for (k = 0; k < solutionCount; k++) {
		//combine initialResults[k] and matrix from child	
		for (i = 0; i < matrixdim; i++) {
			for (j = 0; j < matrixdim; j++) {
				//set to initialResults matrix
				newResults[*solutionNumber][i * matrixdim + j] =
										initialResults[k][i * matrixdim + j];
				//tempMatrix has new info for [i][j], then write
				if (tempMatrix[i * matrixdim + j] != 0 &&
								initialResults[k][i * matrixdim + j] == 0) {
					newResults[*solutionNumber][i * matrixdim + j] =
												tempMatrix[i * matrixdim + j];
				}
			}
		}

		//if solution is valid
		if (checkSolution(newResults[*solutionNumber], matrixdim) == 1) {
			//count it to save it in newResults
			(*solutionNumber) ++;
			//if all solutions were found, stop searching
			if (*solutionNumber >= newSolutionCount) {
				break;
			}
		}
	}
}
