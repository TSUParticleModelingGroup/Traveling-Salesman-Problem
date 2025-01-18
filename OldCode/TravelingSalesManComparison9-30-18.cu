//Maybe add a random pertibation
//I don't like that the wall strength is just an arbitrary large number
//work on pointer in exact
//Move the point off (0,0) before running the extrution program
//Move wall by a set number because the wall will not be a fixed distance apart
//
// Set optimal radius off of average seperation
// Clean up extrution
// If min distance is 1 find max distance then how big would a grid be that fit this. Use this to find optimal radius
//
//nvcc TravelingSalesManComparison9-30-18.cu -o TSPCompare093018 -lglut -lGL -lm
#include <GL/glut.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265359
#define BIG_NUMBER 1000000.0

#define BLOCK 256

#define X_WINDOW 800
#define Y_WINDOW 800

#define X_MAX (1.1)
#define X_MIN (-1.1)

#define Y_MAX (1.1)
#define Y_MIN (-1.1)

#define RANDOM_FILL_SHAPE_CIRCLE -1
#define RANDOM_FILL_SHAPE_SQUARE 1

#define TIME_STEP_SIZE 0.002
#define STEPS_BETWEEN_VIEWING 100

#define NUMBER_OF_WALL_MOVES 1000
#define RADIUS_STEP_SIZE 0.001
#define TIME_BETWEEN_WALL_MOVES 1.0

#define RELAX_TIME 10.0

#define WALL_STRENGTH 20000.0

#define SLOPE_REPULSION -100.0
#define	SLOPE_ATRACTION -1.0

#define DAMP 20.0 
#define MASS 80.0

#define FLOAT_ROUND_OFF 0.000001

#define DRAW_EXHAUSTIVE_PATH 1
#define DRAW_NEAREST_NEIGHBOR_PATH -1
#define DRAW_NBODY_EXTRUSION_PATH 1
#define DELAY_TO_RECORD -1
#define PRINT_EDGE_COST -1
#define PRINT_EXHAUSTIVE_PATHS -1
#define PRINT_PATHS -1
#define PRINT_RAW_DATA_FILE 1

// Globals
FILE *RawDataFile;
FILE *StatsFile;

//Function prototypes
float x_machine_to_x_screen(int x);
float y_machine_to_y_screen(int y);
float x_machine_to_x_world(int x);
float y_machine_to_y_world(int y);
float x_world_to_x_screen(float x);
float y_world_to_y_screen(float y);
void openRawDataFile(int scope, int numberOfRuns);
void placeNodesRandom(float4 *pos, unsigned int srandSeed, int scope, int n);
int checkNodes(float4 *pos, int n);
void placeNodesGrid(float4 *pos, int rows, int columns);
void getNumberOfNodesFromNodeFile(int *numberOfNodes);
float adjustingNodes(float4 *pos, int numberOfNodes);
void placeNodesFromAFile(float4 *pos, int *numberOfNodes);
void printEdgeCosts(float4 *pos, int n, float nodeAdjustmentFactor);
int factorial(int n);
void printPathOrder(int* path, int n);
float getPathCost(int *path, float4 *pos, int type, int n);
void swap(int *path, int i, int j);
void heappermute(int* path, int m, float4 *pos, int* exhaustivePath, float* minCost, int n);
float exhaustiveTSP(float4 *pos, int* exhaustivePath, int n);
float nearestNeighborTSP(float4 *pos, int* path, int n);
void setNbodyInitailConditions(float4 *nodePos, float4 *pos, float4 *vel, float* mass, int n);
float4 adjustSoGeometricCenterIsZero(float4 *pos, int n);
float4 FindGeometericCenter(float4 *pos, int n);
float findAverageSeperation(float4 *pos, int n);
float findOptimalOuterRadius(float4 *pos, int n);
void drawPoints(float4 pos, int n);
void drawNbodyExtrusion(float4 *pos, float4 center, float outerRadius, float innerRadius, float optimalOuterRadius, int inOut, int n);
void getPathNbody(float4 *pos, float4 geometericCenter, int* path, int n);
float NbodyExtrusionTSP(float4 *nodePos, float4 *pos, float4 *vel, float4 *acc, float* mass, int* path, int n);
void drawFInalPicture(float4 *pos, int *pathA, int *pathB, int *pathC, int scope, int n);
void getInputFromUser(int* scope, int* numberOfNodes, int* numberOfRuns, int* maxNumberOfRows, int* maxNumberOfColumns, unsigned int* srandSeed);
void control();

float x_machine_to_x_screen(int x)
{
	return( (2.0*x)/X_WINDOW-1.0 );
}

float y_machine_to_y_screen(int y)
{
	return( -(2.0*y)/Y_WINDOW+1.0 );
}

/*	Takes machine x and y which start in the upper left corner and go from zero to X_WINDOW
	left to right and form zero to Y_WINDOW top to bottom and transslates this into world
	points which are a X_MIN to X_MAX, Y_MIN to Y_MAX window.
*/
float x_machine_to_x_world(int x)
{
	float range;
	range = X_MAX - X_MIN;
	return( (range/X_WINDOW)*x + X_MIN );
}

float y_machine_to_y_world(int y)
{
	float range;
	range = Y_MAX - Y_MIN;
	return(-((range/Y_WINDOW)*y - Y_MAX));
}

/*	Take world  points to screen points
*/
float x_world_to_x_screen(float x)
{
	float range;
	range = X_MAX - X_MIN;
	return( -1.0 + 2.0*(x - X_MIN)/range );
}

float y_world_to_y_screen(float y)
{
	float range;
	range = Y_MAX - Y_MIN;
	return( -1.0 + 2.0*(y - Y_MIN)/range );
}

void openRawDataFile(int scope, int numberOfRuns)
{
	char tagName[50];
	char fileName[256];
	
	strcpy(tagName,"");
	if(scope == 2)
	{
		strcat(tagName,"Random_TSP_Raw_Data_CSV");
	}
	if(scope == 4)
	{
		strcat(tagName,"Grid_TSP_Raw_Data_CSV");
	}
	
	snprintf(fileName, 256, "%s", tagName);
	
	RawDataFile = fopen(fileName, "wb");
  	
	fprintf(RawDataFile, "  Number of runs %d\n",numberOfRuns);
	
	if(scope == 2)
	{
		fprintf(RawDataFile, "  Run  Nodes  Exact  NNeighbor  NBody\n\n");
	}
	if(scope == 4)
	{
		fprintf(RawDataFile, "  Run  Rows  Columns  Exact  NNeighbor  NBody\n\n");
	}
}

void placeNodesRandom(float4 *pos, unsigned int srandSeed, int scope, int n)
{
	time_t t;
	float xStart, xStop, yStart, yStop, radius, maxRadius;
	float mag, seperation;
	int repeatedNodeTest;
	
	xStart = -1.0;
	xStop  = 1.0;
	yStart = -1.0;
	yStop  = 1.0;
	
	maxRadius = 1.0;

	srand((unsigned) time(&t));

	if(scope == 1) srand(srandSeed);
	else srand((unsigned) time(&t));

	if(RANDOM_FILL_SHAPE_SQUARE == 1)
	{
		for(int i = 0; i < n; i++)
		{
			repeatedNodeTest = 0;
			while(repeatedNodeTest == 0)
			{
				pos[i].x = (float)rand()/((float)RAND_MAX)*(xStop - xStart) + xStart;
				pos[i].y = (float)rand()/((float)RAND_MAX)*(yStop - yStart) + yStart;
				repeatedNodeTest = 1;
				for(int j = 0; j < i; j++)
				{
					seperation = sqrt((pos[i].x-pos[j].x)*(pos[i].x-pos[j].x) + (pos[i].y-pos[j].y)*(pos[i].y-pos[j].y));
					if(seperation < FLOAT_ROUND_OFF)
					{
						repeatedNodeTest = 0;
						break;
					}
				}
			}
		}
	}
	
	if(RANDOM_FILL_SHAPE_CIRCLE == 1)
	{
		for(int i = 0; i < n; i++)
		{
			repeatedNodeTest = 0;
			while(repeatedNodeTest == 0)
			{
				pos[i].x = (float)rand()/((float)RAND_MAX)*(xStop - xStart) + xStart;
				pos[i].y = (float)rand()/((float)RAND_MAX)*(yStop - yStart) + yStart;
				mag = sqrt(pos[i].x*pos[i].x + pos[i].y*pos[i].y);
				radius = ((float)rand()/(float)RAND_MAX)*maxRadius;
				pos[i].x *= radius/mag;
				pos[i].y *= radius/mag;
				repeatedNodeTest = 1;
				for(int j = 0; j < i; j++)
				{
					seperation = sqrt((pos[i].x-pos[j].x)*(pos[i].x-pos[j].x) + (pos[i].y-pos[j].y)*(pos[i].y-pos[j].y));
					if(seperation < FLOAT_ROUND_OFF)
					{
						repeatedNodeTest = 0;
						break;
					}
				}
			}
		}
	}
}

int checkNodes(float4 *pos, int n)
{
	float seperation;
	
	for(int i = 0; i < (n - 1); i++)
	{
		for(int j = (i + 1); j < n; j++)
		{
			seperation = sqrt((pos[i].x-pos[j].x)*(pos[i].x-pos[j].x) + (pos[i].y-pos[j].y)*(pos[i].y-pos[j].y));
			if(seperation < FLOAT_ROUND_OFF)
			{
				return(-1);
			}
		}
	}
	return(1);
}

void placeNodesGrid(float4 *pos, int rows, int columns)
{
	int i,j,k;
	float dx, dy;
	float xStart, yStart;
	
	xStart = 0.0;
	yStart = 0.0;
	dx = 1.0;
	dy = 1.0;
	
	k = 0;
	for(i = 0; i < columns; i++)
	{
		for(j = 0; j < rows; j++)
		{
			pos[k].x = xStart + dx*i;
			pos[k].y = yStart + dy*j;
			k++;
		}
	}
}

void getNumberOfNodesFromNodeFile(int *numberOfNodes)
{
	FILE *nodeFile;
	nodeFile = fopen("nodeFile", "rb");
	fscanf(nodeFile,"%d", numberOfNodes);
	fclose(nodeFile);
}

void placeNodesFromAFile(float4 *pos, int *numberOfNodes)
{
	FILE *nodeFile;
	int nodeNumber;
	nodeFile = fopen("nodeFile", "rb");
	fscanf(nodeFile,"%d", numberOfNodes);
	
	for(int i = 0; i < *numberOfNodes; i++)
	{
		fscanf(nodeFile,"%d %f %f", &nodeNumber, &pos[i].x, &pos[i].y );
		//printf("\n\nnode number= %d, %f, %f\n\n", nodeNumber, pos[i].x, pos[i].y);
	}
	
	fclose(nodeFile);
}

//This function adjustes the nodes so that the minimal seperation is 1.0
float adjustingNodes(float4 *pos, int numberOfNodes)
{
	float temp;
	float nodeAdjustmentFactor;
	float minSeperation;
	float dx,dy;
	int i,j;
	
	dx = pos[1].x-pos[2].x;
	dy = pos[1].y-pos[2].y;
	minSeperation = sqrt(dx*dx + dy*dy);
	for(i = 0; i < numberOfNodes; i++)
	{
		for(j = i + 1; j < numberOfNodes; j++)
		{
			dx = pos[i].x-pos[j].x;
			dy = pos[i].y-pos[j].y;
			temp = sqrt(dx*dx + dy*dy);
			if(temp < minSeperation)
			{
				minSeperation = temp;
			}
		}
	}
	
	nodeAdjustmentFactor = 1.0/minSeperation;
	
	for(int i = 0; i < numberOfNodes; i++)
	{
		pos[i].x = pos[i].x*nodeAdjustmentFactor;
		pos[i].y = pos[i].y*nodeAdjustmentFactor;
	}
	
	return(nodeAdjustmentFactor);
}

float findDistanceToOuterMostNode(float4 *pos, int numberOfNodes)
{
	float temp;
	float distanceToOutermostNode = 0.0;
	
	for(int i = 0; i < numberOfNodes; i++)
	{
		temp = sqrt(pos[i].x*pos[i].x + pos[i].y*pos[i].y);
		if(temp > distanceToOutermostNode)
		{
			distanceToOutermostNode = temp;
		}
	}
	return(distanceToOutermostNode);
}

void printEdgeCosts(float4 *pos, int n, float nodeAdjustmentFactor)
{
	float temp;
	for(int i = 0; i < n; i++)
	{
		for(int j = i + 1; j < n; j++)
		{	
			temp = sqrt((pos[i].x-pos[j].x)*(pos[i].x-pos[j].x) + (pos[i].y-pos[j].y)*(pos[i].y-pos[j].y))/nodeAdjustmentFactor;
			printf("edge cost [%d, %d] = %f\n", i, j, temp);
		}
	}
}

int factorial(int n)
{
	int outPut = n;
	
	for(int i = n-1; i > 0; i--)
	{
		outPut *= i;	
	}
	return(outPut);
}

void printPathOrder(int* path, int n)
{
	printf("  ");
	for(int i = 0; i < n-1; i++)
	{
		printf("%d->", path[i]);	
	}
	printf("%d", path[n-1]);
}

float getPathCost(int *path, float4 *pos, int type, int n)
{
	float cost;
	int i, j, k;
	
	//Checking path validaty 
	for(i = 0; i < n; i++)
	{
		if(path[i] < 0 || (n-1) < path[i])
		{
			printf("\n\n  Error -> Path out of range! Type = %d", type);
			printf("\n  path[%d] = %d\n\n", i, path[i]);
			printf("\n\n  Good Bye.  \n\n");
			exit(0);
		}
		
		for(j = 0; j < i; j++)
		{
			if(path[i] == path[j])
			{
				printf("\n\n Error -> Path has a repeated index! Type = %d\n", type);
				printPathOrder(path, n);
				printf("\n\n");
				printf("\n\n  Good Bye.  \n\n");
				exit(0);
			}
		}
	}
	
	cost = 0.0;
	for(k = 0; k < n-1; k++)
	{
		i = path[k];
		j = path[k+1];
		cost += sqrt((pos[i].x-pos[j].x)*(pos[i].x-pos[j].x) + (pos[i].y-pos[j].y)*(pos[i].y-pos[j].y));
	}
	i = path[n-1];
	j = path[0];
	cost += sqrt((pos[i].x-pos[j].x)*(pos[i].x-pos[j].x) + (pos[i].y-pos[j].y)*(pos[i].y-pos[j].y));
	
	return(cost);
}

void swap(int *path, int i, int j)
{
	int temp;
	temp = path[i];
	path[i] = path[j];
	path[j] = temp;
}

void heappermute(int *path, int m, float4 *pos, int *exhaustivePath, float *minCost, int n) 
{
	int i;
	float pathCost;
	int* pathPlus = (int*)malloc(n*sizeof(int));

	if (m == 1) 
	{
		pathPlus[0] = 0;
		for(i = 1; i < n; i++)
		{
			pathPlus[i] = path[i-1];	
		}
		
		pathCost = getPathCost(pathPlus, pos, 1, n);
		
		if(PRINT_EXHAUSTIVE_PATHS == 1)
		{
			printf("\n");
			printPathOrder(pathPlus, n);
			printf(" cost = %f", pathCost);
		}
		
		if(pathCost < minCost[0])
		{
			minCost[0] = pathCost;
			for(i = 0; i < n; i++)
			{
				exhaustivePath[i] = pathPlus[i];	
			}
		}
    	}
	else 
	{
		for (i = 0; i < m; i++) 
		{
			heappermute(path, m-1, pos, exhaustivePath, minCost, n);
			if (m % 2 == 1) 
			{
				swap(path, 0, m-1);
			}
			else 
			{
				swap(path, i, m-1);
			}
		}
	}
	free(pathPlus);
}

float exhaustiveTSP(float4 *pos, int* exhaustivePath, int n)
{
	float cost[1];
	int* path = (int*)malloc((n-1)*sizeof(int));
	
	exhaustivePath[0] = 0;
	for(int i = 1; i < n; i++)
	{
		exhaustivePath[i] = i;
		path[i-1] = i;	
	}
	cost[0] = getPathCost(exhaustivePath, pos, 1, n);
	
	heappermute(path, n-1, pos, exhaustivePath, cost, n);
	free(path);
	return(cost[0]);
}

float nearestNeighborTSP(float4 *pos, int* path, int n)
{
	int i, j, k, nextNode, nodeFound;
	float minCost, pathCost, edgeCost, maxEdgeCost;
	int* used = (int*)malloc(n*sizeof(int));
	
	maxEdgeCost = 0.0;
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			edgeCost = sqrt((pos[i].x-pos[j].x)*(pos[i].x-pos[j].x) + (pos[i].y-pos[j].y)*(pos[i].y-pos[j].y));
			if(edgeCost > maxEdgeCost) 
			{
				maxEdgeCost = edgeCost;
			}	
		}	
	}
	maxEdgeCost += 1.0;
	
	for(i = 0; i < n; i++)
	{
		used[i] = -1;	
	}
	
	path[0] = 0;
	used[0] = 1;
	
	k = 0;
	
	minCost = maxEdgeCost;
	while(k < n-1)
	{
		nodeFound = 0;
		for(j = 0; j < n; j++)
		{
			i = path[k];
			edgeCost = sqrt((pos[i].x-pos[j].x)*(pos[i].x-pos[j].x) + (pos[i].y-pos[j].y)*(pos[i].y-pos[j].y));
			if(edgeCost <= minCost && used[j] == -1)
			{
				minCost = edgeCost;
				nextNode = j;
				nodeFound = 1;
			}	
		}
		if(nodeFound == 0)
		{
			printf("\n\n  There was a problem in the nearest neighbor function. No next node was found.\n\n");
			printf("\n\n  Good Bye.  \n\n");
			exit(0);
		}
		nodeFound = 0;
		
		k++;
		path[k] = nextNode;
		used[nextNode] = 1;
		minCost = maxEdgeCost;
	}
	
	pathCost = getPathCost(path, pos, 2, n);
	free(used);
	return(pathCost);
}

void setNbodyInitailConditions(float4 *nodePos, float4 *pos, float4 *vel, float* mass, int n)
{
	int i;

	for(i = 0; i < n; i++)
	{
		pos[i].x = nodePos[i].x;
		pos[i].y = nodePos[i].y;
		
		vel[i].x = 0.0;
		vel[i].y = 0.0;
		
		mass[i] = MASS;
	}
}

void moveAnyNodeOffDeadCenter(float4 *pos, int n)
{
	int i;

	for(i = 0; i < n; i++)
	{
		if( sqrt(pos[i].x*pos[i].x + pos[i].y*pos[i].y) < 0.001) 
		{
			pos[i].x = 0.001;
			pos[i].y = 0.001;
		}
	}
}

float4 adjustSoGeometricCenterIsZero(float4 *pos, int n)
{
	float4 geometricCenter;
	
	geometricCenter.x = 0.0;
	geometricCenter.y = 0.0;
	geometricCenter.z = 0.0;
	geometricCenter.w = 0.0;
	
	for(int i = 0; i < n; i++)
	{
		geometricCenter.x += pos[i].x;
		geometricCenter.y += pos[i].y;
	}
	
	geometricCenter.x /= (float)n;
	geometricCenter.y /= (float)n;
	
	for(int i = 0; i < n; i++)
	{
		pos[i].x -= geometricCenter.x;
		pos[i].y -= geometricCenter.y;
	}
	return(geometricCenter);
}

float4 FindGeometericCenter(float4 *pos, int n)
{
	float4 massCenter;
	
	massCenter.x = 0.0;
	massCenter.y = 0.0;
	massCenter.z = 0.0;
	massCenter.w = 0.0;
	
	for(int i = 0; i < n; i++)
	{
		massCenter.x += pos[i].x*MASS;
		massCenter.y += pos[i].y*MASS;
	}
	
	massCenter.x /= (float)n*MASS;
	massCenter.y /= (float)n*MASS;
	
	return(massCenter);
}

float findOptimalOuterRadius(float4 *pos, int n)
{
	float maxRadius;
	float perimeter;
	int *nearestNeighborPath;
	
	//The min node seperation is 1. So this lets all node fit on the circle with min seperation.
	//maxRadius = (n*1.0)/(2.0*PI);
	
	// Setting outer radius to be a fraction of the nearest nieghbor length
	nearestNeighborPath = (int*)malloc((n)*sizeof(int));
	perimeter = nearestNeighborTSP(pos, nearestNeighborPath, n);
	maxRadius = 0.9*perimeter/(2.0*PI);
	free(nearestNeighborPath);
	
	
	return(maxRadius);
}

float findAverageSeperation(float4 *pos, int n)
{
	float averageSeperation;
	float totalSeperation = 0.0;
	
	for(int i = 0; i < n; i++)
	{
		for(int j = i + 1; j < n; j++)
		{
			totalSeperation += sqrt((pos[i].x-pos[j].x)*(pos[i].x-pos[j].x) + (pos[i].y-pos[j].y)*(pos[i].y-pos[j].y));
		}
	}
	averageSeperation = ((float)(n*n - n))/2.0;
	
	return(averageSeperation);
}

void drawPoints(float4 *pos, int n)
{
	int i;
	float outerRadius = findDistanceToOuterMostNode(pos, n);
	float normalizingFactor = outerRadius; //((float)n)/IDEAL_NUMBER_OF_NODES;
	
	glClear(GL_COLOR_BUFFER_BIT);
	
	glPointSize(5.0);
	glColor3f(1.0,0.0,0.0);
	for(i = 0; i < n; i++)
	{
		glBegin(GL_POINTS);
		glVertex2f(x_world_to_x_screen(pos[i].x/normalizingFactor),y_world_to_y_screen(pos[i].y/normalizingFactor));
		glEnd();

	}
	glFlush();
}

void drawNbodyExtrusion(float4 *pos, float4 center, float outerRadius, float innerRadius, float optimalOuterRadius, int inOut, int n)
{
	int i;
	int lineAmount = 100;
	float normalizingFactor = outerRadius;
	
	optimalOuterRadius/= normalizingFactor;
	outerRadius /= normalizingFactor;
	innerRadius /= normalizingFactor;

	glClear(GL_COLOR_BUFFER_BIT);
	
	GLfloat twicePi = 2.0f * PI;
	
	glLineWidth(1.0);
	glColor3f(1.0,1.0,1.0);
	glBegin(GL_LINE_LOOP);
		for(i = 0; i <= lineAmount;i++) 
		{ 
			glVertex2f(x_world_to_x_screen(center.x + (optimalOuterRadius*cos(i*twicePi/lineAmount))), 
			           y_world_to_y_screen(center.y + (optimalOuterRadius*sin(i*twicePi/lineAmount))));
		}
	glEnd();
	
	glLineWidth(1.0);
	if(inOut == -1) glColor3f(1.0,0.0,0.0);
	else glColor3f(0.0,1.0,0.0);
	glBegin(GL_LINE_LOOP);
		for(i = 0; i <= lineAmount;i++) 
		{ 
			glVertex2f(x_world_to_x_screen(center.x + (outerRadius*cos(i*twicePi/lineAmount))), 
			           y_world_to_y_screen(center.y + (outerRadius*sin(i*twicePi/lineAmount))));
		}
	glEnd();
	
	glLineWidth(1.0);
	glColor3f(1.0,1.0,0.0);
	glBegin(GL_LINE_LOOP);
		for(i = 0; i <= lineAmount;i++) 
		{ 
			glVertex2f(x_world_to_x_screen(center.x + (innerRadius*cos(i*twicePi/lineAmount))), 
			           y_world_to_y_screen(center.y + (innerRadius*sin(i*twicePi/lineAmount))));
		}
	glEnd();
	
	glPointSize(5.0);
	glColor3f(1.0,0.0,0.0);
	for(i = 0; i < n; i++)
	{
		glBegin(GL_POINTS);
		glVertex2f(x_world_to_x_screen(pos[i].x/normalizingFactor),y_world_to_y_screen(pos[i].y/normalizingFactor));
		glEnd();

	}
	
	glFlush();
}

__global__ void accelerationsNbody(float4 *nodePos, float4 *pos, float4 *vel, float4 *acc, float *mass, float4 geometericCenter, float innerRadius, float outerRadius, int n)
{
	int j,ii;
    //float3 bodyBodyForce, wallForce, forceSum;
    float3 forceSum;
    float4 nodePosMe, posMe;
    float dx, dy, d, edgeLength; 
    //float dc, temp;
    float radius, forceMag;
    __shared__ float4 shNodePos[BLOCK], shPos[BLOCK];
    int id = threadIdx.x + blockDim.x*blockIdx.x;
    
    forceSum.x = 0.0;
	forceSum.y = 0.0;
	
	nodePosMe.x = nodePos[id].x;
	nodePosMe.y = nodePos[id].y;
	posMe.x = pos[id].x;
	posMe.y = pos[id].y;
		    
    for(j=0; j < gridDim.x; j++)
    {
    	if(threadIdx.x + blockDim.x*j < n)
    	{
    		shPos[threadIdx.x] = pos[threadIdx.x + blockDim.x*j];
    		shNodePos[threadIdx.x] = nodePos[threadIdx.x + blockDim.x*j];
    	}
    	__syncthreads();
   
		#pragma unroll 32
        for(int i = 0; i < blockDim.x; i++)	
        {
        	ii = i + blockDim.x*j;
		    if(ii != id && ii < n) 
		    {
				dx = shPos[i].x - posMe.x;
				dy = shPos[i].y - posMe.y;
				d = sqrtf(dx*dx + dy*dy);
				//dc = (sqrtf(posMe.x*posMe.x + posMe.y*posMe.y) + sqrtf(shPos[i].x*shPos[i].x + shPos[i].y*shPos[i].y))/2.0;
				//temp = (outerRadius - dc)/outerRadius;
				//dc = abs(shPos[i].x*posMe.y - shPos[i].y*posMe.x)/d;
				edgeLength = sqrtf((shNodePos[i].x - nodePosMe.x)*(shNodePos[i].x - nodePosMe.x) + (shNodePos[i].y - nodePosMe.y)*(shNodePos[i].y - nodePosMe.y));
				
				//if(dc < 0.5*innerRadius) forceMag = 0.0;
				if(d <= edgeLength)
				{
					forceMag = (edgeLength - d)*SLOPE_REPULSION;
	
				}
				else if(d <= 1.2*edgeLength)
				{
					//forceMag = 25.0*temp*(edgeLength - d)*SLOPE_ATRACTION;
					forceMag = (edgeLength - d)*SLOPE_ATRACTION;
					//forceMag = 0.0;
				}
				else
				{
					forceMag = 0.0;
				}
				
				forceSum.x += forceMag*dx/d;
				forceSum.y += forceMag*dy/d;
		    }
		}
	}
	
	if(id < n)
	{
		// Forces between node and the walls
		dx = posMe.x - geometericCenter.x;
		dy = posMe.y - geometericCenter.y; 
		radius = sqrtf(dx*dx + dy*dy);
	
		if(radius < innerRadius) // Inside inner wall
		{
			forceMag = WALL_STRENGTH*(innerRadius - radius);
			forceSum.x += forceMag*dx/radius;
			forceSum.y += forceMag*dy/radius;
		}
		else if(radius > outerRadius) // Outside outer wall
		{
			forceMag = WALL_STRENGTH*(outerRadius - radius);
			forceSum.x += forceMag*dx/radius;
			forceSum.y += forceMag*dy/radius;
		}
		
		// Adding on damping force.
		forceSum.x += -DAMP*vel[id].x;
		forceSum.y += -DAMP*vel[id].y;
		
		// Creating the accelerations.
	    acc[id].x = forceSum.x/mass[id];
	    acc[id].y = forceSum.y/mass[id];
    }
}

__global__ void moveNbody(float4 *pos, float4 *vel, float4 *acc, float dt, int n)
{
    int id = threadIdx.x + blockDim.x*blockIdx.x;
    if(id < n)
    {
	    vel[id].x += acc[id].x*dt;
		vel[id].y += acc[id].y*dt;
		
		pos[id].x  += vel[id].x*dt;
		pos[id].y  += vel[id].y*dt;
    }
}

void getPathNbody(float4 *pos, float4 geometericCenter, int* path, int n)
{
	int i;
	float minValue;
	float *angle = (float*)malloc(n*sizeof(float));
	int *used = (int*)malloc(n*sizeof(int));
	
	for(i = 0; i < n; i++)
	{
		if(pos[i].x == 0 && pos[i].y == 0)
		{
			angle[i] = 0.0;
		}
		else if(pos[i].x >= 0 && pos[i].y >= 0)
		{
			if(pos[i].x == 0) angle[i] = 90.0;
			else angle[i] = atan(pos[i].y/pos[i].x)*180.0/PI;
		}
		else if(pos[i].x < 0 && pos[i].y >= 0)
		{
			angle[i] = 180.0 - atan(pos[i].y/(-pos[i].x))*180.0/PI;
		}
		else if(pos[i].x <= 0 && pos[i].y < 0)
		{
			if(pos[i].x == 0) angle[i] = 270.0;
			else angle[i] = 180.0 + atan(pos[i].y/pos[i].x)*180.0/PI;
		}
		else
		{
			angle[i] = 360.0 - atan(-pos[i].y/pos[i].x)*180.0/PI;
		}
	}
	
	for(i = 0; i < n; i++)
	{
		used[i] = 0;
	}
	
	for(int k = 0; k < n; k++)
	{
		minValue = 400.0;
		for(i = 0; i < n; i++)
		{
			if(angle[i] < minValue && used[i] == 0)
			{
				minValue = angle[i];
				path[k] = i;
			}
		}
		used[path[k]] = 1;
		//printf("path[%d] = %d\n", k, path[k]);
	}
	
	free(angle);
	free(used);
}

float NbodyExtrusionTSP(float4 *nodePos, float4 *pos, float4 *vel, float4 *acc, float* mass, int* path, int n)
{
	int draw_count;
	int outerWallDirection;
	float dr;
	float dt = TIME_STEP_SIZE;
	float pathCost;
	float time;
	float4 geometericCenter;
	float innerRadius, outerRadius, optimalOuterRadius;
	float averageSeperation;
	
	dim3 block, grid;
	float4 *posGPU, *velGPU, *accGPU; 
	float *massGPU;
	float4 *nodePosGPU;
	int wallMovesLeft;
	float nodeDensity, optimalDensity;
	
	// Setting up GPU parrellel structure.
	block.x = BLOCK;
	block.y = 1;
	block.z = 1;
	
	grid.x = (n-1)/block.x + 1;
	grid.y = 1;
	grid.z = 1;
	
	// Allocating memory.
	cudaMalloc( (void**)&nodePosGPU, n *sizeof(float4));
	cudaMalloc( (void**)&posGPU, n *sizeof(float4));
	cudaMalloc( (void**)&velGPU, n *sizeof(float4));
	cudaMalloc( (void**)&accGPU, n *sizeof(float4));
	cudaMalloc( (void**)&massGPU, n *sizeof(float));
	
	// This is used to pause the program so you can setup to take a video of a run.
	if(DELAY_TO_RECORD == 1)
	{
		printf("\n\n  Enter a character to start\n\n"); getchar();
	}
	
	// Copying information up to the GPU.
	cudaMemcpy( nodePosGPU, nodePos, n *sizeof(float4), cudaMemcpyHostToDevice );
	cudaMemcpy( posGPU, pos, n *sizeof(float4), cudaMemcpyHostToDevice );
    cudaMemcpy( velGPU, vel, n *sizeof(float4), cudaMemcpyHostToDevice );
    cudaMemcpy( massGPU, mass, n *sizeof(float), cudaMemcpyHostToDevice );
    
    // Moving the nodes into a circle.
    averageSeperation = findAverageSeperation(nodePos, n);
    optimalOuterRadius = findOptimalOuterRadius(nodePos, n);
	dr = (optimalOuterRadius - 0.0)/(float)NUMBER_OF_WALL_MOVES;
    outerWallDirection = -1;
    innerRadius = 0.0;
    outerRadius = findDistanceToOuterMostNode(pos, n);
    geometericCenter = FindGeometericCenter(pos, n);
	drawNbodyExtrusion(pos, geometericCenter, outerRadius, innerRadius, optimalOuterRadius, outerWallDirection, n);
	//optimalDensity = 2.23; //based off of 4X4 grid
	//nodeDensity = ((float)n)/(PI*outerRadius*outerRadius);
	time = 0.0;
	draw_count = 0;
	while(optimalOuterRadius/2.0 < outerRadius)
	{
		outerRadius -= dr;
		while(time < TIME_BETWEEN_WALL_MOVES)
		{		
			accelerationsNbody<<<grid, block>>>(nodePosGPU, posGPU, velGPU, accGPU, massGPU, geometericCenter, innerRadius, outerRadius, n);
			moveNbody<<<grid, block>>>(posGPU, velGPU, accGPU, dt, n);
			
			if(draw_count == STEPS_BETWEEN_VIEWING)
			{
				cudaMemcpy( pos, posGPU, n *sizeof(float4), cudaMemcpyDeviceToHost );
				drawNbodyExtrusion(pos, geometericCenter, outerRadius, innerRadius, optimalOuterRadius, outerWallDirection, n);
				draw_count = 0;
			}
			draw_count++;
			time += dt;
		}
		
		//nodeDensity = ((float)n)/(PI*outerRadius*outerRadius);
		time = 0.0;
		outerRadius -= dr;
	}
	
	wallMovesLeft = NUMBER_OF_WALL_MOVES;
	outerWallDirection = -1;
	while(innerRadius + 0.01 < optimalOuterRadius)
	{
		while(time < TIME_BETWEEN_WALL_MOVES)
		{		
			accelerationsNbody<<<grid, block>>>(nodePosGPU, posGPU, velGPU, accGPU, massGPU, geometericCenter, innerRadius, outerRadius, n);
			moveNbody<<<grid, block>>>(posGPU, velGPU, accGPU, dt, n);
			
			if(draw_count == STEPS_BETWEEN_VIEWING)
			{
				cudaMemcpy( pos, posGPU, n *sizeof(float4), cudaMemcpyDeviceToHost );
				drawNbodyExtrusion(pos, geometericCenter, outerRadius, innerRadius, optimalOuterRadius, outerWallDirection, n);
				draw_count = 0;
			}
			draw_count++;
			time += dt;
		}
		
		// Reseting run conditions.
		time = 0.0;
		innerRadius += dr;
		//wallMovesLeft--;
		//printf("\n Number of moves left = %d", wallMovesLeft);
		if((outerRadius - innerRadius) < 1.0) //averageSeperation/2.0);
		{ 
			outerWallDirection = 1.0;
		}
		if(outerRadius < optimalOuterRadius) outerRadius += dr*outerWallDirection;
	}
	
	getPathNbody(pos, geometericCenter, path, n);
	pathCost = getPathCost(path, nodePos, 3, n);
	
	return(pathCost);
}

void drawFInalPicture(float4 *pos, int *pathA, int *pathB, int *pathC, int scope, int n)
{	
	int i;
	float outerRadius = findDistanceToOuterMostNode(pos, n);
	float normalizingFactor = outerRadius; //((float)n)/IDEAL_NUMBER_OF_NODES;

	glClear(GL_COLOR_BUFFER_BIT);
	
	//exhuastivePath path
	if(scope == 1 || scope == 2)
	{
		if(DRAW_EXHAUSTIVE_PATH == 1)
		{
			glLineWidth(4.0);
			glColor3f(0.0,0.0,1.0);
			glBegin(GL_LINE_LOOP);
				for(i = 0; i < n; i++)
				{
					glVertex2f(x_world_to_x_screen(pos[pathA[i]].x/normalizingFactor),y_world_to_y_screen(pos[pathA[i]].y/normalizingFactor));
				}
			glEnd();
		}
	}
	
	//Nearest Neighbor path
	if(DRAW_NEAREST_NEIGHBOR_PATH == 1)
	{
		glLineWidth(6.0);
		glColor3f(0.0,1.0,0.0);
		glBegin(GL_LINE_LOOP);
			for(i = 0; i < n; i++)
			{
				glVertex2f(x_world_to_x_screen(pos[pathB[i]].x/normalizingFactor),y_world_to_y_screen(pos[pathB[i]].y/normalizingFactor));
			}
		glEnd();
	}
	
	//Nbody Extrusion path
	if(DRAW_NBODY_EXTRUSION_PATH == 1)
	{
		glLineWidth(4.0);
		glColor3f(1.0,0.0,0.0);
		glBegin(GL_LINE_LOOP);
			for(i = 0; i < n; i++)
			{
				glVertex2f(x_world_to_x_screen(pos[pathC[i]].x/normalizingFactor),y_world_to_y_screen(pos[pathC[i]].y/normalizingFactor));
			}
		glEnd();
	}
	
	//Placing nodes
	glPointSize(8.0);
	glColor3f(1.0,1.0,1.0);
	for(i = 0; i < n; i++)
	{
		glBegin(GL_POINTS);
			glVertex2f(x_world_to_x_screen(pos[i].x/normalizingFactor),y_world_to_y_screen(pos[i].y/normalizingFactor));
		glEnd();
	}
	
	//Nearest neighbor start node 
	if(DRAW_NEAREST_NEIGHBOR_PATH == 1)
	{
		glPointSize(10.0);
		glColor3f(0.0,0.0,1.0);
		glBegin(GL_POINTS);
			glVertex2f(x_world_to_x_screen(pos[pathB[0]].x/normalizingFactor),y_world_to_y_screen(pos[pathB[0]].y/normalizingFactor));
		glEnd();
	}
	
	//Nbody extrution start and stop nodes
	if(DRAW_NBODY_EXTRUSION_PATH == 1)
	{
		glPointSize(10.0);
		glColor3f(0.0,1.0,0.0);
		glBegin(GL_POINTS);
			glVertex2f(x_world_to_x_screen(pos[pathC[0]].x/normalizingFactor),y_world_to_y_screen(pos[pathC[0]].y/normalizingFactor));
		glEnd();
	
		glColor3f(1.0,0.0,0.0);
		glBegin(GL_POINTS);
			glVertex2f(x_world_to_x_screen(pos[pathC[n-1]].x/normalizingFactor),y_world_to_y_screen(pos[pathC[n-1]].y/normalizingFactor));
		glEnd();
	}
	
	glFlush();
}

void getInputFromUser(int* scope, int* numberOfNodes, int* numberOfRuns, int* maxNumberOfRows, int* maxNumberOfColumns, unsigned int* srandSeed)
{
	*scope = -1;
	*numberOfNodes = -1;
	*numberOfRuns = -1;
	*maxNumberOfRows = -1;
	*maxNumberOfColumns = -1;
	
	printf("\n\n  What type run would you like to perform?");
	printf("\n  1 for one small randomly generated run.");
	printf("\n  2 for a series of small randomly generated run.");
	printf("\n  3 for one on on a grid.");
	printf("\n  4 for a series of runs on randomly generated sized grids.");
	printf("\n  5 for a large randomly generated run.");
	printf("\n  6 to read nodes from nodeFile.");
	printf("\n\n  Inter an integer value: ");
	scanf("%d", scope);
	
	if(*scope == 1)
	{
		printf("\n\n  You will need to enter the number of nodes (Cities) wound you like to generate?");
		printf("\n  Note: If you choose a number bigger than 13 you may lock your computer up.");
		printf("\n        This is because we will be testing against the exaustive algorithim which has ((n-1)!)/2 paths.");
		
		printf("\n\n  You will also need to enter a seed to generate the random placed nodes");
		printf("\n  By intering the same seed on different exicutions of the program");
		printf("\n  you will be able to run multiple tests on the same node configuration.");
		
		printf("\n  Enter two positive integer values seperated by a space: \n  (number of nodes) (seed)");
		printf("\n\n  Inter your values: ");
		scanf("%d %d", numberOfNodes, srandSeed);
		
		*numberOfRuns = 1;
	}
	else if(*scope == 2)
	{
		printf("\n\n  You will need to enter the number nodes (Cities) wound like to generate?");
		printf("\n  Note: If you choose a nuber bigger than 13 you may lock your computer up.");
		
		printf("\n\n  You will also need to enter how many randomly generated runs would you like to perform?");
		
		printf("\n  Enter two positive integer values seperated by a space: \n  (number of nodes) (number of runs)");
		printf("\n\n  Inter your values: ");
		scanf("%d %d", numberOfNodes, numberOfRuns);
	}
	else if(*scope == 3)
	{
		printf("\n  You will need to enter the number of rows and columns for the grid you will generate.");
	
		printf("\n  Enter two positive integer values seperated by a space: \n  (number of rows) (number of columns)");
		printf("\n\n  Inter your values: ");
		scanf("%d %d", maxNumberOfRows, maxNumberOfColumns);
		
		*numberOfRuns = 1;
	}
	else if(*scope == 4)
	{
		printf("\n  You will need to enter the number of rows and columns for the grid you will generate");
		printf("\n  and the number of runs");
	
		printf("\n  Enter three positive integer values seperated by a space: \n  (max number of rows) (max number of columns) (number of runs)");
		printf("\n\n  Inter your values: ");
		scanf("%d %d %d", maxNumberOfRows, maxNumberOfColumns, numberOfRuns);
	}
	else if(*scope == 5)
	{
		printf("\n\n  You will need to enter the number of nodes (Cities) wound you like to generate?");
		printf("\n  Enter a positive integer values: \n  (number of nodes)");
		printf("\n\n  Inter your values: ");
		scanf("%d", numberOfNodes);
		
		*numberOfRuns = 1;
	}
	else if(*scope == 6)
	{	
		*numberOfRuns = 1;
	}
	else
	{
		printf("\n\n  Your input for type of run was invalide.");
		printf("\n\n  Good Bye.  \n\n");
		exit(0);
	}
}

void control()
{
	time_t t;
	int scope, numberOfNodes, numberOfRuns, maxNumberOfRows, maxNumberOfColumns;
	unsigned int srandSeed;
	int rows, columns, done;
	float4 *pos;
	float nodeAdjustmentFactor;
	float4 geometricCenter;
	float distanceToOutermostNode;
	int *exhaustivePath, *nearestNeighborPath, *NbodyExtrusionPath;
	float4 *posNbody, *velNbody, *accNbody; 
	float *massNbody;
	float exhaustiveCost, nearestNeighborCost, NbodyExtrusionCost;
	int nodeCheck;
	float temp;
	
	getInputFromUser(&scope, &numberOfNodes, &numberOfRuns, &maxNumberOfRows, &maxNumberOfColumns, &srandSeed);
	
	if(scope == 2 || scope == 4 && PRINT_RAW_DATA_FILE == 1)
	{
		openRawDataFile(scope, numberOfRuns);
	}
	
	float totalNearestNeighborCost = 0.0;
	float totalNbodyExtrusionCost = 0.0;
	float totalPercentErrorNearestNeighbor = 0.0;
	float totalPercentErrorNbodyExtrusion = 0.0;
	float NbodyExtrusionVSNearestNeighbor = 0.0;
	
	for(int i = 0; i < numberOfRuns; i++)
	{	
		printf("\n\n\n  ********************* Intermediate Run %d ********************* ", i+1);
		
		exhaustiveCost = BIG_NUMBER;
		nearestNeighborCost = BIG_NUMBER;
		NbodyExtrusionCost = BIG_NUMBER;
		
		//Alocating memory
		if(scope == 3)
		{
			rows = maxNumberOfRows;
			columns = maxNumberOfColumns;
			numberOfNodes = rows*columns;
		}	
		if(scope == 4)
		{
			srand((unsigned) time(&t));
			done = -1;
			while(done == -1)
			{
				rows = 1 + (float)rand()/((float)RAND_MAX)*maxNumberOfRows;
				columns = 1 + (float)rand()/((float)RAND_MAX)*maxNumberOfColumns;
				printf("\n\n  rows = %d columns = %d", rows, columns);
				numberOfNodes = rows*columns;
				if(rows == 1 && columns == 1) done = -1;
				else done = 1;
			}
		}
		if(scope == 6)
		{
			getNumberOfNodesFromNodeFile(&numberOfNodes);
		}
		
		pos = (float4*)malloc((numberOfNodes)*sizeof(float4));
	
		exhaustivePath = (int*)malloc((numberOfNodes)*sizeof(int));
		nearestNeighborPath = (int*)malloc((numberOfNodes)*sizeof(int));
		NbodyExtrusionPath = (int*)malloc((numberOfNodes)*sizeof(int));

		posNbody = (float4*)malloc((numberOfNodes)*sizeof(float4));
		velNbody = (float4*)malloc((numberOfNodes)*sizeof(float4));
		accNbody = (float4*)malloc((numberOfNodes)*sizeof(float4));
		massNbody = (float*)malloc((numberOfNodes)*sizeof(float4));
		
		//Creating nodes
		if(scope == 1 || scope == 2 || scope == 5)
		{	
			placeNodesRandom(pos, srandSeed, scope, numberOfNodes);
		}
		else if(scope == 3 || scope == 4)
		{
			placeNodesGrid(pos, rows, columns);
		}
		else if(scope == 6)
		{
			placeNodesFromAFile(pos, &numberOfNodes);
		}
		
		//Adjusting nodes
		geometricCenter = adjustSoGeometricCenterIsZero(pos, numberOfNodes);
		printf("\n\n  The geometric center of the nodes = (%f, %f)", geometricCenter.x, geometricCenter.y);
		
		distanceToOutermostNode = findDistanceToOuterMostNode(pos, numberOfNodes);
		printf("\n  The distance to the outermost node from the geometric center pre adjustment is %f", distanceToOutermostNode);
		
		nodeAdjustmentFactor = adjustingNodes(pos, numberOfNodes);
		printf("\n  The node adjustment factor = %f", nodeAdjustmentFactor);
		
		distanceToOutermostNode = findDistanceToOuterMostNode(pos, numberOfNodes);
		printf("\n  The distance to the outermost node from the geometric center post adjustment is %f", distanceToOutermostNode);
		
		//Checking to see if a node is repeated
		nodeCheck = checkNodes(pos, numberOfNodes);
		if(nodeCheck == -1)
		{
			printf("\n\n  There is a repeated node. Check your data set.");
			printf("\n\n  Good Bye.  \n\n");
			exit(0);
		}
		
		//Drawing the adjusted nodes on the screen.
		drawPoints(pos, numberOfNodes); 
		
		//Printing the edge costs (lengths in this case)
		if(PRINT_EDGE_COST == 1)
		{
			printEdgeCosts(pos, numberOfNodes, nodeAdjustmentFactor);
		}
		
		//Finding exact cost
		printf("\n\n  Determining the exact cost.");
		if(scope == 1 || scope == 2)
		{	
			exhaustiveCost = exhaustiveTSP(pos, exhaustivePath, numberOfNodes);
		}
		else if(scope == 3 || scope == 4)
		{
			//Assuming all edges are the same length. So just get the length of the first edge.
			temp = sqrt((pos[0].x-pos[1].x)*(pos[0].x-pos[1].x) + (pos[0].y-pos[1].y)*(pos[0].y-pos[1].y));
			if(rows == 1 || columns == 1)
			{
				exhaustiveCost = temp*2.0*(numberOfNodes - 1);
			}
			else if(rows%2 == 0 || columns%2 == 0)
			{
				exhaustiveCost = temp*numberOfNodes;
			}
			else
			{
				exhaustiveCost = temp*(numberOfNodes - 1.0 + sqrt(2.0));
			}
		}
		else if(scope == 5)
		{
			exhaustiveCost = -1.0;
		}
		else if(scope == 6)
		{
			if(numberOfNodes < 14)
			{
				exhaustiveCost = exhaustiveTSP(pos, exhaustivePath, numberOfNodes);
			}
			else
			{
				exhaustiveCost = -1.0;
			}
		}
		printf("\n  Determining the exact cost is done.");
		
		//Finding nearest neighbor cost
		printf("\n\n  Running the nearest nieghbor algorithm.");
		nearestNeighborCost = nearestNeighborTSP(pos, nearestNeighborPath, numberOfNodes);
		printf("\n  The nearest nieghbor algorithm is done.");
		
		//Running n-body extrusion code
		printf("\n\n  Running the N-body extrusion algorithm."); 
		printf("  \n"); //I had to enter this carage return so it would print the line above before it started the algorithm
		setNbodyInitailConditions(pos, posNbody, velNbody, massNbody, numberOfNodes);
		moveAnyNodeOffDeadCenter(posNbody, numberOfNodes);
		NbodyExtrusionCost = NbodyExtrusionTSP(pos, posNbody, velNbody, accNbody, massNbody, NbodyExtrusionPath, numberOfNodes);
		printf("  The N-body extrusion algorithm is done.");
		
		//Unadjusting costs
		exhaustiveCost /= nodeAdjustmentFactor;
		nearestNeighborCost /= nodeAdjustmentFactor;
		NbodyExtrusionCost /= nodeAdjustmentFactor;
		
		totalNearestNeighborCost += nearestNeighborCost;
		totalNbodyExtrusionCost += NbodyExtrusionCost;
		
		//Sanity check
		if(nearestNeighborCost < exhaustiveCost - FLOAT_ROUND_OFF)
		{
			printf("\n\n  Nearest Neighbor cost (%f) is smaller than exhaustive cost (%f). Something is wrong!\n",nearestNeighborCost, exhaustiveCost);
			printf("\n\n  Good Bye.  \n\n");
			exit(0);
		}
		if(NbodyExtrusionCost < exhaustiveCost - FLOAT_ROUND_OFF)
		{
			printf("\n\n  Nbody Extrution cost (%f) is smaller than exhaustive cost (%f). Something is wrong!\n",NbodyExtrusionCost, exhaustiveCost);
			printf("\n\n  Good Bye.  \n\n");
			exit(0);
		}
		
		printf("\n\n  --------------------- Intermediate Run Results --------------------- ");
		
		// This is for debugging
		if(PRINT_PATHS == 1)
		{
			if(scope == 1 || scope == 2)
			{
				printf("\n\n  The exhaustive  path is     : "); 
				printPathOrder(exhaustivePath, numberOfNodes); 
				printf(" cost = %f", exhaustiveCost);
			}
			
			printf("\n\n  The nearest neighbor path is: "); 
			printPathOrder(nearestNeighborPath, numberOfNodes); 
			printf(" cost = %f", nearestNeighborCost);
			
			printf("\n\n  The Nbody extrusion path is : "); 
			printPathOrder(NbodyExtrusionPath, numberOfNodes); 
			printf(" cost = %f", NbodyExtrusionCost);
		}
		
		// Printing out the single run stats and acumulating the multiple run info to create final stats.
		// Stephen your stat collection should go here.
		if(exhaustiveCost < 0.0)
		{
			printf("\n\n  The minimum cost is unknown");
			printf("\n\n  The nearest neighbor cost is : %f", nearestNeighborCost);
			printf("\n\n  The Nbody extrusion cost is  : %f", NbodyExtrusionCost);
			
			NbodyExtrusionVSNearestNeighbor += (nearestNeighborCost - NbodyExtrusionCost)/nearestNeighborCost;
		}
		else
		{
			printf("\n\n  The minimum cost is          : %f ", exhaustiveCost);
			printf("\n\n  The nearest neighbor cost is : %f the precent error = %f", nearestNeighborCost, 100.0*(nearestNeighborCost - exhaustiveCost)/exhaustiveCost);
			printf("\n\n  The Nbody extrusion cost is  : %f the precent error = %f", NbodyExtrusionCost, 100.0*(NbodyExtrusionCost - exhaustiveCost)/exhaustiveCost);
			
			totalPercentErrorNearestNeighbor += 100.0*(nearestNeighborCost - exhaustiveCost)/exhaustiveCost;
			totalPercentErrorNbodyExtrusion  += 100.0*(NbodyExtrusionCost  - exhaustiveCost)/exhaustiveCost;
			
			NbodyExtrusionVSNearestNeighbor += (nearestNeighborCost - NbodyExtrusionCost)/nearestNeighborCost;
		}
		
		if(scope == 2 && PRINT_RAW_DATA_FILE == 1)
		{
			fprintf(RawDataFile, "  %d, %d, %f, %f, %f\n", i+1, numberOfNodes, exhaustiveCost, nearestNeighborCost, NbodyExtrusionCost);
		}
		if(scope == 4 && PRINT_RAW_DATA_FILE == 1)
		{
			fprintf(RawDataFile, "  %d, %d, %d, %f, %f, %f\n", i+1, rows, columns, exhaustiveCost, nearestNeighborCost, NbodyExtrusionCost);
		}
	
		drawFInalPicture(pos, exhaustivePath, nearestNeighborPath, NbodyExtrusionPath, scope, numberOfNodes);
		
		free(pos);
		free(exhaustivePath);
		free(nearestNeighborPath);
		free(NbodyExtrusionPath);
		free(posNbody);
		free(velNbody);
		free(accNbody);
		free(massNbody);
	}
	
	printf("\n\n\n  $$$$$$$$$$$$$$$$$$$$$$$$$ Final results $$$$$$$$$$$$$$$$$$$$$$$$$$$$$");
	
	// Printing out the final acumulated stats.
	// Stephen your stat final stats should go here.
	if(exhaustiveCost < 0.0)
	{
		printf("\n\n  The average value of the nearest neighbor method was %f on %d run(s).", totalNearestNeighborCost/numberOfRuns, numberOfRuns);
		printf("\n\n  The average value of the Nbody extrution method was %f on %d run(s).", totalNbodyExtrusionCost/numberOfRuns, numberOfRuns);
	}
	else
	{
		printf("\n\n  The average percent error of the nearest neighbor method was %f on %d runs.", totalPercentErrorNearestNeighbor/(float)numberOfRuns, numberOfRuns);
		printf("\n  The average percent error of the Nbody extrution method was %f on %d runs", totalPercentErrorNbodyExtrusion/(float)numberOfRuns, numberOfRuns);
	}
	
	NbodyExtrusionVSNearestNeighbor = 100.0*NbodyExtrusionVSNearestNeighbor/(float)numberOfRuns;
	if(NbodyExtrusionVSNearestNeighbor >= 0)
	{
		printf("\n\n  The Nbody ectrusion method was on average %f percent better than the nearest neighbor method on %d run(s).", NbodyExtrusionVSNearestNeighbor, numberOfRuns);
	}
	else
	{
		printf("\n\n  The Nbody ectrusion method was on average %f percent worse than the nearest neighbor method on %d run(s).", -NbodyExtrusionVSNearestNeighbor, numberOfRuns);
	}
	
	if(scope == 2 || scope == 4 && PRINT_RAW_DATA_FILE == 1)
	{
		fclose(RawDataFile);
	}
	
	printf("\n\nDone\n");
	while(1);
}

int main(int argc, char** argv)
{
	glutInit(&argc,argv);
	glutInitWindowSize(X_WINDOW,Y_WINDOW);
	glutInitWindowPosition(0,0);
	glutCreateWindow("Traveling Salesman Problem");
	glutDisplayFunc(control);
	glutMainLoop();
}



    

