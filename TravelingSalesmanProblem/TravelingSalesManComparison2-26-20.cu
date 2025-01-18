//Maybe add a random pertibation
//I don't like that the wall strength is just an arbitrary large number
//work on pointer in exact
//Move the point off (0,0) before running the extrution program
//Move wall by a set number because the wall will not be a fixed distance apart
//
// Clean up extrution
// 
// Think about setting atraction as gravity and repultion as linear. WIll need epsilon to remove sengularity
//
// ** Here I am removing wall pressure and creating stress to move the wall (sum of stress divided by 0.5*n*(n-1)).
// In this version I am going to try and add anealing

//nvcc TravelingSalesManComparison2-26-20.cu -o temp -lglut -lGL -lm -use_fast_math
#include <GL/glut.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// General constants
#define PI 3.14159265359
#define BIG_NUMBER 1000000.0
#define SMALL_NUMBER 0.0000001

// GPU settings
#define BLOCK 256

// OpenGL drawing settings
#define X_WINDOW 800
#define Y_WINDOW 800
#define X_MAX (1.1)
#define X_MIN (-1.1)
#define Y_MAX (1.1)
#define Y_MIN (-1.1)

// Error cheaking settings
#define DRAW_EXHAUSTIVE_PATH 1
#define DRAW_NEAREST_NEIGHBOR_PATH -1
#define DRAW_NBODY_EXTRUSION_PATH 1
#define DELAY_TO_RECORD -1
#define PRINT_EDGE_COST -1
#define PRINT_EXHAUSTIVE_PATHS -1
#define PRINT_PATHS -1
#define PRINT_RAW_DATA_FILE 1

// Node creation settings
#define RANDOM_FILL_SHAPE_CIRCLE -1
#define RANDOM_FILL_SHAPE_SQUARE 1
#define NODES_TOO_CLOSE 0.000001

// General N-body settings
#define TIME_STEP_SIZE 0.002
#define DAMP 20.0 
#define MASS 80.0
#define STEPS_BETWEEN_VIEWING 100

// Force function settings
#define FORCE_FUNTION_TYPE 1 // 0 for LJ-HPQ. 1 for LJ-MPQ. 2 Linear. 3 for billards.

#define H 1.0
#define P 7.0 
#define Q 11.5  // P < Q 
#define M -0.04

#define SLOPE_REPULSION 100.0
#define	MAG_ATRACTION 10.0 //5.0
#define FORCE_CUTOFF 10.0 //2.0

#define NODE_RADIUS 0.5

// Wall move general
#define STARTING_POINT_FOR_NUMBER_OF_MOVES 1000
#define RADIUS_STEP_SIZE 0.001
#define TIME_BETWEEN_WALL_MOVES 1.0
#define WALL_STRENGTH 200000.0

// Wall move type and selection settings
#define WALL_MOVE_SCOPE 1 // 0 for nothing. 1 for pressure. 2 for average stress. 3 for average absolute stress. 4 for average positive stress. 5 for max stress

#define LOWER_PRESSURE_LIMIT 20
#define UPPER_PRESSURE_LIMIT 1000

#define LOWER_AVERAGE_STRESS_LIMIT 10
#define UPPER_AVERAGE_STRESS_LIMIT 30

#define LOWER_AVERAGE_ABSOLUTE_STRESS_LIMIT 10
#define UPPER_AVERAGE_ABSOLUTE_STRESS_LIMIT 30

#define LOWER_AVERAGE_POSITIVE_STRESS_LIMIT 10
#define UPPER_AVERAGE_POSITIVE_STRESS_LIMIT 30

#define LOWER_MAX_POSITIVE_STRESS_FROM_OUTER_WALL_LIMIT 30
#define UPPER_MAX_POSITIVE_STRESS_FROM_OUTER_WALL_LIMIT 3000

// Annealing settings
#define NUMBER_OF_ANEALINGS  0
#define ANNEALING_TIME  50.0

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
void placeNodesRandom(float4 *node, unsigned int srandSeed, int scope, int n);
int checkNodes(float4 *node, int n);
void placeNodesGrid(float4 *node, int rows, int columns);
double setAverageSeperationToOne(float4 *node, int numberOfNodes);
double setMinimumSeperationToOne(float4 *node, int numberOfNodes);
float4 setGeometricCenterToZero(float4 *node, int n);
void getNumberOfNodesFromNodeFile(int *numberOfNodes, double *exhaustiveCost, char *nodeFileName);
void placeNodesFromAFile(float4 *node, int *numberOfNodes, char *nodeFileName);
void getNumberOfNodesFromASelfDefinedFunction(int *numberOfNodes);
void placeNodesFromASelfDefinedFunction(float4 *node, int numberOfNodes);
void printEdgeCosts(float4 *node, int n, double nodeAdjustmentFactor);
int factorial(int n);
void printPathOrder(int* path, int n);
double getPathCost(int *path, float4 *node, int type, int n);
void swap(int *path, int i, int j);
void heappermute(int* path, int m, float4 *node, int* exhaustivePath, double* minCost, int n);
double exhaustiveTSP(float4 *node, int* exhaustivePath, int n);
double nearestNeighborTSP(float4 *node, int* path, int n);
void setNbodyInitailConditions(float4 *node, float4 *pos, float4 *vel, float* mass, int n);
void drawPoints(float4 *pos, int n);
void drawNbodyExtrusion(float4 *pos, float innerRadius, float outerRadius, int innerWallDirection, int outerWallDirection, int n);
void getPathNbody(float4 *pos, int* path, int n);
double findMinimumSeperation(float4 *pos, int n);
int findWallMoveDirections(float4 *node, float4 *pos, int n, int scope, float outerRadius, int *innerWallDirection, int *outerWallDirection);
double NbodyExtrusionTSP(float4 *node, float4 *pos, float4 *vel, float4 *acc, float* mass, int* path, int n);
void drawFInalPicture(float4 *node, int *pathA, int *pathB, int *pathC, int scope, int n);
void getInputFromUser(int* scope, int* numberOfNodes, int* numberOfRuns, int* maxNumberOfRows, int* maxNumberOfColumns, unsigned int* srandSeed, char *nodeFileName);
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
	if(scope == 2 || scope == 6)
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
	
	if(scope == 2 || scope == 6)
	{
		fprintf(RawDataFile, "  Run  Nodes  Exact  NNeighbor  NBody\n\n");
	}
	if(scope == 4)
	{
		fprintf(RawDataFile, "  Run  Rows  Columns  Exact  NNeighbor  NBody\n\n");
	}
}

void placeNodesRandom(float4 *node, unsigned int srandSeed, int scope, int n)
{
	time_t t;
	double xStart, xStop, yStart, yStop, radius, maxRadius;
	double mag, seperation;
	int repeatedNodeTest;
	
	xStart = -sqrt(n)/2.0;
	xStop  = sqrt(n)/2.0;
	yStart = -sqrt(n)/2.0;
	yStop  = sqrt(n)/2.0;
	
	maxRadius = sqrt(n/2.0);

	srand((unsigned) time(&t));

	if(scope == 1 || scope == 5) srand(srandSeed);
	else srand((unsigned) time(&t));

	if(RANDOM_FILL_SHAPE_SQUARE == 1)
	{
		for(int i = 0; i < n; i++)
		{
			repeatedNodeTest = 0;
			while(repeatedNodeTest == 0)
			{
				node[i].x = (double)rand()/((double)RAND_MAX)*(xStop - xStart) + xStart;
				node[i].y = (double)rand()/((double)RAND_MAX)*(yStop - yStart) + yStart;
				repeatedNodeTest = 1;
				for(int j = 0; j < i; j++)
				{
					seperation = sqrt((node[i].x-node[j].x)*(node[i].x-node[j].x) + (node[i].y-node[j].y)*(node[i].y-node[j].y));
					if(seperation < NODES_TOO_CLOSE)
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
				node[i].x = (double)rand()/((double)RAND_MAX)*(xStop - xStart) + xStart;
				node[i].y = (double)rand()/((double)RAND_MAX)*(yStop - yStart) + yStart;
				mag = sqrt(node[i].x*node[i].x + node[i].y*node[i].y);
				radius = ((double)rand()/(double)RAND_MAX)*maxRadius;
				node[i].x *= radius/mag;
				node[i].y *= radius/mag;
				repeatedNodeTest = 1;
				for(int j = 0; j < i; j++)
				{
					seperation = sqrt((node[i].x-node[j].x)*(node[i].x-node[j].x) + (node[i].y-node[j].y)*(node[i].y-node[j].y));  
					if(seperation < NODES_TOO_CLOSE)
					{
						repeatedNodeTest = 0;
						break;
					}
				}
			}
		}
	}
	//printf("\n %f  %f  %f  %f", node[63].x, node[63].y, node[2322].x, node[2322].y);
}

int checkNodes(float4 *node, int n)
{
	double seperation;
	
	for(int i = 0; i < (n - 1); i++)
	{
		for(int j = (i + 1); j < n; j++)
		{
			seperation = sqrt((node[i].x-node[j].x)*(node[i].x-node[j].x) + (node[i].y-node[j].y)*(node[i].y-node[j].y));
			if(seperation < NODES_TOO_CLOSE)
			{
			//printf("\n i = %d j = %d", i,j);
				return(-1);
			}
		}
	}
	return(1);
}

void placeNodesGrid(float4 *node, int rows, int columns)
{
	int i,j,k;
	double dx, dy;
	double xStart, yStart;
	
	xStart = 0.0;
	yStart = 0.0;
	dx = 1.0;
	dy = 1.0;
	
	k = 0;
	for(i = 0; i < columns; i++)
	{
		for(j = 0; j < rows; j++)
		{
			node[k].x = xStart + dx*i;
			node[k].y = yStart + dy*j;
			k++;
		}
	}
}

void getNumberOfNodesFromNodeFile(int *numberOfNodes, double *exhaustiveCost, char *nodeFileName)
{
	FILE *nodeFile;
	nodeFile = fopen(nodeFileName, "rb");
	fscanf(nodeFile,"%d", numberOfNodes);
	fscanf(nodeFile,"%lf", exhaustiveCost);
	fclose(nodeFile);
}

void placeNodesFromAFile(float4 *node, int *numberOfNodes, char *nodeFileName)
{
	FILE *nodeFile;
	double exhaustiveCost ;
	nodeFile = fopen(nodeFileName, "rb");
	fscanf(nodeFile,"%d", numberOfNodes);
	fscanf(nodeFile,"%lf", &exhaustiveCost);
	
	//printf("\n\n numberOfNodes = %d, exhaustiveCost = %lf\n\n", *numberOfNodes, exhaustiveCost);
	
	for(int i = 0; i < *numberOfNodes; i++)
	{
		//fscanf(nodeFile,"%d %f %f", &nodeNumber, &node[i].x, &node[i].y );
		fscanf(nodeFile,"%f %f", &node[i].x, &node[i].y );
		//printf("\n\n node number= %d, %f, %f\n\n", i, node[i].x, node[i].y);
	}
	
	fclose(nodeFile);
}

void getNumberOfNodesFromASelfDefinedFunction(int *numberOfNodes)
{
	*numberOfNodes = 12*16;
}

void placeNodesFromASelfDefinedFunction(float4 *node, int numberOfNodes)
{
	int i, j, k, count;
	float delta, shiftX, shiftY;
	
	delta = 8.4;
	count = 0;
	for(j = 0; j <2; j++)
	{
		shiftX = -4.2 + j*delta;
		for(k = 0; k <2; k++)
		{
			shiftY = -4.2 + k*delta;
			for(i = 0; i < 12; i++)
			{
				node[count].x = 2.0*cos(i*2.0*PI/(12.0)) + 2.1 + shiftX;
				node[count].y = 2.0*sin(i*2.0*PI/(12.0)) + 2.1 + shiftY;
				count++;
			}
			for(i = 0; i < 12; i++)
			{
				node[count].x = 2.0*cos(i*2.0*PI/(12.0)) - 2.1 + shiftX;
				node[count].y = 2.0*sin(i*2.0*PI/(12.0)) + 2.1 + shiftY;
				count++;
			}
			for(i = 0; i < 12; i++)
			{
				node[count].x = 2.0*cos(i*2.0*PI/(12.0)) - 2.1 + shiftX;
				node[count].y = 2.0*sin(i*2.0*PI/(12.0)) - 2.1 + shiftY;
				count++;
			}
			for(i = 0; i < 12; i++)
			{
				node[count].x = 2.0*cos(i*2.0*PI/(12.0)) + 2.1 + shiftX;
				node[count].y = 2.0*sin(i*2.0*PI/(12.0)) - 2.1 + shiftY;
				count++;
			}
		}
	}
	
	/*
	for(i = 0; i < 12; i++)
	{
		node[count].x = 2.0*cos(i*2.0*PI/(12.0)) + 2.1 + 8.1;
		node[count].y = 2.0*sin(i*2.0*PI/(12.0)) + 2.1;
		count++;
	}
	for(i = 0; i < 12; i++)
	{
		node[count].x = 2.0*cos(i*2.0*PI/(12.0)) - 2.1 + 8.1;
		node[count].y = 2.0*sin(i*2.0*PI/(12.0)) + 2.1;
		count++;
	}
	for(i = 0; i < 12; i++)
	{
		node[count].x = 2.0*cos(i*2.0*PI/(12.0)) - 2.1 + 8.1;
		node[count].y = 2.0*sin(i*2.0*PI/(12.0)) - 2.1;
		count++;
	}
	for(i = 0; i < 12; i++)
	{
		node[count].x = 2.0*cos(i*2.0*PI/(12.0)) + 2.1 + 8.1;
		node[count].y = 2.0*sin(i*2.0*PI/(12.0)) - 2.1;
		count++;
	}
	*/
}

//This function adjustes the nodes so that the average seperation is 1.0
double setAverageSeperationToOne(float4 *node, int numberOfNodes)
{
	double nodeAdjustmentFactor;
	double sum;
	int numberOfEdges;
	double dx,dy;
	int i,j;
	
	sum = 0.0;
	for(i = 0; i < numberOfNodes; i++)
	{
		for(j = i + 1; j < numberOfNodes; j++)
		{
			dx = node[i].x-node[j].x;
			dy = node[i].y-node[j].y;
			sum += sqrt(dx*dx + dy*dy);
		}
	}
	
	numberOfEdges = ((numberOfNodes)*(numberOfNodes - 1))/2;
	nodeAdjustmentFactor = sum/numberOfEdges;
	
	for(int i = 0; i < numberOfNodes; i++)
	{
		node[i].x = node[i].x/nodeAdjustmentFactor;
		node[i].y = node[i].y/nodeAdjustmentFactor;
	}
	
	return(nodeAdjustmentFactor);
}

//This function adjustes the nodes so that the minimum seperation is 1.0
double setMinimumSeperationToOne(float4 *node, int numberOfNodes)
{
	double nodeAdjustmentFactor;
	double minimum;
	double dx,dy, d;
	int i,j;
	
	minimum = BIG_NUMBER;
	for(i = 0; i < numberOfNodes; i++)
	{
		for(j = i + 1; j < numberOfNodes; j++)
		{
			dx = node[i].x-node[j].x;
			dy = node[i].y-node[j].y;
			d = sqrt(dx*dx + dy*dy);
			if(d < minimum) minimum = d;
		}
	}
	
	nodeAdjustmentFactor = minimum;
	
	for(int i = 0; i < numberOfNodes; i++)
	{
		node[i].x = node[i].x/nodeAdjustmentFactor;
		node[i].y = node[i].y/nodeAdjustmentFactor;
	}
	
	return(nodeAdjustmentFactor);
}

float4 setGeometricCenterToZero(float4 *node, int n)
{
	float4 geometricCenter;
	
	geometricCenter.x = 0.0;
	geometricCenter.y = 0.0;
	
	for(int i = 0; i < n; i++)
	{
		geometricCenter.x += node[i].x;
		geometricCenter.y += node[i].y;
	}
	
	geometricCenter.x /= (float)n;
	geometricCenter.y /= (float)n;
	
	for(int i = 0; i < n; i++)
	{
		node[i].x -= geometricCenter.x;
		node[i].y -= geometricCenter.y;
	}
	return(geometricCenter);
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

double findMinimumSeperation(float4 *pos, int n)
{
	int i,j;
	double min;
	double temp;
	
	min = BIG_NUMBER;
	
	for(i = 0; i < n; i++)
	{
		for(j = i + 1; j < n; j++)
		{
			temp = sqrt((pos[i].x-pos[j].x)*(pos[i].x-pos[j].x) + (pos[i].y-pos[j].y)*(pos[i].y-pos[j].y));
			if( temp < min) 
			{
				min = temp;
			}
		}
	}
	if(min < SMALL_NUMBER)
	{
		printf("\n  You have at least two nodes in the same location\n");
		printf("\n 	Minimum seperation = %f \n", min);
		printf("\n  Good Bye\n");
		exit(0);
	}
	return(min);
}

double findDistanceToOuterMostElement(float4 *element, int numberOfElements)
{
	double temp;
	double distanceToOutermostElement = 0.0;
	
	for(int i = 0; i < numberOfElements; i++)
	{
		temp = sqrt(element[i].x*element[i].x + element[i].y*element[i].y);
		if(temp > distanceToOutermostElement)
		{
			distanceToOutermostElement = temp;
		}
	}
	return(distanceToOutermostElement);
}

double findForceCPU(float4 node1, float4 node2, float4 pos1, float4 pos2)
{
	double naturalLength, actualLength, deltaLength, forceMag;
	
	naturalLength = sqrt((node1.x - node2.x)*(node1.x - node2.x) + (node1.y - node2.y)*(node1.y - node2.y));
	actualLength  = sqrt((pos1.x - pos2.x)*(pos1.x - pos2.x) + (pos1.y - pos2.y)*(pos1.y - pos2.y));
	deltaLength = actualLength - naturalLength;
	
	if(deltaLength <= 0.0)
	{
		forceMag = -deltaLength*SLOPE_REPULSION;
	}
	else if(actualLength < FORCE_CUTOFF)
	{
		forceMag =  MAG_ATRACTION/naturalLength;
	}
	else
	{
		forceMag = 0.0;
	}
	return(forceMag);
}

int findWallMoveDirections(float4 *node, float4 *pos, int n, int scope, float outerRadius, int *innerWallDirection, int *outerWallDirection)
{
	int i,j;
	int count;
	double sum, temp, maxValue, testValue;
	float lowerLimit, upperLimit;
	
	if(scope == 0) // No adjustment inner wall just moves out.
	{
		*innerWallDirection = 1;
		*outerWallDirection = 0;
		return(1);
	}

	if(scope == 1) // Move walls by pressure on the outer wall.
	{
		lowerLimit = LOWER_PRESSURE_LIMIT;
		upperLimit = UPPER_PRESSURE_LIMIT;
		sum = 0.0;
		count = 0;
		for(i = 0; i < n; i++)
		{
			temp = sqrt(pos[i].x*pos[i].x + pos[i].y*pos[i].y) - outerRadius;
			if(0.0 < temp) 
			{
				sum += temp;
			}
		}
		testValue = sum*WALL_STRENGTH/(2.0*PI*outerRadius);
	}
	
	if(scope == 2) // Move walls by average stress on the nodes.
	{
		lowerLimit = LOWER_AVERAGE_STRESS_LIMIT;
		upperLimit = UPPER_AVERAGE_STRESS_LIMIT;
		sum = 0.0;
		count = 0;
		for(i = 0; i < n; i++)
		{
			for(j = i + 1; j < n; j++)
			{
				sum += findForceCPU(node[i], node[j], pos[i], pos[j]);
				count++;
			}
		}
		testValue = sum/count;
	}
	
	if(scope == 3) // Move walls by absolute average stress.
	{
		lowerLimit = LOWER_AVERAGE_ABSOLUTE_STRESS_LIMIT;
		upperLimit = UPPER_AVERAGE_ABSOLUTE_STRESS_LIMIT;
		sum = 0.0;
		count = 0;
		for(i = 0; i < n; i++)
		{
			for(j = i + 1; j < n; j++)
			{
				sum += abs(findForceCPU(node[i], node[j], pos[i], pos[j]));
				count++;
			}
		}
		testValue = sum/count;
	}
	
	if(scope == 4)
	{
		lowerLimit = LOWER_AVERAGE_POSITIVE_STRESS_LIMIT;
		upperLimit = UPPER_AVERAGE_POSITIVE_STRESS_LIMIT;
		sum = 0.0;
		count = 0;
		for(i = 0; i < n; i++)
		{
			for(j = i + 1; j < n; j++)
			{
				temp = findForceCPU(node[i], node[j], pos[i], pos[j]);
				if(0.0 < temp) 
				{
					sum += temp;
					count++;
				}
			}
		}
		testValue = sum/count;
	}
	
	if(scope == 5)
	{
		lowerLimit = LOWER_MAX_POSITIVE_STRESS_FROM_OUTER_WALL_LIMIT;
		upperLimit = UPPER_MAX_POSITIVE_STRESS_FROM_OUTER_WALL_LIMIT;
		maxValue = 0.0;
		for(i = 0; i < n; i++)
		{
			if(0.0 <= sqrt(pos[i].x*pos[i].x + pos[i].y*pos[i].y) - outerRadius) 
			{
				for(j = i + 1; j < n; j++)
				{
					temp = findForceCPU(node[i], node[j], pos[i], pos[j]);
					if(maxValue < temp) 
					{
						maxValue = temp;
					}
				}
			}
		}
		testValue = maxValue;
	}
	
	if(testValue < lowerLimit)
	{
		*innerWallDirection = 0;
		*outerWallDirection = -1;
	}
	else if(testValue < upperLimit)
	{
		*innerWallDirection = 1;
		*outerWallDirection = 0;
	}
	else
	{
		*innerWallDirection = 0;
		*outerWallDirection = 1;
	}
		
	return(1);
}

void printEdgeCosts(float4 *node, int n, float nodeAdjustmentFactor)
{
	double temp;
	for(int i = 0; i < n; i++)
	{
		for(int j = i + 1; j < n; j++)
		{	
			temp = sqrt((node[i].x-node[j].x)*(node[i].x-node[j].x) + (node[i].y-node[j].y)*(node[i].y-node[j].y))*nodeAdjustmentFactor;
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

double getPathCost(int *path, float4 *node, int type, int n)
{
	double cost;
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
		cost += sqrt((node[i].x-node[j].x)*(node[i].x-node[j].x) + (node[i].y-node[j].y)*(node[i].y-node[j].y));
	}
	i = path[n-1];
	j = path[0];
	cost += sqrt((node[i].x-node[j].x)*(node[i].x-node[j].x) + (node[i].y-node[j].y)*(node[i].y-node[j].y));
	
	return(cost);
}

void swap(int *path, int i, int j)
{
	int temp;
	temp = path[i];
	path[i] = path[j];
	path[j] = temp;
}

void heappermute(int *path, int m, float4 *node, int *exhaustivePath, double *minCost, int n) 
{
	int i;
	double pathCost;
	int* pathPlus = (int*)malloc(n*sizeof(int));

	if (m == 1) 
	{
		pathPlus[0] = 0;
		for(i = 1; i < n; i++)
		{
			pathPlus[i] = path[i-1];	
		}
		
		pathCost = getPathCost(pathPlus, node, 1, n);
		
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
			heappermute(path, m-1, node, exhaustivePath, minCost, n);
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

double exhaustiveTSP(float4 *node, int* exhaustivePath, int n)
{
	double cost[1];
	int* path = (int*)malloc((n-1)*sizeof(int));
	
	exhaustivePath[0] = 0;
	for(int i = 1; i < n; i++)
	{
		exhaustivePath[i] = i;
		path[i-1] = i;	
	}
	cost[0] = getPathCost(exhaustivePath, node, 1, n);
	
	heappermute(path, n-1, node, exhaustivePath, cost, n);
	free(path);
	return(cost[0]);
}

double nearestNeighborTSP(float4 *node, int* path, int n)
{
	int i, j, k, nextNode, nodeFound;
	double minCost, pathCost, edgeCost, maxEdgeCost;
	int* used = (int*)malloc(n*sizeof(int));
	
	maxEdgeCost = 0.0;
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			edgeCost = sqrt((node[i].x-node[j].x)*(node[i].x-node[j].x) + (node[i].y-node[j].y)*(node[i].y-node[j].y));
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
			edgeCost = sqrt((node[i].x-node[j].x)*(node[i].x-node[j].x) + (node[i].y-node[j].y)*(node[i].y-node[j].y));
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
	
	pathCost = getPathCost(path, node, 2, n);
	free(used);
	return(pathCost);
}

void setNbodyInitailConditions(float4 *node, float4 *pos, float4 *vel, float* mass, int n)
{
	int i;

	for(i = 0; i < n; i++)
	{
		pos[i].x = node[i].x;
		pos[i].y = node[i].y;
		
		vel[i].x = 0.0;
		vel[i].y = 0.0;
		
		mass[i] = MASS;
	}
}

__device__ float getBodyBodyForceMag(float edgeLength, float d)
{
	float forceMag;
	float c, h;
	
	if(FORCE_FUNTION_TYPE == 0)  // Lenard Jones type force
	{
		forceMag = -H*powf(edgeLength, Q - P)/powf(d, Q) + H/powf(d, P);
	}
	if(FORCE_FUNTION_TYPE == 1)  // Lenard Jones type force
	{
		//d = dist(shPos[i], pos) + 0.01;

		h = M*(powf(powf(Q/P, 1/(Q-P))*edgeLength, P))/(1 - P/Q);
		
		c = powf(edgeLength/d, Q - P);
		forceMag = (c-1)*h/powf(d, P);
	}
	else if(FORCE_FUNTION_TYPE == 2)  // Linear function
	{
		if(d <= edgeLength)
		{
			forceMag = -(edgeLength - d)*SLOPE_REPULSION;

		}
		else if(edgeLength < d && d < FORCE_CUTOFF)
		{
			forceMag =  MAG_ATRACTION/edgeLength;
		}
		else
		{
			forceMag = 0.0;
		}
	}
	else if(FORCE_FUNTION_TYPE == 3)  // Billards function
	{
		if(d <= NODE_RADIUS)
		{
			forceMag = -edgeLength*(NODE_RADIUS - d);

		}
		else
		{
			forceMag = 0.0;
		}
	}
	else
	{
		forceMag = 0.0;
	}
	
	return(forceMag);
	
	
}

__global__ void accelerationsNbody(float4 *node, float4 *pos, float4 *vel, float4 *acc, float *mass, float innerRadius, float outerRadius, int n)
{
	int j,ii;
    float3 forceSum;
    float4 nodeMe, posMe;
    float dx, dy, d, edgeLength; 
    float radius, forceMag;
    __shared__ float4 shnode[BLOCK], shPos[BLOCK];
    int id = threadIdx.x + blockDim.x*blockIdx.x;
    
    forceSum.x = 0.0;
	forceSum.y = 0.0;
	
	nodeMe.x = node[id].x;
	nodeMe.y = node[id].y;
	posMe.x = pos[id].x;
	posMe.y = pos[id].y;
		    
    for(j=0; j < gridDim.x; j++)
    {
    	if(threadIdx.x + blockDim.x*j < n)
    	{
    		shPos[threadIdx.x] = pos[threadIdx.x + blockDim.x*j];
    		shnode[threadIdx.x] = node[threadIdx.x + blockDim.x*j];
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
				
				edgeLength = sqrtf((shnode[i].x - nodeMe.x)*(shnode[i].x - nodeMe.x) + (shnode[i].y - nodeMe.y)*(shnode[i].y - nodeMe.y));
				
				forceMag = getBodyBodyForceMag(edgeLength, d);
				forceSum.x += forceMag*dx/d;
				forceSum.y += forceMag*dy/d;
		    }
		}
	}
	
	if(id < n)
	{
		// Forces between node and the walls
		dx = posMe.x;
		dy = posMe.y; 
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

void getPathNbody(float4 *pos, int* path, int n)
{
	int i;
	double minValue;
	double *angle = (double*)malloc(n*sizeof(double));
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

double NbodyExtrusionTSP(float4 *node, float4 *pos, float4 *vel, float4 *acc, float* mass, int* path, int n)
{
	int draw_count;
	int innerWallDirection, outerWallDirection;
	int annealingCount;
	double dr;
	float dt = TIME_STEP_SIZE;
	double pathCost;
	double time;
	float innerRadius, outerRadius;
	double stopSeperation;
	
	dim3 block, grid;
	float4 *posGPU, *velGPU, *accGPU; 
	float *massGPU;
	float4 *nodeGPU;
	
	// Setting up GPU parrellel structure.
	block.x = BLOCK;
	block.y = 1;
	block.z = 1;
	
	grid.x = (n-1)/block.x + 1;
	grid.y = 1;
	grid.z = 1;
	
	// Allocating memory.
	cudaMalloc( (void**)&nodeGPU, n *sizeof(float4));
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
	cudaMemcpy( nodeGPU, node, n *sizeof(float4), cudaMemcpyHostToDevice );
	cudaMemcpy( posGPU, pos, n *sizeof(float4), cudaMemcpyHostToDevice );
    cudaMemcpy( velGPU, vel, n *sizeof(float4), cudaMemcpyHostToDevice );
    cudaMemcpy( massGPU, mass, n *sizeof(float), cudaMemcpyHostToDevice );
	
	annealingCount = 0;
	while(annealingCount <= NUMBER_OF_ANEALINGS)
	{
		moveAnyNodeOffDeadCenter(pos, n);
		stopSeperation = findMinimumSeperation(pos, n)/1.0;
		innerRadius = 0.0;
		outerRadius = findDistanceToOuterMostElement(pos, n);
		drawNbodyExtrusion(pos, innerRadius, outerRadius, 0, 0, n);
		
		outerWallDirection = 0;
		innerWallDirection = 0;
		dr = outerRadius/STARTING_POINT_FOR_NUMBER_OF_MOVES;
		draw_count = 0;
	
		while(innerRadius + stopSeperation < outerRadius)
		//while(innerRadius + stopSeperation < dr)
		{
			outerRadius += dr*outerWallDirection;
			innerRadius += dr*innerWallDirection;
			time = 0.0;
			while(time < TIME_BETWEEN_WALL_MOVES)
			{		
				accelerationsNbody<<<grid, block>>>(nodeGPU, posGPU, velGPU, accGPU, massGPU, innerRadius, outerRadius, n);
				moveNbody<<<grid, block>>>(posGPU, velGPU, accGPU, dt, n);
			
				if(draw_count == STEPS_BETWEEN_VIEWING)
				{
					cudaMemcpy( pos, posGPU, n *sizeof(float4), cudaMemcpyDeviceToHost );
					drawNbodyExtrusion(pos, innerRadius, outerRadius, innerWallDirection, outerWallDirection, n);
					draw_count = 0;
				}
				draw_count++;
				time += dt;
			}
			cudaMemcpy( pos, posGPU, n *sizeof(float4), cudaMemcpyDeviceToHost );
			findWallMoveDirections(node, pos, n, WALL_MOVE_SCOPE, outerRadius, &innerWallDirection, &outerWallDirection);
			dr = outerRadius/STARTING_POINT_FOR_NUMBER_OF_MOVES;
		}
		
		if(annealingCount < NUMBER_OF_ANEALINGS)
		{
			innerRadius = 0.0;
			outerRadius = outerRadius + outerRadius*0.7;
			time = 0.0;
			while(time < ANNEALING_TIME)
			{		
				accelerationsNbody<<<grid, block>>>(nodeGPU, posGPU, velGPU, accGPU, massGPU, innerRadius, outerRadius, n);
				moveNbody<<<grid, block>>>(posGPU, velGPU, accGPU, dt, n);
			
				if(draw_count == STEPS_BETWEEN_VIEWING)
				{
					cudaMemcpy( pos, posGPU, n *sizeof(float4), cudaMemcpyDeviceToHost );
					drawNbodyExtrusion(pos, innerRadius, outerRadius, innerWallDirection, outerWallDirection, n);
					draw_count = 0;
				}
				draw_count++;
				time += dt;
			}
		}
		annealingCount++;
	}
	getPathNbody(pos, path, n);
	pathCost = getPathCost(path, node, 3, n);
	
	return(pathCost);
}

void drawPoints(float4 *pos, int n)
{
	int i;
	float outerRadius = findDistanceToOuterMostElement(pos, n);
	float normalizingFactor = outerRadius;
	
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

void drawNbodyExtrusion(float4 *pos, float innerRadius, float outerRadius, int innerWallDirection, int outerWallDirection, int n)
{
	int i;
	int lineAmount = 100;
	float normalizingFactor = outerRadius;
	
	outerRadius /= normalizingFactor;
	innerRadius /= normalizingFactor;

	glClear(GL_COLOR_BUFFER_BIT);
	
	GLfloat twicePi = 2.0f * PI;
	
	glLineWidth(1.0);
	if(innerWallDirection == -1) glColor3f(1.0,0.0,0.0);
	else if(innerWallDirection == 0) glColor3f(1.0,1.0,0.0);
	else glColor3f(0.0,0.0,1.0);
	glBegin(GL_LINE_LOOP);
		for(i = 0; i <= lineAmount;i++) 
		{ 
			glVertex2f(x_world_to_x_screen(innerRadius*cos(i*twicePi/lineAmount)), 
			           y_world_to_y_screen(innerRadius*sin(i*twicePi/lineAmount)));
		}
	glEnd();
	
	glLineWidth(1.0);
	if(outerWallDirection == -1) glColor3f(1.0,0.0,0.0);
	else if(outerWallDirection == 0) glColor3f(1.0,1.0,0.0);
	else glColor3f(0.0,0.0,1.0);
	glBegin(GL_LINE_LOOP);
		for(i = 0; i <= lineAmount;i++) 
		{ 
			glVertex2f(x_world_to_x_screen(outerRadius*cos(i*twicePi/lineAmount)), 
			           y_world_to_y_screen(outerRadius*sin(i*twicePi/lineAmount)));
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

void drawFInalPicture(float4 *node, int *pathA, int *pathB, int *pathC, int scope, int n)
{	
	int i;
	float outerRadius = findDistanceToOuterMostElement(node, n);
	float normalizingFactor = outerRadius; //((float)n)/IDEAL_NUMBER_OF_NODES;

	glClear(GL_COLOR_BUFFER_BIT);
	
	//exhuastivePath path
	if(scope == 1 || scope == 2)
	{
		if(DRAW_EXHAUSTIVE_PATH == 1)
		{
			glLineWidth(6.0);
			glColor3f(0.0,0.0,1.0);
			glBegin(GL_LINE_LOOP);
				for(i = 0; i < n; i++)
				{
					glVertex2f(x_world_to_x_screen(node[pathA[i]].x/normalizingFactor),y_world_to_y_screen(node[pathA[i]].y/normalizingFactor));
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
				glVertex2f(x_world_to_x_screen(node[pathB[i]].x/normalizingFactor),y_world_to_y_screen(node[pathB[i]].y/normalizingFactor));
			}
		glEnd();
	}
	
	//Nbody Extrusion path
	if(DRAW_NBODY_EXTRUSION_PATH == 1)
	{
		glLineWidth(3.0);
		glColor3f(1.0,0.0,0.0);
		glBegin(GL_LINE_LOOP);
			for(i = 0; i < n; i++)
			{
				glVertex2f(x_world_to_x_screen(node[pathC[i]].x/normalizingFactor),y_world_to_y_screen(node[pathC[i]].y/normalizingFactor));
			}
		glEnd();
	}
	
	//Placing nodes
	glPointSize(8.0);
	glColor3f(1.0,1.0,1.0);
	for(i = 0; i < n; i++)
	{
		glBegin(GL_POINTS);
			glVertex2f(x_world_to_x_screen(node[i].x/normalizingFactor),y_world_to_y_screen(node[i].y/normalizingFactor));
		glEnd();
	}
	
	//Nearest neighbor start node 
	if(DRAW_NEAREST_NEIGHBOR_PATH == 1)
	{
		glPointSize(10.0);
		glColor3f(0.0,0.0,1.0);
		glBegin(GL_POINTS);
			glVertex2f(x_world_to_x_screen(node[pathB[0]].x/normalizingFactor),y_world_to_y_screen(node[pathB[0]].y/normalizingFactor));
		glEnd();
	}
	
	//Nbody extrution start and stop nodes
	if(DRAW_NBODY_EXTRUSION_PATH == 1)
	{
		glPointSize(10.0);
		glColor3f(0.0,1.0,0.0);
		glBegin(GL_POINTS);
			glVertex2f(x_world_to_x_screen(node[pathC[0]].x/normalizingFactor),y_world_to_y_screen(node[pathC[0]].y/normalizingFactor));
		glEnd();
	
		glColor3f(1.0,0.0,0.0);
		glBegin(GL_POINTS);
			glVertex2f(x_world_to_x_screen(node[pathC[n-1]].x/normalizingFactor),y_world_to_y_screen(node[pathC[n-1]].y/normalizingFactor));
		glEnd();
	}
	
	glFlush();
}

void getInputFromUser(int* scope, int* numberOfNodes, int* numberOfRuns, int* maxNumberOfRows, int* maxNumberOfColumns, unsigned int* srandSeed, char *nodeFileName)
{
	*scope = -1;
	*numberOfNodes = -1;
	*numberOfRuns = -1;
	*maxNumberOfRows = -1;
	*maxNumberOfColumns = -1;
	
	printf("\n\n  What type run would you like to perform?");
	printf("\n  1 for one small randomly generated run.");
	printf("\n  2 for a series of small randomly generated runs.");
	printf("\n  3 for one on on a grid.");
	printf("\n  4 for a series of runs on randomly generated sized grids.");
	printf("\n  5 for one large randomly generated run.");
	printf("\n  6 for a series of large randomly generated runs.");
	printf("\n  7 to read nodes from nodeFile.");
	printf("\n  8 for a set created an self created function.");
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
		printf("\n  Note: If you choose a number bigger than 13 you may lock your computer up.");
		
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
		
		printf("\n\n  You will also need to enter a seed to generate the random placed nodes");
		printf("\n  By intering the same seed on different exicutions of the program");
		printf("\n  you will be able to run multiple tests on the same node configuration.");
		
		printf("\n  Enter two positive integer values seperated by a space: \n  (number of nodes) (seed)");
		printf("\n\n  Inter your values: ");
		scanf("%d %d", numberOfNodes, srandSeed);
		
		*numberOfRuns = 1;
	}
	else if(*scope == 6)
	{
		printf("\n\n  You will need to enter the number nodes (Cities) wound like to generate?");
		
		printf("\n\n  You will also need to enter how many randomly generated runs would you like to perform?");
		
		printf("\n  Enter two positive integer values seperated by a space: \n  (number of nodes) (number of runs)");
		printf("\n\n  Inter your values: ");
		scanf("%d %d", numberOfNodes, numberOfRuns);
	}
	else if(*scope == 7)
	{	
		*numberOfRuns = 1;
		
		printf("\n\n  You will need to enter the name of the file containing the nodes");
		
		printf("\n\n  The first element of the file should be an intiger containing the number of nodes.");
		printf("\n  The second element should be the value of the minumum path. This should be -1 if the value is unknown.");
		printf("\n  The following lines will contain the x y positions of all the nodes.");
		printf("\n\n  ");
		scanf("%s" , nodeFileName);
	
	}
	else if(*scope == 8)
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
	float4 *node;
	double nodeAdjustmentFactor;
	float4 geometricCenter;
	double distanceToOutermostNode;
	int *exhaustivePath, *nearestNeighborPath, *NbodyExtrusionPath;
	float4 *posNbody, *velNbody, *accNbody; 
	float *massNbody;
	double exhaustiveCost, nearestNeighborCost, NbodyExtrusionCost;
	int nodeCheck;
	double temp;
	char nodeFileName[100];
	
	getInputFromUser(&scope, &numberOfNodes, &numberOfRuns, &maxNumberOfRows, &maxNumberOfColumns, &srandSeed, nodeFileName);
	
	if(scope == 2 || scope == 4 || scope == 6 && PRINT_RAW_DATA_FILE == 1)
	{
		openRawDataFile(scope, numberOfRuns);
	}
	
	double totalNearestNeighborCost = 0.0;
	double totalNbodyExtrusionCost = 0.0;
	double totalPercentErrorNearestNeighbor = 0.0;
	double totalPercentErrorNbodyExtrusion = 0.0;
	double NbodyExtrusionVSNearestNeighbor = 0.0;
	
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
		if(scope == 7)
		{
			getNumberOfNodesFromNodeFile(&numberOfNodes, &exhaustiveCost, nodeFileName);
		}
		if(scope == 8)
		{
			getNumberOfNodesFromASelfDefinedFunction(&numberOfNodes);
		}
		
		node = (float4*)malloc((numberOfNodes)*sizeof(float4));
	
		exhaustivePath = (int*)malloc((numberOfNodes)*sizeof(int));    // !!!!!!!!!!!! only needed in scope 1 and 2
		nearestNeighborPath = (int*)malloc((numberOfNodes)*sizeof(int));
		NbodyExtrusionPath = (int*)malloc((numberOfNodes)*sizeof(int));

		posNbody = (float4*)malloc((numberOfNodes)*sizeof(float4));
		velNbody = (float4*)malloc((numberOfNodes)*sizeof(float4));
		accNbody = (float4*)malloc((numberOfNodes)*sizeof(float4));
		massNbody = (float*)malloc((numberOfNodes)*sizeof(float4));
		
		//Creating nodes
		if(scope == 1 || scope == 2 || scope == 5 || scope == 6)
		{	
			placeNodesRandom(node, srandSeed, scope, numberOfNodes);
		}
		else if(scope == 3 || scope == 4)
		{
			placeNodesGrid(node, rows, columns);
		}
		else if(scope == 7)
		{
			placeNodesFromAFile(node, &numberOfNodes, nodeFileName);
		}
		else if(scope == 8)
		{
			placeNodesFromASelfDefinedFunction(node, numberOfNodes);
		}
		
		//Adjusting nodes
		geometricCenter = setGeometricCenterToZero(node, numberOfNodes);
		printf("\n\n  The geometric center of the nodes = (%f, %f)", geometricCenter.x, geometricCenter.y);
		
		distanceToOutermostNode = findDistanceToOuterMostElement(node, numberOfNodes);
		printf("\n  The distance to the outermost node from the geometric center pre adjustment is %f", distanceToOutermostNode);
		
		//nodeAdjustmentFactor = setAverageSeperationToOne(node, numberOfNodes);
		//nodeAdjustmentFactor = setMinimumSeperationToOne(node, numberOfNodes);
		nodeAdjustmentFactor = 1.0;
		printf("\n  The node adjustment factor = %f", nodeAdjustmentFactor);
		
		distanceToOutermostNode = findDistanceToOuterMostElement(node, numberOfNodes);
		printf("\n  The distance to the outermost node from the geometric center post adjustment is %f", distanceToOutermostNode);
		
		//Checking to see if a node is repeated
		nodeCheck = checkNodes(node, numberOfNodes);
		if(nodeCheck == -1)
		{
			printf("\n\n  There is a repeated node. Check your data set.");
			printf("\n\n  Good Bye.  \n\n");
			exit(0);
		}
		
		//Drawing the adjusted nodes on the screen.
		drawPoints(node, numberOfNodes); 
		
		//Printing the edge costs (lengths in this case)
		if(PRINT_EDGE_COST == 1)
		{
			printEdgeCosts(node, numberOfNodes, nodeAdjustmentFactor);
		}
		
		//Finding exact cost
		printf("\n\n  Determining the exact cost.");
		if(scope == 1 || scope == 2)
		{	
			exhaustiveCost = exhaustiveTSP(node, exhaustivePath, numberOfNodes);
		}
		else if(scope == 3 || scope == 4)
		{
			//Assuming all edges are the same length. So just get the length of the first edge.
			temp = sqrt((node[0].x-node[1].x)*(node[0].x-node[1].x) + (node[0].y-node[1].y)*(node[0].y-node[1].y));
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
		else if(scope == 5 || scope == 6)
		{
			exhaustiveCost = -1.0;
		}
		else if(scope == 7)
		{
			if(numberOfNodes < 13)
			{
				exhaustiveCost = exhaustiveTSP(node, exhaustivePath, numberOfNodes);
			}
		}
		else if(scope == 8)
		{
			if(numberOfNodes < 13)
			{
				exhaustiveCost = exhaustiveTSP(node, exhaustivePath, numberOfNodes);
			}
			else
			{
				exhaustiveCost = -1.0;
			}
		}
		printf("\n  Determining the exact cost is done.");
		
		//Finding nearest neighbor cost
		printf("\n\n  Running the nearest nieghbor algorithm.");
		nearestNeighborCost = nearestNeighborTSP(node, nearestNeighborPath, numberOfNodes);
		printf("\n  The nearest nieghbor algorithm is done.");
		
		//Running n-body extrusion code
		printf("\n\n  Running the N-body extrusion algorithm."); 
		printf("  \n"); //I had to enter this carage return so it would print the line above before it started the algorithm
		setNbodyInitailConditions(node, posNbody, velNbody, massNbody, numberOfNodes);
		NbodyExtrusionCost = NbodyExtrusionTSP(node, posNbody, velNbody, accNbody, massNbody, NbodyExtrusionPath, numberOfNodes);
		printf("  The N-body extrusion algorithm is done.");
		
		//Unadjusting costs
		exhaustiveCost *= nodeAdjustmentFactor;
		nearestNeighborCost *= nodeAdjustmentFactor;
		NbodyExtrusionCost *= nodeAdjustmentFactor;
		
		totalNearestNeighborCost += nearestNeighborCost;
		totalNbodyExtrusionCost += NbodyExtrusionCost;
		
		//Sanity check
		if(nearestNeighborCost < exhaustiveCost - SMALL_NUMBER)
		{
			printf("\n\n  Nearest Neighbor cost (%f) is smaller than exhaustive cost (%f). Something is wrong!\n",nearestNeighborCost, exhaustiveCost);
			printf("\n\n  Good Bye.  \n\n");
			exit(0);
		}
		if(NbodyExtrusionCost < exhaustiveCost - SMALL_NUMBER)
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
		
		if(scope == 2 || scope == 6 && PRINT_RAW_DATA_FILE == 1)
		{
			fprintf(RawDataFile, "  %d, %d, %f, %f, %f\n", i+1, numberOfNodes, exhaustiveCost, nearestNeighborCost, NbodyExtrusionCost);
		}
		if(scope == 4 && PRINT_RAW_DATA_FILE == 1)
		{
			fprintf(RawDataFile, "  %d, %d, %d, %f, %f, %f\n", i+1, rows, columns, exhaustiveCost, nearestNeighborCost, NbodyExtrusionCost);
		}
	
		drawFInalPicture(node, exhaustivePath, nearestNeighborPath, NbodyExtrusionPath, scope, numberOfNodes);
		
		free(node);
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



    

