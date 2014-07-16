//
//  CVConcentricStructure.c
//  CVNetwork
//
//  Created by Filipi Nascimento Silva on 26/04/14.
//  Copyright (c) 2014 Filipi Nascimento Silva. All rights reserved.
//

#include "CVConcentricStructure.h"


CVConcentricStructure* CVNewConcentricStructureForNetwork(const CVNetwork* network, CVBool trackNodes){
	CVConcentricStructure* structure = calloc(1, sizeof(CVConcentricStructure));
	structure->distances = calloc(network->verticesCount, sizeof(CVSize));
	structure->_visitedNodes = CVNewBitArray(network->verticesCount);
	structure->network = network;
	if(trackNodes){
		structure->totalVerticesCapacity = 2;
		structure->vertices = calloc(structure->totalVerticesCapacity, sizeof(CVIndex));
		structure->concentricIndex = calloc(network->verticesCount, sizeof(CVSize));
	}
	
	structure->levelsCapacity = 2;
	structure->levelsCount = 1;
	structure->levelsIndices = calloc(structure->levelsCapacity+1, sizeof(CVIndex));
	structure->levelsIndices[0] = 0;
	structure->levelsIndices[1] = 0;
	return structure;
}

void _CV_ConcentricStructureAddVertexAtDistance(CVIndex vertexIndex, CVSize distance, CVConcentricStructure* structure){
	if(structure->vertices){
		structure->totalVerticesCount++;
		if(structure->totalVerticesCapacity<structure->totalVerticesCount){
			structure->totalVerticesCapacity = CVCapacityGrow(structure->totalVerticesCount+1);
			structure->vertices = realloc(structure->vertices, sizeof(CVIndex)*(structure->totalVerticesCapacity));
		}
		structure->vertices[structure->totalVerticesCount-1] = vertexIndex;
	}
	if(distance>structure->levelsCount-1){
		structure->levelsCount++;
		if(structure->levelsCapacity < structure->levelsCount){
			structure->levelsCapacity = CVCapacityGrow(structure->levelsCount+1);
			structure->levelsIndices = realloc(structure->levelsIndices, sizeof(CVIndex)*(structure->levelsCapacity+1));
		}
		structure->levelsIndices[structure->levelsCount] = structure->levelsIndices[structure->levelsCount-1];
	}
	if(structure->concentricIndex) {
		structure->concentricIndex[vertexIndex]= structure->levelsIndices[structure->levelsCount]-structure->levelsIndices[structure->levelsCount-1];
	}
	structure->levelsIndices[structure->levelsCount]++;
	structure->distances[vertexIndex] = distance;
}

void CVConcentricStructureSetReferenceVertex(CVIndex referenceVertex, CVSize maximumConcentricLevel, CVConcentricStructure* structure){
	const CVNetwork* network = structure->network;
	memset(structure->distances,0,network->verticesCount*sizeof(CVSize));
	
	
	//CVSize* distances = structure->distances;
	structure->levelsCount = 1;
	structure->levelsIndices[0] = 0;
	structure->levelsIndices[1] = 0;
	structure->totalVerticesCount = 0;
	structure->referenceVertex = referenceVertex;
	
	CVQueue toVisitQueue = CVQueueCreate();
	CVInteger currentVertex = referenceVertex;
	CVQueuePush(&toVisitQueue, currentVertex);
	CVBitArraySet(structure->_visitedNodes, currentVertex);
	_CV_ConcentricStructureAddVertexAtDistance(currentVertex, 0, structure);
	while (CVQueueDequeue(&toVisitQueue,&currentVertex)) {
		CVSize currentDistance = structure->distances[currentVertex];
		if(currentDistance<maximumConcentricLevel){
			CVSize vertexEdgesCount = network->vertexNumOfEdges[currentVertex];
			CVIndex* vertexEdgesList = network->vertexEdgesLists[currentVertex];
			CVIndex ni;
			for(ni=0;ni<vertexEdgesCount;ni++){
				CVIndex neighborVertex = vertexEdgesList[ni];
				if(!CVBitArrayTest(structure->_visitedNodes, neighborVertex)){
					_CV_ConcentricStructureAddVertexAtDistance(neighborVertex, currentDistance+1, structure);
					CVBitArraySet(structure->_visitedNodes, neighborVertex);
					CVQueuePush(&toVisitQueue, neighborVertex);
				}
			}
		}
	}
	CVQueueDestroy(&toVisitQueue);
	CVBitArrayClearAll(structure->_visitedNodes, network->verticesCount);
}


void CVConcentricStructurePrint(const CVConcentricStructure* structure){
	CVIndex l;
	CVSize levelsCount = structure->levelsCount;
	printf("Vertex %"CVUIntegerScan":\n",structure->referenceVertex);
	printf("\t v = [");
	for (l=0; l<structure->totalVerticesCount;l++){
		if (l) {
			printf(" ");
		}
		printf("%"CVSizeScan,structure->vertices[l]);
	}
	printf("]\n");
	
	
	if(structure->concentricIndex){
		printf("\t c = [");
		for (l=0; l<structure->totalVerticesCount;l++){
			if (l) {
				printf(" ");
			}
			printf("%"CVSizeScan,structure->concentricIndex[structure->vertices[l]]);
		}
		printf("]\n");
		
	}
	
	printf("\t h = [");
	for (l=0; l<structure->levelsCount+1;l++){
		if (l) {
			printf(" ");
		}
		printf("%"CVSizeScan,structure->levelsIndices[l]);
	}
	printf("]\n");
	if(structure->vertices){
		for (l=0; l<levelsCount;l++){
			CVSize concentricCount = 0;
			const CVIndex* concentricVertices = CVConcentricVerticesAtLevel(l, &concentricCount, structure);
			CVIndex i;
			printf("\t- l%"CVIndexScan": ",l);
			for (i=0; i<concentricCount; i++) {
				if (i) {
					printf(" ");
				}
				printf("%"CVIndexScan"",concentricVertices[i]);
			}
			printf("\n");
		}
	}else{
		for (l=0; l<levelsCount;l++){
			CVSize concentricCount = CVConcentricCountVerticesAtLevel(l, structure);
			printf("\t- l%"CVSizeScan": %"CVSizeScan" vertices.",l,concentricCount);
			printf("\n");
		}
	}
}


void _CV_ConcentricBackbonePropagateProbabilities(CVFloat* probabilities, CVSize* pahtsCount, CVSize level,CVSize* deadEndsCount, const CVConcentricStructure* structure){
	CVSize concentricCount = 0;
	const CVIndex* concentricVertices = CVConcentricVerticesAtLevel(level, &concentricCount, structure);
	CVIndex i;
	const CVNetwork* network = structure->network;
	const CVSize* distances = structure->distances;
	for (i=0; i<concentricCount; i++) {
		CVIndex vertexIndex = concentricVertices[i];
		CVSize vertexEdgesCount = network->vertexNumOfEdges[vertexIndex];
		CVIndex* vertexEdgesList = network->vertexEdgesLists[vertexIndex];
		CVIndex ni;
		CVSize betweenLevelConnections = 0;
		for(ni=0;ni<vertexEdgesCount;ni++){
			CVIndex neighborVertex = vertexEdgesList[ni];
			if(distances[neighborVertex]>level){
				betweenLevelConnections++;
			}
		}
		if(betweenLevelConnections==0){ //No connections to next level (dead end)
			(*deadEndsCount)++;
		}else{
			CVFloat propagateFactor = probabilities[vertexIndex]/betweenLevelConnections;
			CVSize pathsFactor = pahtsCount[vertexIndex];
			for(ni=0;ni<vertexEdgesCount;ni++){
				CVIndex neighborVertex = vertexEdgesList[ni];
				if(distances[neighborVertex]>level){
					probabilities[neighborVertex] += propagateFactor;
					pahtsCount[neighborVertex] += pathsFactor;
				}
			}
		}
	}
}

void CVConcentricBackboneGetProbabilities(CVFloat* probabilities, CVSize* pahtsCount, CVSize* deadEndCounts, CVSize maximumConcentricLevel, const CVConcentricStructure* structure){
	memset(probabilities,0, sizeof(CVFloat)*structure->network->verticesCount);
	memset(pahtsCount,0, sizeof(CVSize)*structure->network->verticesCount);
	
	probabilities[structure->referenceVertex] = 1.0f;
	pahtsCount[structure->referenceVertex] = 1;
	
	CVIndex l;
	CVSize levelCount = CVMIN(structure->levelsCount, maximumConcentricLevel+1);
	if(deadEndCounts){
		deadEndCounts[0]=0;
	}
	for (l=0; l<levelCount-1; l++) {
		CVSize deadEndCount = 0;
		_CV_ConcentricBackbonePropagateProbabilities(probabilities, pahtsCount, l, &deadEndCount, structure);
		if(deadEndCounts){
			deadEndCounts[l+1] = deadEndCounts[l]+deadEndCount;
		}
	}
}


void CVConcentricStructureDestroy(CVConcentricStructure* structure){
	free(structure->distances);
	free(structure->concentricIndex);
	
	free(structure->vertices);
	free(structure->levelsIndices);
	
	CVBitArrayDestroy(structure->_visitedNodes);
	
	free(structure);
}




/*
 edges					[0,1,2,0,1,2]
 l: 0
 0-> 0,1,2
 l: 1
 0-> []
 1-> [0,1]
 2-> [2]
 
 edgesIndices;			[[0,3],[3,3,5,6]]
 levelsEdgesIndices;	[0,2,6]
 
 mergedTranslate
 */

/*
 edges					[0,1,2,0,1,2]
 l: 0
 0-> 0,1,2
 l: 1
 0-> []
 1-> [0,1]
 2-> [2]
 
 edgesIndices;			[[0,3],[3,3,5,6]]
 levelsEdgesIndices;	[0,2,6]
 
 mergedTranslate
 */


void CVConcentricMergedInformationPrint(CVConcentricMergedInformation* mergedInformation){
	CVIndex l;
	printf("Merged Info for Vertex %"CVUIntegerScan":\n",mergedInformation->concentricStructure->referenceVertex);
	/*
	 printf("\t m = [");
	 for(i=0;i<mergedInformation->mergedVertices.count;i++){
	 if (i) {
	 printf(" ");
	 }
	 printf("%"CVIndexScan"",mergedInformation->mergedVertices.data[i]);
	 }
	 printf("]\n");
	 
	 printf("\tmi = [");
	 for(i=0;i<mergedInformation->mergedIndices.count;i++){
	 if (i) {
	 printf(" ");
	 }
	 printf("%"CVIndexScan"",mergedInformation->mergedIndices.data[i]);
	 }
	 printf("]\n");
	 
	 printf("\t h = [");
	 for (i=0; i<mergedInformation->concentricStructure->levelsCount+1;i++){
	 if (i) {
	 printf(" ");
	 }
	 printf("%"CVSizeScan"",mergedInformation->concentricStructure->levelsIndices[i]);
	 }
	 printf("]\n");
	 
	 for (l=0; l<mergedInformation->levelsCount; l++) {
	 printf("Level : %"CVIndexScan":\n",l);
	 CVSize concentricCount = CVConcentricCountVerticesAtLevel(l, mergedInformation->concentricStructure);
	 CVIndex concentricIndex;
	 for (concentricIndex=0; concentricIndex<concentricCount; concentricIndex++) {
	 printf("\t - [");
	 CVSize inMergedCount = 0;
	 CVIndex* groupsIndices = CVConcentricMergedConcentricIndices(concentricIndex, l, &inMergedCount, mergedInformation);
	 CVIndex mergedConcentricIndex;
	 for (mergedConcentricIndex=0; mergedConcentricIndex<inMergedCount; mergedConcentricIndex++) {
	 if (mergedConcentricIndex) {
	 printf(" ");
	 }
	 printf("%"CVIndexScan"", groupsIndices[mergedConcentricIndex]);
	 }
	 printf("]\n");
	 }
	 }
	 */
	
	for (l=0; l<mergedInformation->levelsCount-1; l++) {
		printf("Edges l%"CVIndexScan"\n",l);
		CVSize concentricCount = CVConcentricCountVerticesAtLevel(l, mergedInformation->concentricStructure);
		CVIndex concentricIndex;
		for (concentricIndex=0; concentricIndex<concentricCount; concentricIndex++) {
			CVSize inMergedCount = 0;
			CVIndex* groupsIndices = CVConcentricMergedConcentricIndices(concentricIndex, l, &inMergedCount, mergedInformation);
			if(inMergedCount){
				printf("\t - [");
				CVIndex mergedConcentricIndex;
				for (mergedConcentricIndex=0; mergedConcentricIndex<inMergedCount; mergedConcentricIndex++) {
					if (mergedConcentricIndex) {
						printf(" ");
					}
					printf("%"CVIndexScan"", groupsIndices[mergedConcentricIndex]);
					
				}
				printf("]\n");
				CVSize edgesCount = 0;
				CVIndex* edges = CVConcentricMergedEdgesConcentricIndices(concentricIndex, l, &edgesCount, mergedInformation);
				//printf(" (N=%"CVSizeScan" l=%"CVSizeScan" lCount=%"CVSizeScan")\n",edgesCount,l,mergedInformation->levelsCount);
				
				CVIndex edgeIndex;
				for (edgeIndex=0; edgeIndex<edgesCount; edgeIndex++) {
					printf("\t\t-> ");
					CVSize edgeInMergeCount = 0;
					CVIndex* edgeGroupIndices = CVConcentricMergedConcentricIndices(edges[edgeIndex], l+1, &edgeInMergeCount, mergedInformation);
					CVIndex mergedEdgeConcentricIndex;
					if(edgeInMergeCount){
						printf("[");
						for (mergedEdgeConcentricIndex=0; mergedEdgeConcentricIndex<edgeInMergeCount; mergedEdgeConcentricIndex++) {
							if (mergedEdgeConcentricIndex) {
								printf(" ");
							}
							printf("%"CVIndexScan"", edgeGroupIndices[mergedEdgeConcentricIndex]);
							
						}
						printf("]\n");
					}
				}
			}
		}
	}
}



void CVConcentricMergedInformationUpdateWithConcentricStructure(CVConcentricMergedInformation* mergedInformation, CVSize level,CVBool mergeLastLevel, const CVConcentricStructure* concentricStructure){
	CVConcentricMergedInformation* info = mergedInformation;
	info->concentricStructure = concentricStructure;
	
	if(level<concentricStructure->levelsCount){
		info->levelsCount=level+1;
	}else{
		info->levelsCount=concentricStructure->levelsCount;
	}
	
	const CVNetwork* network = concentricStructure->network;
	
	CVUIntegerArraySetCount(concentricStructure->totalVerticesCount, &(info->mergedTranslate));
	
	CVUIntegerArray* mergedVertices = &(info->mergedVertices);
	CVUIntegerArray* mergedIndices = &(info->mergedIndices);
	//CVUIntegerArray* levelsMergedIndices = &(info->levelsMergedIndices);
	memset(info->mergedTranslate.data,0, concentricStructure->totalVerticesCount*sizeof(CVUInteger));
	
	/*CVIndex ii;
	 
	 for(ii=0;ii<concentricStructure->totalVerticesCount;ii++){
	 info->mergedTranslate.data[ii] = -1;
	 }
	 printf("\nlevel=%"CVSizeScan" info->levelsCount=%"CVSizeScan"\n",level,info->levelsCount);
	 printf("\n\n\n|||||||\n");
	 for (l=0; l<info->levelsCount+1;l++){
	 CVSize concentricCount = 0;
	 const CVIndex* concentricVertices = CVConcentricVerticesAtLevel(l, &concentricCount, concentricStructure);
	 CVIndex concentricIndex;
	 for (concentricIndex=0; concentricIndex<concentricCount; concentricIndex++) {
	 info->mergedTranslate.data[concentricStructure->levelsIndices[l]+concentricIndex] = 9999999999;
	 printf(" %"CVSizeScan,concentricStructure->levelsIndices[l]+concentricIndex);
	 }
	 }
	 printf("\n||||||||\n");
	 */
	mergedVertices->count = 0;
	mergedIndices->count = 0;
	//levelsMergedIndices->count = 0;
	
	//CVUIntegerArrayAdd(0, levelsMergedIndices);
	CVUIntegerStackPush(0, mergedIndices);
	CVIndex l;
	for (l=0; l<info->levelsCount;l++){
		CVSize concentricCount = 0;
		const CVIndex* concentricVertices = CVConcentricVerticesAtLevel(l, &concentricCount, concentricStructure);
		
		CVGrowBitArray(concentricCount, &(info->bitArrayBufferCount), &(info->bitArrayBufferCapacity), &(info->bitArrayBuffer));
		
		CVUIntegerStack* toVisit =&(info->arrayBuffer);
		CVBitArray visited = info->bitArrayBuffer;
		CVBitArrayClearAll(info->bitArrayBuffer, concentricCount);
		toVisit->count=0;
		
		//printf("Level: %"CVIndexScan"\n",l);
		CVIndex concentricIndex;
		for (concentricIndex=0; concentricIndex<concentricCount; concentricIndex++) {
			CVSize mergedCount=0;
			if(!CVBitArrayTest(visited, concentricIndex)){
				info->mergedTranslate.data[concentricStructure->levelsIndices[l]+concentricIndex] = concentricIndex;
				CVUIntegerStackPush(concentricIndex, mergedVertices);
				mergedCount++;
				CVUIntegerArrayAdd(concentricIndex, toVisit);
				CVBitArraySet(visited, concentricIndex);
				//printf("\t -%"CVIndexScan"(%"CVIndexScan"): ",concentricVertices[concentricIndex],concentricIndex);
				while(!CVUIntegerStackIsEmpty(toVisit)){
					CVIndex currentConcentricIndex =CVUIntegerStackPop(toVisit);
					CVIndex currentVertex = concentricVertices[currentConcentricIndex];
					CVSize  currentDistance = l;
					
					CVSize vertexEdgesCount = network->vertexNumOfEdges[currentVertex];
					CVIndex* vertexEdgesList = network->vertexEdgesLists[currentVertex];
					
					CVIndex ni;
					for(ni=0;ni<vertexEdgesCount;ni++){
						CVIndex neighborVertex = vertexEdgesList[ni];
						if(concentricStructure->distances[neighborVertex]==currentDistance&&(l!=level || mergeLastLevel)){// Last level->do not merge!
							CVIndex concentricNeighborIndex = concentricStructure->concentricIndex[neighborVertex];
							if(!CVBitArrayTest(visited, concentricNeighborIndex)){
								CVUIntegerArrayAdd(concentricNeighborIndex, toVisit);
								CVBitArraySet(visited, concentricNeighborIndex);
								//printf("%"CVIndexScan"(%"CVIndexScan") ",neighborVertex,concentricNeighborIndex);
								CVIndex mergedDataIndex = concentricStructure->levelsIndices[currentDistance]+concentricNeighborIndex;
								info->mergedTranslate.data[mergedDataIndex] = concentricIndex;
								CVUIntegerStackPush(concentricNeighborIndex, mergedVertices);
								mergedCount++;
							}
						}
					}
				}
				//printf("\n");
			}
			CVUIntegerStackPush(CVUIntegerStackTop(mergedIndices)+mergedCount, mergedIndices);
		}
		//CVUIntegerStackPush(CVUIntegerStackTop(mergedIndices)+concentricCount, mergedIndices);
		//CVIndex i;
		
		/*printf("\tt = [");
		 for(i=0;i<info->mergedTranslate.count;i++){
		 printf("%"CVIndexScan" ",info->mergedTranslate.data[i]);
		 }
		 printf("]\n");
		 printf("\tth = [");
		 for (i=0; i<concentricCount; i++) {
		 printf("%"CVIndexScan" ",CVConcentricMergedTranslateIndex(i, l, info));
		 }
		 printf("]\n");*/
		/*printf("\t m = [");
		 for(i=0;i<mergedVertices->count;i++){
		 printf("%"CVIndexScan" ",mergedVertices->data[i]);
		 }
		 printf("]\n");
		 printf("\tmi = [");
		 for(i=0;i<mergedIndices->count;i++){
		 printf("%"CVIndexScan" ",mergedIndices->data[i]);
		 }
		 printf("]\n");
		 printf("\t h = [");
		 for (i=0; i<concentricStructure->levelsCount+1;i++){
		 printf("%"CVSizeScan" ",concentricStructure->levelsIndices[i]);
		 }
		 printf("]\n");*/
	}
	
	//FIX Edges
	
	CVUIntegerArray* edges = &(info->edges);
	CVUIntegerArray* edgesIndices = &(info->edgesIndices);
	//CVUIntegerArray* levelsEdgesIndices = &(info->levelsEdgesIndices);
	
	edges->count = 0;
	edgesIndices->count = 0;
	//	levelsEdgesIndices->count = 0;
	CVUIntegerStackPush(0, edgesIndices);
	for (l=0; l<info->levelsCount-1;l++){
		CVSize concentricCount = 0;
		CVSize nextLevelConcentricCount = 0;
		
		const CVIndex* concentricVertices = CVConcentricVerticesAtLevel(l, &concentricCount, concentricStructure);
		nextLevelConcentricCount = CVConcentricCountVerticesAtLevel(l+1, concentricStructure);
		
		CVGrowBitArray(nextLevelConcentricCount, &(info->bitArrayBufferCount), &(info->bitArrayBufferCapacity), &(info->bitArrayBuffer));
		
		CVBitArray nextLevelVisited = info->bitArrayBuffer;
		
		CVIndex concentricIndex;
		for (concentricIndex=0; concentricIndex<concentricCount; concentricIndex++) {
			CVBitArrayClearAll(nextLevelVisited, nextLevelConcentricCount);
			
			CVSize mergedCount = 0;
			const CVIndex* mergedVertices = CVConcentricMergedConcentricIndices(concentricIndex, l, &mergedCount, info);
			
			CVSize edgesCount = 0;
			
			CVIndex mergedIndex;
			for (mergedIndex=0; mergedIndex<mergedCount; mergedIndex++) {
				CVIndex mergedConcentricIndex = mergedVertices[mergedIndex];
				
				CVIndex vertexIndex = concentricVertices[mergedConcentricIndex];
				CVSize vertexEdgesCount = network->vertexNumOfEdges[vertexIndex];
				CVIndex* vertexEdgesList = network->vertexEdgesLists[vertexIndex];
				
				CVIndex ni;
				for(ni=0;ni<vertexEdgesCount;ni++){
					CVIndex neighborVertex = vertexEdgesList[ni];
					if(concentricStructure->distances[neighborVertex]==l+1){
						CVIndex neighborConcentricIndex = concentricStructure->concentricIndex[neighborVertex];
						//printf("\n%"CVSizeScan"/%"CVSizeScan,neighborConcentricIndex,nextLevelConcentricCount);
						
						CVIndex mergedConcentricIndex = CVConcentricMergedTranslateIndex(neighborConcentricIndex, l+1, info);
						
						/*CVSize tamanho = CVBitArraySlot(mergedConcentricIndex);
						 if(tamanho>network->verticesCount){
						 printf("\nLevelsCount:%"CVSizeScan" level:%"CVSizeScan" neighConcentricIndex:%"CVSizeScan" concentricMax:%"CVSizeScan" mergedConcentricIndex:%"CVSizeScan" tamanho:%"CVSizeScan"\n",
						 concentricStructure->levelsCount,l,neighborConcentricIndex,CVConcentricCountVerticesAtLevel(l+1, concentricStructure),mergedConcentricIndex,tamanho);
						 
						 printf("mergedInformation->mergedTranslate.count=%"CVSizeScan"\nconcentricStructure->levelsIndices[level+1]=%"CVSizeScan"\n",mergedInformation->mergedTranslate.count,concentricStructure->levelsIndices[level+1]);
						 CVIndex j;
						 for (j=0; j<mergedInformation->mergedTranslate.count; j++) {
						 if(j){
						 printf(" ");
						 }
						 printf("%"CVSizeScan,mergedInformation->mergedTranslate.data[j]);
						 }
						 printf("\n");
						 CVConcentricStructurePrint(concentricStructure);
						 CVConcentricMergedInformationPrint(mergedInformation);
						 CVBitArrayTest(nextLevelVisited, mergedConcentricIndex)?printf("0\n"):printf("1\n");
						 }*/
						
						
						//NEWVERSION
						//if(!CVBitArrayTest(nextLevelVisited, mergedConcentricIndex)){
						//	CVBitArraySet(nextLevelVisited, mergedConcentricIndex);
							CVUIntegerStackPush(mergedConcentricIndex, edges);
							edgesCount++;
						//}
					}
				}
			}
			CVUIntegerStackPush(CVUIntegerStackTop(edgesIndices)+edgesCount, edgesIndices);
		}
	}
}

/*
 edges					[0,1,2,0,1,2]
 l: 0
 0-> [0,1,2]
 l: 1
 0-> []
 1-> [0,1]
 2-> [2]
 
 edgesIndices;			[[0,3],[3,3,5,6]]
 levelsEdgesIndices;	[0,2,6]
 
 mergedTranslate
 */



void _CV_ConcentricMergedPropagateProbabilities(CVFloatArray* mergedProbabilities, CVUIntegerArray* mergedPathsCount, CVSize level, CVSize* deadEndsCount, const CVConcentricMergedInformation* info){
	const CVConcentricStructure* structure = info->concentricStructure;
	
	CVSize concentricCount = CVConcentricCountVerticesAtLevel(level, structure);
	
	CVIndex concentricIndex;
	for (concentricIndex=0; concentricIndex<concentricCount; concentricIndex++) {
		CVSize inMergedCount = 0;
		CVConcentricMergedConcentricIndices(concentricIndex, level, &inMergedCount, info);
		if(inMergedCount){
			CVSize edgesCount = 0;
			CVIndex* edges = CVConcentricMergedEdgesConcentricIndices(concentricIndex, level, &edgesCount, info);
			//printf(" (N=%"CVSizeScan" l=%"CVSizeScan" lCount=%"CVSizeScan")\n",edgesCount,l,mergedInformation->levelsCount);
			
			CVIndex edgeIndex;
			CVSize betweenLevelConnections = edgesCount;
			CVFloat propagateFactor = mergedProbabilities->data[structure->levelsIndices[level]+concentricIndex]/betweenLevelConnections;
			CVSize pathsFactor = mergedPathsCount->data[structure->levelsIndices[level]+concentricIndex];
			for (edgeIndex=0; edgeIndex<edgesCount; edgeIndex++) {
				mergedProbabilities->data[structure->levelsIndices[level+1]+edges[edgeIndex]]+= propagateFactor;
				mergedPathsCount->data[structure->levelsIndices[level+1]+edges[edgeIndex]]+= pathsFactor;
			}
			if(betweenLevelConnections==0){
				(*deadEndsCount)++;
			}
		}
	}
}

void CVConcentricMergedGetProbabilities(CVFloatArray* mergedProbabilities,CVUIntegerArray* mergedPathsCount, CVSize* deadEndCounts, CVIndex level, const CVConcentricMergedInformation* info){
	const CVConcentricStructure* structure = info->concentricStructure;
	CVFloatArraySetCount(structure->totalVerticesCount, mergedProbabilities);
	CVUIntegerArraySetCount(structure->totalVerticesCount, mergedPathsCount);
	
	memset(mergedProbabilities->data,0, sizeof(CVFloat)*structure->totalVerticesCount);
	memset(mergedPathsCount->data,0, sizeof(CVSize)*structure->totalVerticesCount);
	
	mergedProbabilities->data[0] = 1.0f;
	mergedPathsCount->data[0] = 1;
	
	CVIndex l;
	CVSize levelCount = CVMIN(info->levelsCount, level+1);
	if(deadEndCounts){
		deadEndCounts[0]=0;
	}
	for (l=0; l<levelCount-1; l++) {
		CVSize deadEndCount = 0;
		_CV_ConcentricMergedPropagateProbabilities(mergedProbabilities,mergedPathsCount, l, &deadEndCount, info);
		if(deadEndCounts){
			deadEndCounts[l+1] = deadEndCounts[l]+deadEndCount;
		}
	}
	
}

void CVConcentricMergedInformationDestroy(CVConcentricMergedInformation* information){
	CVUIntegerArrayDestroy(&information->edges);
	CVUIntegerArrayDestroy(&information->edgesIndices);
	CVUIntegerArrayDestroy(&information->mergedVertices);
	CVUIntegerArrayDestroy(&information->mergedIndices);
	
	CVUIntegerArrayDestroy(&information->mergedTranslate);
	
	CVBitArrayDestroy(information->bitArrayBuffer);
	
	CVUIntegerArrayDestroy(&information->arrayBuffer);
	CVUIntegerArrayDestroy(&information->arrayBuffer2);
	
	free(information);
}



CVNetwork* CVNewNetworkFromConcentricStructure(CVConcentricStructure* concentricStructure, CVSize maximumLevel){
	CVSize verticesCount = 0;
	CVSize levelsCount = CVMIN(concentricStructure->levelsCount,maximumLevel+1);
	const CVNetwork* originalNetwork = concentricStructure->network;
	CVIndex level;
	for (level=0; level<levelsCount; level++) {
		verticesCount += CVConcentricCountVerticesAtLevel(level, concentricStructure);
	}
	CVNetwork* theNetwork = CVNewNetwork(verticesCount, CVFalse, CVFalse);
	CVFloat* levelsArray = calloc(verticesCount, sizeof(CVFloat));
	for (level=0; level<levelsCount; level++){
		CVSize concentricCount = 0;
		const CVIndex* concentricVertices = CVConcentricVerticesAtLevel(level, &concentricCount, concentricStructure);
		CVIndex concentricIndex;
		CVIndex levelStartIndex = concentricStructure->levelsIndices[level];
		for (concentricIndex=0; concentricIndex<concentricCount; concentricIndex++) {
			CVIndex fromNewIndex = levelStartIndex+concentricIndex;
			levelsArray[fromNewIndex] = level;
			CVIndex fromIndex = concentricVertices[concentricIndex];
			
			CVSize vertexEdgesCount = originalNetwork->vertexNumOfEdges[fromIndex];
			CVIndex* vertexEdgesList = originalNetwork->vertexEdgesLists[fromIndex];
			
			CVIndex ni;
			for(ni=0;ni<vertexEdgesCount;ni++){
				CVIndex neighborVertex = vertexEdgesList[ni];
				CVSize neighborVertexLevel = concentricStructure->distances[neighborVertex];
				if(neighborVertexLevel==level){ // Same Level
					CVIndex concentricToIndex = concentricStructure->concentricIndex[neighborVertex];
					CVIndex toNewIndex = levelStartIndex+concentricToIndex;
					CVNetworkAddNewEdge(theNetwork, fromNewIndex, toNewIndex, 1.0f);
				}else if(neighborVertexLevel==level+1 && (level+1)<levelsCount){
					CVIndex concentricToIndex = concentricStructure->concentricIndex[neighborVertex];
					CVIndex toNewIndex = concentricStructure->levelsIndices[level+1]+concentricToIndex;
					CVNetworkAddNewEdge(theNetwork, fromNewIndex, toNewIndex, 1.0f);
				}
			}
		}
	}
	CVNetworkAppendProperty(theNetwork, "level", CVNumberPropertyType, levelsArray);
	free(levelsArray);
	return theNetwork;
}










