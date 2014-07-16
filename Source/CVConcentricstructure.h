//
//  CVConcentricStructure.h
//  CVNetwork
//
//  Created by Filipi Nascimento Silva on 26/04/14.
//  Copyright (c) 2014 Filipi Nascimento Silva. All rights reserved.
//

#ifndef CVNetwork_CVConcentricStructure_h
#define CVNetwork_CVConcentricStructure_h

#include "CVNetwork.h"
#include "CVBasicArrays.h"
#include "CVSimpleQueue.h"




typedef struct{
	CVSize* distances;
	CVIndex* concentricIndex;
	CVIndex referenceVertex;
	
	CVIndex* vertices;
	CVSize totalVerticesCount;
	CVSize totalVerticesCapacity;
	
	CVIndex* levelsIndices;
	CVSize levelsCapacity;
	CVSize levelsCount;
	const CVNetwork* network;
	CVBitArray _visitedNodes;
	
} CVConcentricStructure;

CVConcentricStructure* CVNewConcentricStructureForNetwork(const CVNetwork* network, CVBool trackNodes);

void _CV_ConcentricStructureAddVertexAtDistance(CVIndex vertexIndex, CVSize distance, CVConcentricStructure* structure);

void CVConcentricStructureSetReferenceVertex(CVIndex referenceVertex, CVSize maximumConcentricLevel, CVConcentricStructure* structure);

CV_INLINE const CVIndex* CVConcentricVerticesAtLevel(CVSize level, CVSize* count, const CVConcentricStructure* structure){
	if(structure->vertices && level < structure->levelsCount){
		*count = structure->levelsIndices[level+1]-structure->levelsIndices[level];
		return structure->vertices + structure->levelsIndices[level];
	}else{
		*count = 0;
		return NULL;
	}
}

CV_INLINE CVSize CVConcentricCountVerticesAtLevel(CVSize level, const CVConcentricStructure* structure){
	if(level < structure->levelsCount){
		return structure->levelsIndices[level+1]-structure->levelsIndices[level];
	}else{
		return 0;
	}
}

void CVConcentricStructurePrint(const CVConcentricStructure* structure);

void _CV_ConcentricBackbonePropagateProbabilities(CVFloat* probabilities, CVSize* pahtsCount, CVSize level,CVSize* deadEndsCount, const CVConcentricStructure* structure);

void CVConcentricBackboneGetProbabilities(CVFloat* probabilities, CVSize* pahtsCount, CVSize* deadEndCounts, CVSize maximumConcentricLevel, const CVConcentricStructure* structure);

void CVConcentricStructureDestroy(CVConcentricStructure* structure);

typedef struct{
	CVUIntegerArray edges;
	CVUIntegerArray edgesIndices;
	//CVUIntegerArray levelsEdgesIndices;
	
	
	CVUIntegerArray mergedVertices;
	CVUIntegerArray mergedIndices;
	//CVUIntegerArray levelsMergedIndices;
	
	CVUIntegerArray mergedTranslate;
	
	CVBitArray bitArrayBuffer;
	CVSize bitArrayBufferCapacity;
	CVSize bitArrayBufferCount;
	
	CVUIntegerArray arrayBuffer;
	CVUIntegerArray arrayBuffer2;
	
	
	const CVConcentricStructure* concentricStructure;
	CVSize levelsCount;
} CVConcentricMergedInformation;



CV_INLINE CVConcentricMergedInformation* CVNewConcentricMergedInformation(){
	CVConcentricMergedInformation* mergedInformation = calloc(1, sizeof(CVConcentricMergedInformation));
	return mergedInformation;
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

CV_INLINE CVIndex* _CV_ConcentricNeighVerticesInNextLevel(CVIndex referenceConcentricIndex, CVIndex level, CVSize neighVerticesCount, CVConcentricMergedInformation* info){
	//CVConcentricMergedInformation* mergedInformation = calloc(1, sizeof(CVConcentricMergedInformation));
	return 0;
}

CV_INLINE CVIndex CVConcentricMergedTranslateIndex(CVIndex concentricIndex, CVSize level, const CVConcentricMergedInformation* info){
	const CVConcentricStructure* structure = info->concentricStructure;
	return info->mergedTranslate.data[structure->levelsIndices[level]+concentricIndex];
}

CV_INLINE CVIndex* CVConcentricMergedConcentricIndices(CVIndex concentricIndex, CVSize level, CVSize* mergedCount, const CVConcentricMergedInformation* info){
	const CVConcentricStructure* structure = info->concentricStructure;
	if(level < info->levelsCount && concentricIndex < structure->levelsIndices[level+1]-structure->levelsIndices[level]){
		CVIndex* mergedIndices = info->mergedIndices.data + structure->levelsIndices[level];
		*mergedCount = mergedIndices[concentricIndex+1]-mergedIndices[concentricIndex];
		return info->mergedVertices.data + mergedIndices[concentricIndex];
	}else{
		*mergedCount = 0;
		return NULL;
	}
}


CV_INLINE CVIndex* CVConcentricMergedEdgesConcentricIndices(CVIndex concentricIndex, CVSize level, CVSize* edgesCount, const CVConcentricMergedInformation* info){
	const CVConcentricStructure* structure = info->concentricStructure;
	if(level < info->levelsCount-1 && concentricIndex<structure->levelsIndices[level+1]-structure->levelsIndices[level]){
		CVIndex* edgesIndices = info->edgesIndices.data + structure->levelsIndices[level];
		*edgesCount = edgesIndices[concentricIndex+1]-edgesIndices[concentricIndex];
		return info->edges.data + edgesIndices[concentricIndex];
	}else{
		*edgesCount = 0;
		return NULL;
	}
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


void CVConcentricMergedInformationPrint(CVConcentricMergedInformation* mergedInformation);


void CVConcentricMergedInformationUpdateWithConcentricStructure(CVConcentricMergedInformation* mergedInformation, CVSize level,CVBool mergeLastLevel, const CVConcentricStructure* concentricStructure);

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

void _CV_ConcentricMergedPropagateProbabilities(CVFloatArray* mergedProbabilities, CVUIntegerArray* mergedPathsCount, CVSize level, CVSize* deadEndsCount, const CVConcentricMergedInformation* info);

void CVConcentricMergedGetProbabilities(CVFloatArray* mergedProbabilities,CVUIntegerArray* mergedPathsCount, CVSize* deadEndCounts, CVIndex level, const CVConcentricMergedInformation* info);
void CVConcentricMergedInformationDestroy(CVConcentricMergedInformation* information);

CVNetwork* CVNewNetworkFromConcentricStructure(CVConcentricStructure* concentricStructure, CVSize level);

#endif
