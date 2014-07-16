//
//  CVNetworkSymmetry.h
//  CVNetwork
//
//  Created by Filipi Nascimento Silva on 05/05/14.
//  Copyright (c) 2014 Filipi Nascimento Silva. All rights reserved.
//

#ifndef __CVNetwork__CVNetworkSymmetry__
#define __CVNetwork__CVNetworkSymmetry__

#include "CVNetwork.h"
#include "CVConcentricStructure.h"
#include "CVBasicArrays.h"

typedef struct{
	const CVNetwork* network;
	CVIndex referenceIndex;
	CVSize maximumLevel;
	CVBool saveProbabilities;
	CVBool saveExtraInformation;
	CVBool mergeLastLevel;
	
	CVConcentricStructure* concentricStructure;
	CVConcentricMergedInformation* mergedInformation;
	
	CVFloat* ___backboneProbabilitiesBuffer;
	CVFloatArray ___mergedProbabilitiesBuffer;
	
	CVSize* ___backbonePathsCountBuffer;
	CVUIntegerArray ___mergedPathsCountBuffer;
	
	CVSize* ___backboneDeadEndCounts;
	CVSize* ___mergedDeadEndCounts;
} CVSymmetryInputParameters;


typedef struct{
	CVSize levelsCount;
	CVFloatArray* backboneProbabilities;
	CVFloatArray* mergedProbabilities;
	CVUIntegerArray* lastLevelIndices;
	
	CVFloat* backboneAccessibility;
	CVFloat* mergedAccessibility;
	
	CVFloat* normalizedBackboneAccessibility;
	CVFloat* normalizedMergedAccessibility;
	
	CVFloat* normalizedBackbonePathAccessibility;
	CVFloat* normalizedMergedPathAccessibility;
	
	CVSize* accessedVertices;
	CVSize* mergedAccessedVertices;
	
	CVSize* backboneDeadEnds;
	CVSize* mergedDeadEnds;
} CVSymmetryOutputParameters;


CV_INLINE CVSymmetryInputParameters* CVNewSymmetryInput(const CVNetwork* network, CVIndex referenceIndex, CVSize maximumLevel, CVBool saveExtraInformation, CVBool saveProbabilities, CVBool mergeLastLevel){
	CVSymmetryInputParameters* inputParameters = calloc(1, sizeof(CVSymmetryInputParameters));
	
	inputParameters->network = network;
	inputParameters->maximumLevel = maximumLevel;
	inputParameters->referenceIndex=referenceIndex;
	inputParameters->saveProbabilities=saveProbabilities;
	inputParameters->saveExtraInformation = saveExtraInformation;
	inputParameters->mergeLastLevel = mergeLastLevel;
	
	inputParameters->___backboneProbabilitiesBuffer = calloc(network->verticesCount, sizeof(CVFloat));
	inputParameters->___backbonePathsCountBuffer = calloc(network->verticesCount, sizeof(CVSize));
	
	inputParameters->___backboneDeadEndCounts = calloc(maximumLevel+1, sizeof(CVSize));
	inputParameters->___mergedDeadEndCounts = calloc(maximumLevel+1, sizeof(CVSize));
	
	CVFloatArrayInitWithCapacity(20, &inputParameters->___mergedProbabilitiesBuffer);
	CVUIntegerArrayInitWithCapacity(20, &inputParameters->___mergedPathsCountBuffer);
	
	inputParameters->concentricStructure = CVNewConcentricStructureForNetwork(network, CVTrue);
	inputParameters->mergedInformation = CVNewConcentricMergedInformation();
	
	CVConcentricStructureSetReferenceVertex(inputParameters->referenceIndex, inputParameters->maximumLevel, inputParameters->concentricStructure);
	
	return inputParameters;
}

CV_INLINE void CVSymmetryInputDestroy(CVSymmetryInputParameters* inputParameters){
	free(inputParameters->___backboneProbabilitiesBuffer);
	free(inputParameters->___backbonePathsCountBuffer);
	free(inputParameters->___backboneDeadEndCounts);
	free(inputParameters->___mergedDeadEndCounts);
	CVFloatArrayDestroy(&inputParameters->___mergedProbabilitiesBuffer);
	CVUIntegerArrayDestroy(&inputParameters->___mergedPathsCountBuffer);
	CVConcentricStructureDestroy(inputParameters->concentricStructure);
	CVConcentricMergedInformationDestroy(inputParameters->mergedInformation);
	
	free(inputParameters);
}

CV_INLINE CVSymmetryOutputParameters* CVNewSymmetryOutput(){
	return calloc(1, sizeof(CVSymmetryOutputParameters));
}


CV_INLINE void CVSymmetryOutputInitFromInput(CVSymmetryOutputParameters* outputParameters, const CVSymmetryInputParameters* inputParameters){
	outputParameters->levelsCount = inputParameters->maximumLevel+1;
	if(inputParameters->saveProbabilities){
		outputParameters->backboneProbabilities = calloc(outputParameters->levelsCount, sizeof(CVFloatArray));
		outputParameters->mergedProbabilities = calloc(outputParameters->levelsCount, sizeof(CVFloatArray));
		outputParameters->lastLevelIndices = calloc(outputParameters->levelsCount, sizeof(CVUIntegerArray));
		CVIndex l;
		for(l=0; l<outputParameters->levelsCount; l++) {
			CVFloatArrayInitWithCapacity(5, outputParameters->backboneProbabilities+l);
			CVFloatArrayInitWithCapacity(5, outputParameters->mergedProbabilities+l);
			CVUIntegerArrayInitWithCapacity(5, outputParameters->lastLevelIndices+l);
		}
	}
	
	if(inputParameters->saveExtraInformation){
		outputParameters->accessedVertices = calloc(outputParameters->levelsCount, sizeof(CVSize));
		outputParameters->mergedAccessedVertices = calloc(outputParameters->levelsCount, sizeof(CVSize));
		outputParameters->backboneAccessibility = calloc(outputParameters->levelsCount, sizeof(CVFloat));
		outputParameters->mergedAccessibility = calloc(outputParameters->levelsCount, sizeof(CVFloat));
		
		outputParameters->backboneDeadEnds = calloc(outputParameters->levelsCount, sizeof(CVSize));
		outputParameters->mergedDeadEnds = calloc(outputParameters->levelsCount, sizeof(CVSize));
	}
	
	outputParameters->normalizedBackboneAccessibility = calloc(outputParameters->levelsCount, sizeof(CVFloat));
	outputParameters->normalizedMergedAccessibility = calloc(outputParameters->levelsCount, sizeof(CVFloat));
	
	outputParameters->normalizedBackbonePathAccessibility = calloc(outputParameters->levelsCount, sizeof(CVFloat));
	outputParameters->normalizedMergedPathAccessibility = calloc(outputParameters->levelsCount, sizeof(CVFloat));
}

CV_INLINE void CVSymmetryOutputDestroy(CVSymmetryOutputParameters* outputParameters){
	if(outputParameters->backboneProbabilities){
		CVIndex l;
		for(l=0; l<outputParameters->levelsCount; l++) {
			CVFloatArrayDestroy(outputParameters->backboneProbabilities+l);
		}
		free(outputParameters->backboneProbabilities);
	}
	
	if(outputParameters->mergedProbabilities){
		CVIndex l;
		for(l=0; l<outputParameters->levelsCount; l++) {
			CVFloatArrayDestroy(outputParameters->mergedProbabilities+l);
		}
		free(outputParameters->mergedProbabilities);
	}
	
	if(outputParameters->lastLevelIndices){
		CVIndex l;
		for(l=0; l<outputParameters->levelsCount; l++) {
			CVUIntegerArrayDestroy(outputParameters->lastLevelIndices+l);
		}
		free(outputParameters->lastLevelIndices);
	}
	
	free(outputParameters->accessedVertices);
	free(outputParameters->mergedAccessedVertices);
	
	free(outputParameters->backboneAccessibility);
	free(outputParameters->mergedAccessibility);
	
	free(outputParameters->normalizedBackboneAccessibility);
	free(outputParameters->normalizedMergedAccessibility);
	
	free(outputParameters->normalizedBackbonePathAccessibility);
	free(outputParameters->normalizedMergedPathAccessibility);
	
	free(outputParameters->backboneDeadEnds);
	free(outputParameters->mergedDeadEnds);
	free(outputParameters);
}

CV_INLINE void CVSymmetryInputChangeReferenceVertex(CVIndex newReferenceVertex, CVSymmetryInputParameters* symmetryInputParameters){
	if(symmetryInputParameters->referenceIndex!=newReferenceVertex){
		symmetryInputParameters->referenceIndex = newReferenceVertex;
		CVConcentricStructureSetReferenceVertex(newReferenceVertex, symmetryInputParameters->maximumLevel, symmetryInputParameters->concentricStructure);
	}
}


CVBool CVNetworkCalculateSymmetryOfVertex(CVSymmetryInputParameters* inputParameters, CVSymmetryOutputParameters* outputParameters);


CVBool CVNetworkCalculateSymmetry(const CVNetwork* network, CVSize maximumLevel, CVBool saveExtraInformation, CVBool saveProbabilities,CVBool mergeLastLevel, CVSymmetryOutputParameters** outputArray, CVOperationControl* operationControl);

#endif /* defined(__CVNetwork__CVNetworkSymmetry__) */
