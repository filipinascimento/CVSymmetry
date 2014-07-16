//
//  CVNetworkSymmetry.cpp
//  CVNetwork
//
//  Created by Filipi Nascimento Silva on 05/05/14.
//  Copyright (c) 2014 Filipi Nascimento Silva. All rights reserved.
//

#include "CVNetworkSymmetry.h"
#include "CVConcentricStructure.h"
#include "CVNetworkSymmetry.h"



CVBool CVNetworkCalculateSymmetryOfVertex(CVSymmetryInputParameters* inputParameters, CVSymmetryOutputParameters* outputParameters){
	if(inputParameters->maximumLevel<2){
		return CVFalse;
	}
	CVSymmetryOutputInitFromInput(outputParameters, inputParameters);
	
	CVBool saveExtraInformation = inputParameters->saveExtraInformation;
	CVBool saveProbabilities = inputParameters->saveProbabilities;
	
	CVConcentricStructure* concentricStructure = inputParameters->concentricStructure;
	CVConcentricMergedInformation* mergedInformation = mergedInformation = inputParameters->mergedInformation;
	
	CVFloat* backboneProbabilitiesBuffer = inputParameters->___backboneProbabilitiesBuffer;
	CVSize* backbonePathsCountBuffer = inputParameters->___backbonePathsCountBuffer;
	
	CVSize* backboneDeadEndCounts = inputParameters->___backboneDeadEndCounts;
	
	CVFloatArray* mergedProbabilitiesBuffer = &(inputParameters->___mergedProbabilitiesBuffer);
	CVUIntegerArray* mergedPathsCountBuffer = &(inputParameters->___mergedPathsCountBuffer);
	
	CVSize* mergedDeadEndCounts = inputParameters->___mergedDeadEndCounts;
	
	CVConcentricBackboneGetProbabilities(inputParameters->___backboneProbabilitiesBuffer,
										 inputParameters->___backbonePathsCountBuffer,
										 inputParameters->___backboneDeadEndCounts,
										 inputParameters->maximumLevel,
										 concentricStructure);
	
	CVIndex level;
	for (level=2; level<inputParameters->maximumLevel+1;level++){
		CVSize concentricCount = 0;
		const CVIndex* concentricVertices = CVConcentricVerticesAtLevel(level, &concentricCount, concentricStructure);
		CVIndex concentricIndex;
		
		double backboneSumProbability = 0.0;
		double backboneEntropy = 0.0;
		
		CVSize backboneSumPathsCount = 0.0;
		double backbonePathEntropy = 0.0;
		
		for (concentricIndex=0; concentricIndex<concentricCount; concentricIndex++) {
			double probability = backboneProbabilitiesBuffer[concentricVertices[concentricIndex]];
			CVSize pathsCount = backbonePathsCountBuffer[concentricVertices[concentricIndex]];
			
			backboneSumProbability += probability;
			backboneSumPathsCount += pathsCount;
		}
		
		for (concentricIndex=0; concentricIndex<concentricCount; concentricIndex++) {
			CVIndex vertexIndex = concentricVertices[concentricIndex];
			
			double probability = backboneProbabilitiesBuffer[vertexIndex];
			CVSize pathsCount = backbonePathsCountBuffer[concentricVertices[concentricIndex]];
			
			if(probability>0.0&&pathsCount>0){
				double normProbability = probability/backboneSumProbability;
				backboneEntropy -= normProbability*log(normProbability);
				
				double normPathsCount = pathsCount/(double)backboneSumPathsCount;
				backbonePathEntropy -= normPathsCount*log(normPathsCount);
				
				if(saveProbabilities){
					CVFloatArrayAdd(normProbability, outputParameters->backboneProbabilities+level);
					CVUIntegerArrayAdd(vertexIndex, outputParameters->lastLevelIndices+level);
				}
			}
			//backboneProbabilitiesBuffer[concentricVertices[concentricIndex]]=0.0f;
		}
		
		if((concentricCount)>0){
			outputParameters->normalizedBackboneAccessibility[level] = exp(backboneEntropy)/(concentricCount+backboneDeadEndCounts[level]);
			outputParameters->normalizedBackbonePathAccessibility[level] = exp(backbonePathEntropy)/(concentricCount+backboneDeadEndCounts[level]);
		}else{
			outputParameters->normalizedBackboneAccessibility[level] = 0.0f;
			outputParameters->normalizedBackbonePathAccessibility[level] = 0.0f;
		}

			
		if(saveExtraInformation){
			if((concentricCount)>0){
				outputParameters->backboneAccessibility[level] = exp(backboneEntropy);
			}
			outputParameters->backboneDeadEnds[level] = backboneDeadEndCounts[level];
			outputParameters->accessedVertices[level] = concentricCount;
		}
		
		CVConcentricMergedInformationUpdateWithConcentricStructure(mergedInformation, level, inputParameters->mergeLastLevel, concentricStructure);
		CVConcentricMergedGetProbabilities(mergedProbabilitiesBuffer,mergedPathsCountBuffer, mergedDeadEndCounts, level, mergedInformation);
		
		double mergedSumProbability = 0.0;
		double mergedEntropy = 0.0;
		CVSize mergedAccessedCount = 0;
		
		CVSize mergedSumPathsCount = 0.0;
		double mergedPathEntropy = 0.0;
		
		for (concentricIndex=0; concentricIndex<concentricCount; concentricIndex++) {
			double probability = mergedProbabilitiesBuffer->data[concentricStructure->levelsIndices[level]+concentricIndex];
			CVSize pathsCount = mergedPathsCountBuffer->data[concentricStructure->levelsIndices[level]+concentricIndex];
			
			if(probability>0.0&&pathsCount>0.0){
				mergedSumProbability+=probability;
				mergedSumPathsCount+=pathsCount;
				mergedAccessedCount++;
			}
		}
		
		for (concentricIndex=0; concentricIndex<concentricCount; concentricIndex++) {
			//CVInteger vertexIndex = concentricVertices[concentricIndex];
			double probability = mergedProbabilitiesBuffer->data[concentricStructure->levelsIndices[level]+concentricIndex];
			CVSize pathsCount = mergedPathsCountBuffer->data[concentricStructure->levelsIndices[level]+concentricIndex];
			if(probability>0.0&&pathsCount>0){
				double normProbability = probability/mergedSumProbability;
				double normPathsCount = pathsCount/(double)mergedSumPathsCount;
				mergedEntropy-=normProbability*log(normProbability);
				mergedPathEntropy-=normPathsCount*log(normPathsCount);
				if(saveProbabilities){
					CVFloatArrayAdd(normProbability, outputParameters->mergedProbabilities+level);
				}
			}
		}
		if(mergedAccessedCount>0){
			outputParameters->normalizedMergedAccessibility[level] = exp(mergedEntropy)/(mergedAccessedCount+mergedDeadEndCounts[level]);
			outputParameters->normalizedMergedPathAccessibility[level] = exp(mergedPathEntropy)/(mergedAccessedCount+mergedDeadEndCounts[level]);
		}else{
			outputParameters->normalizedMergedAccessibility[level] = 0.0f;
			outputParameters->normalizedMergedPathAccessibility[level] = 0.0f;
		}
			
		if(saveExtraInformation){
			if((mergedAccessedCount)>0){
				outputParameters->mergedAccessibility[level] = exp(mergedEntropy);
			}
			outputParameters->mergedDeadEnds[level] = mergedDeadEndCounts[level];
			outputParameters->mergedAccessedVertices[level] = mergedAccessedCount;
		}
	}

	
	return CVTrue;
}


CVBool _CVNetworkCalculateSymmetry_parallel_implementation(const CVNetwork* network, CVSize maximumLevel, CVBool saveExtraInformation, CVBool saveProbabilities,CVBool mergeLastLevel, CVSymmetryOutputParameters** outputArray, CVOperationControl* operationControl){
	CVSize verticesCount = network->verticesCount;
	CVSize unrolledLoops = kCVDefaultParallelBlocks;
	//__block BOOL stop = NO;
	CVInteger* currentProgress = NULL;
	
	void (*updateCallback)(CVOperationControl*)  = NULL;
	void (*streamCallback)(CVOperationControl*, CVIndex index, const char* format,...)  = NULL;
	if(operationControl){
		operationControl->maxProgress = verticesCount;
		operationControl->currentProgress = 0;
		currentProgress = &(operationControl->currentProgress);
		if(operationControl->maxParallelBlocks>0){
			unrolledLoops = operationControl->maxParallelBlocks;
		}
		updateCallback = operationControl->updateCallback;
		streamCallback = operationControl->streamCallback;
	}
	if(updateCallback){
		updateCallback(operationControl);
	}
	CVSize unrolledSize = 1 + ((verticesCount - 1) / unrolledLoops);
	//printf("unrolledSize:%"CVSizeScan" UnrolledLoops%"CVSizeScan"\n",unrolledSize,unrolledLoops);
	
	
	CVParallelForStart(symmetryLoop, blockIndex, unrolledLoops){
		CVIndex referenceVertex = blockIndex*unrolledSize;
		CVSize maxIndex = CVMIN((blockIndex+1)*unrolledSize, verticesCount);
		//printf("CALCULATING: %"CVSizeScan"-%"CVSizeScan"\n",blockIndex*unrolledSize,maxIndex);
		
		if(referenceVertex<maxIndex){
			CVSymmetryInputParameters* inputParameters;
			
			inputParameters= CVNewSymmetryInput(network, referenceVertex, maximumLevel, saveExtraInformation, saveProbabilities,mergeLastLevel);
			
			
			for(referenceVertex=blockIndex*unrolledSize;referenceVertex<maxIndex;referenceVertex++){
				CVSymmetryInputChangeReferenceVertex(referenceVertex, inputParameters);
				
				outputArray[referenceVertex] = CVNewSymmetryOutput();
				
				CVNetworkCalculateSymmetryOfVertex(inputParameters, outputArray[referenceVertex]);
				
				if(currentProgress){
					CVAtomicIncrementInteger(currentProgress);
					if(updateCallback){
						CVParallelLoopCriticalRegionStart(symmetryLoop){
							updateCallback(operationControl);
						}
						CVParallelLoopCriticalRegionEnd(symmetryLoop);
					}
					if(streamCallback){
						CVParallelLoopCriticalRegionStart(symmetryLoop){
							streamCallback(operationControl,referenceVertex,"%"CVFloatScan,0.1);//Value Streamming
						}
						CVParallelLoopCriticalRegionEnd(symmetryLoop);
					}
				}
			}
			//ClearFiles
			CVSymmetryInputDestroy(inputParameters);
			
			
		}
	}CVParallelForEnd(symmetryLoop);
	
	return CVTrue;
}

CVBool _CVNetworkCalculateSymmetry_implementation(const CVNetwork* network, CVSize maximumLevel, CVBool saveExtraInformation, CVBool saveProbabilities,CVBool mergeLastLevel, CVSymmetryOutputParameters** outputArray, CVOperationControl* operationControl){
	CVSize verticesCount = network->verticesCount;
	//__block BOOL stop = NO;
	CVInteger* currentProgress = NULL;
	
	void (*updateCallback)(CVOperationControl*)  = NULL;
	void (*streamCallback)(CVOperationControl*, CVIndex index, const char* format,...)  = NULL;
	if(operationControl){
		operationControl->maxProgress = verticesCount;
		operationControl->currentProgress = 0;
		currentProgress = &(operationControl->currentProgress);
		updateCallback = operationControl->updateCallback;
		streamCallback = operationControl->streamCallback;
	}
	if(updateCallback){
		updateCallback(operationControl);
	}
	CVSymmetryInputParameters* inputParameters;
	if(verticesCount>0){
		inputParameters = CVNewSymmetryInput(network, 0, maximumLevel, saveExtraInformation, saveProbabilities,mergeLastLevel);
	}
	CVIndex referenceVertex;
	
	for(referenceVertex=0;referenceVertex<verticesCount;referenceVertex++){
		CVSymmetryInputChangeReferenceVertex(referenceVertex, inputParameters);
		
		//CVConcentricMergedInformationDestroy(inputParameters->mergedInformation);
		//inputParameters->mergedInformation = CVNewConcentricMergedInformation();
		
		outputArray[referenceVertex] = CVNewSymmetryOutput();
		
		CVNetworkCalculateSymmetryOfVertex(inputParameters, outputArray[referenceVertex]);
		
		if(currentProgress){
			CVAtomicIncrementInteger(currentProgress);
			if(updateCallback){
				updateCallback(operationControl);
			}
		}
		if(streamCallback){
			streamCallback(operationControl,referenceVertex,"%"CVFloatScan,0.1);//Value Streamming
		}
	}
	//ClearFiles
	
	if(verticesCount>0){
		CVSymmetryInputDestroy(inputParameters);
	}
	return CVTrue;
}

CVBool CVNetworkCalculateSymmetry(const CVNetwork* network, CVSize maximumLevel, CVBool saveExtraInformation, CVBool saveProbabilities, CVBool mergeLastLevel, CVSymmetryOutputParameters** outputArray, CVOperationControl* operationControl){
	CV_BenchmarkPrepare(CVNetworkCalculateAccessibility);
	CV_BenchmarkStart(CVNetworkCalculateAccessibility);
	
	CVBool returnValue;
	
#if CV_ENABLE_PARALLELISM
	CVInteger maxParallelBlocksCount = kCVDefaultParallelBlocks;
	CVSize problemSize = network->verticesCount;
	if(operationControl){
		maxParallelBlocksCount = operationControl->maxParallelBlocks;
	}
	if(network&&problemSize>=128&&maxParallelBlocksCount>1){
		returnValue = _CVNetworkCalculateSymmetry_parallel_implementation(network, maximumLevel, saveExtraInformation, saveProbabilities,mergeLastLevel, outputArray, operationControl);
	}else{
		returnValue = _CVNetworkCalculateSymmetry_implementation(network, maximumLevel, saveExtraInformation, saveProbabilities,mergeLastLevel, outputArray, operationControl);
	}
#else
	returnValue = _CVNetworkCalculateSymmetry_implementation(network, maximumLevel, saveExtraInformation, saveProbabilities,mergeLastLevel, outputArray, operationControl);
#endif //CV_ENABLE_PARALLELISM
	
	CV_BenchmarkStop(CVNetworkCalculateAccessibility);
	CV_BenchmarkPrint(CVNetworkCalculateAccessibility);
	return returnValue;
}