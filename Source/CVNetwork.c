//
//  CVNetwork.c
//
//
//  Created by Filipi Nascimento Silva on 9/22/12.
//
//
#include "CVNetwork.h"
#include "uthash.h"
#include "CVGridUtilities.h"


CVBool CVNetworkAddNewEdges(CVNetwork* network, CVIndex* fromIndices, CVIndex* toIndices, CVFloat* weights, CVSize count){
	CVIndex initialEdge = network->edgesCount;
	CVSize verticesCount = network->verticesCount;
	CVIndex i;

	CVNetworkGrowEdgesCount(network,count);

	for(i=0; i<count; i++){
		CVIndex toIndex   = toIndices[i];
		CVIndex fromIndex = fromIndices[i];
		if (toIndex >= verticesCount||fromIndex >= verticesCount) {
			return CVFalse;
		}
		network->edgeFromList[i+initialEdge]=fromIndex;
		network->edgeToList[i+initialEdge]=toIndex;

		CVNetworkGrowVertexSetEdgeForVertex(network,i+initialEdge,fromIndex,toIndex);

		if(network->edgeWeighted&&weights==NULL){
			network->edgesWeights[i+initialEdge]=1.0f;
		}else if(network->edgeWeighted){
			network->edgesWeights[i+initialEdge]=weights[i];
		}
		if(!network->directed){
			CVNetworkGrowVertexSetEdgeForVertex(network,i+initialEdge,toIndex,fromIndex);
		}else{
			CVNetworkGrowVertexSetInEdgeForVertex(network,i+initialEdge,toIndex,fromIndex);
		}
	}
	return CVTrue;
}



void CVNetworkPrint(const CVNetwork* network){
	printf("Vertices:" "%" CVSizeScan "\n", network->verticesCount);
	printf("Edges: " "%" CVSizeScan "\n", network->edgesCount);
	CVIndex fromVertex;
	for (fromVertex=0; fromVertex<=network->verticesCount; fromVertex++) {
		CVIndex i;
		CVSize toVerticesCount = network->vertexNumOfEdges[fromVertex];
		CVIndex* toVertices = network->vertexEdgesLists[fromVertex];
		printf("%"CVIndexScan"\t:",fromVertex);
		for (i=0; i<toVerticesCount; i++) {
			CVIndex toVertex = toVertices[i];
			printf("\t" "%"CVIndexScan,toVertex);
		}
		printf("\n");
	}
}

CVNetwork* CV_NewAllocationNetwork(CVSize verticesCount){
	CVNetwork* newNetwork = malloc(sizeof(CVNetwork));
	newNetwork->vertexNumOfEdges = calloc(verticesCount, sizeof(CVSize));
	newNetwork->vertexCapacityOfEdges = calloc(verticesCount, sizeof(CVSize));
	newNetwork->vertexEdgesLists = calloc(verticesCount, sizeof(CVSize*));
	newNetwork->vertexEdgesIndices = calloc(verticesCount, sizeof(CVSize*));

	newNetwork->vertexNumOfInEdges = calloc(verticesCount, sizeof(CVSize));
	newNetwork->vertexCapacityOfInEdges = calloc(verticesCount, sizeof(CVSize));
	newNetwork->vertexInEdgesLists = calloc(verticesCount, sizeof(CVSize*));
	newNetwork->vertexInEdgesIndices = calloc(verticesCount, sizeof(CVSize*));
	newNetwork->verticesWeights = calloc(verticesCount, sizeof(CVFloat));
	newNetwork->verticesEnabled = calloc(verticesCount, sizeof(CVBool));
	CVIndex i;
	for(i=0;i<verticesCount;i++){
		newNetwork->verticesWeights[i] = 1.0f;
		newNetwork->verticesEnabled[i] = CVTrue;
	}

	newNetwork->edgeFromList = NULL;
	newNetwork->edgeToList = NULL;

	newNetwork->edgesWeights = NULL;

	newNetwork->vertexCapacity = verticesCount;
	newNetwork->edgesCapacity = 0;

	newNetwork->edgesCount = 0;
	newNetwork->verticesCount = verticesCount;

	newNetwork->editable = CVFalse;
	newNetwork->directed = CVFalse;
	newNetwork->edgeWeighted = CVFalse;
	newNetwork->vertexWeighted = CVFalse;

	newNetwork->vertexNames = NULL;
	newNetwork->propertiesData = NULL;
	newNetwork->propertiesTypes = NULL;
	newNetwork->propertiesCount = 0;
	newNetwork->propertiesNames = NULL;
	
	return newNetwork;
}



void CV_NetworkDestroyProperties(CVNetwork* theNetwork){
	CVIndex i;
	for (i=0; i<theNetwork->propertiesCount; i++) {
		CVPropertyType type = theNetwork->propertiesTypes[i];
		
		if(type==CVStringPropertyType){
			CVIndex j;
			CVString* values = theNetwork->propertiesData[i];
			for (j=0; j<theNetwork->verticesCount; j++) {
				free(values[j]);
			}
		}
		free(theNetwork->propertiesData[i]);
		free(theNetwork->propertiesNames[i]);
	}
	free(theNetwork->propertiesData);
	free(theNetwork->propertiesNames);
	free(theNetwork->propertiesTypes);
	theNetwork->propertiesCount = 0;
	theNetwork->propertiesData= NULL;
	theNetwork->propertiesNames= NULL;
	theNetwork->propertiesTypes= NULL;
}


void CVNetworkAppendProperty(CVNetwork* theNetwork, CVString name, CVPropertyType type, void* data){
	theNetwork->propertiesCount++;
	CVIndex newIndex = theNetwork->propertiesCount-1;
	theNetwork->propertiesData  = realloc(theNetwork->propertiesData, sizeof(void*)*theNetwork->propertiesCount);
	theNetwork->propertiesNames = realloc(theNetwork->propertiesNames, sizeof(CVString)*theNetwork->propertiesCount);
	theNetwork->propertiesTypes = realloc(theNetwork->propertiesTypes, sizeof(CVPropertyType)*theNetwork->propertiesCount);
	
	theNetwork->propertiesTypes[newIndex] =type;
	theNetwork->propertiesNames[newIndex] = calloc(strlen(name)+1, sizeof(CVChar));
	strncpy(theNetwork->propertiesNames[newIndex], name, strlen(name));
	
	switch (type) {
		case CVNumberPropertyType:{
			CVFloat* values = calloc(theNetwork->verticesCount, sizeof(CVFloat));
			memcpy(values, data, sizeof(CVFloat)*theNetwork->verticesCount);
			theNetwork->propertiesData[newIndex] = values;
			break;
		}
		case CVStringPropertyType:{
			CVString* values = calloc(theNetwork->verticesCount, sizeof(CVString));
			CVString* stringValues = data;
			CVIndex i;
			for (i=0; i<theNetwork->verticesCount; i++) {
				values[i] = calloc(strlen(stringValues[i])+1, sizeof(CVChar));
				strncpy(values[i], stringValues[i], strlen(stringValues[i]));
			}
			theNetwork->propertiesData[newIndex] = values;
			break;
		}
		case CVVector2DPropertyType:{
			CVFloat* values = calloc(theNetwork->verticesCount*2, sizeof(CVFloat));
			memcpy(values, data, sizeof(CVFloat)*theNetwork->verticesCount*2);
			theNetwork->propertiesData[newIndex] = values;
			break;
		}
		case CVVector3DPropertyType:{
			CVFloat* values = calloc(theNetwork->verticesCount*3, sizeof(CVFloat));
			memcpy(values, data, sizeof(CVFloat)*theNetwork->verticesCount*3);
			theNetwork->propertiesData[newIndex] = values;
			break;
		}
		default:
			break;
	}
}

void* CVNetworkPropertyWithName(CVNetwork* network, CVString name, CVPropertyType* outType){
	CVIndex propertyIndex;
	void* data = NULL;
	for(propertyIndex=0;propertyIndex<network->propertiesCount;propertyIndex++){
		CVPropertyType type = network->propertiesTypes[propertyIndex];
		
		if(strcmp(network->propertiesNames[propertyIndex], name)==0){
			if(outType){
				*outType = type;
			}
			data = network->propertiesData[propertyIndex];
		}
	}
	return data;
}



/*
 if((token = strsep(&lineSegment, " ")) != NULL){
 }
 */

CVNetwork* CVNewNetwork(CVSize verticesCount, CVBool edgeWeighted, CVBool directed){
	CVNetwork * theNetwork = NULL;
	theNetwork = CV_NewAllocationNetwork(verticesCount);
	theNetwork->vertexWeighted = CVFalse;
	theNetwork->edgeWeighted = edgeWeighted;
	theNetwork->directed = directed;
	return theNetwork;
}

CVNetwork* CVNewNetworkWithNetwork(CVNetwork* originalNetwork, CVBool edgeWeighted, CVBool directed){
	CVNetwork * theNetwork = NULL;
	theNetwork = CV_NewAllocationNetwork(originalNetwork->verticesCount);
	theNetwork->vertexWeighted = CVFalse;
	theNetwork->edgeWeighted = edgeWeighted;
	theNetwork->directed = directed;
	CVIndex i;
	for(i=0;i<originalNetwork->edgesCount;i++){
		CVIndex from,to;
		from = originalNetwork->edgeFromList[i];
		to = originalNetwork->edgeToList[i];
		CVFloat weight = 1.0f;
		if(originalNetwork->edgeWeighted){
			weight = originalNetwork->edgesWeights[i];
		}
		CVNetworkAddNewEdge(theNetwork, from, to, weight);
	}
	CVIndex propertyIndex;
	for(propertyIndex=0;propertyIndex<originalNetwork->propertiesCount;propertyIndex++){
		CVPropertyType type = originalNetwork->propertiesTypes[propertyIndex];
		void* data = originalNetwork->propertiesData[propertyIndex];
		CVString name = originalNetwork->propertiesNames[propertyIndex];
		CVNetworkAppendProperty(theNetwork, name, type, data);
	}
	
	return theNetwork;
}

void CVNetworkWriteToFile(CVNetwork* theNetwork, FILE* networkFile){
	fprintf(networkFile,"#vertices ""%"CVSizeScan" nonweighted\n",theNetwork->verticesCount);
	if(theNetwork->vertexNames){
		CVIndex i;
		for(i=0;i<theNetwork->verticesCount;i++){
			fprintf(networkFile,"\"%s\"\n",theNetwork->vertexNames[i]);
		}
	}
	fprintf(networkFile,"#edges %s %s\n",theNetwork->edgeWeighted?"weighted":"nonweighted",theNetwork->directed?"directed":"undirected");
	CVIndex edgeIndex;
	
	CVIndex* edgesFrom = theNetwork->edgeFromList;
	CVIndex* edgesTo = theNetwork->edgeToList;
	CVFloat* edgesWeights = theNetwork->edgesWeights;
	for (edgeIndex=0; edgeIndex<theNetwork->edgesCount; edgeIndex++) {
		CVIndex fromVertex = edgesFrom[edgeIndex];
		CVIndex toVertex = edgesTo[edgeIndex];
		if(theNetwork->edgeWeighted){
			CVFloat weight = edgesWeights[edgeIndex];
			fprintf(networkFile,"%"CVIndexScan" ""%"CVIndexScan" ""%"CVFloatScan"\n",fromVertex,toVertex,weight);
		}else{
			fprintf(networkFile,"%"CVIndexScan" ""%"CVIndexScan"\n",fromVertex,toVertex);
		}
	}
	CVIndex propertyIndex;
	for(propertyIndex=0;propertyIndex<theNetwork->propertiesCount;propertyIndex++){
		CVPropertyType type = theNetwork->propertiesTypes[propertyIndex];
		void* data = theNetwork->propertiesData[propertyIndex];
		CVString name = theNetwork->propertiesNames[propertyIndex];
		switch (type) {
			case CVNumberPropertyType:{
				CVFloat* floatData = data;
				fprintf(networkFile,"#v \"%s\" n\n",name);
				CVIndex i;
				for(i=0;i<theNetwork->verticesCount;i++){
					fprintf(networkFile,"%"CVFloatScan"\n",floatData[i]);
				}
				break;
			}
			case CVStringPropertyType:{
				CVString* stringData = data;
				fprintf(networkFile,"#v \"%s\" s\n",name);
				CVIndex i;
				for(i=0;i<theNetwork->verticesCount;i++){
					fprintf(networkFile,"\"%s\"\n",stringData[i]);
				}
				break;
			}
			case CVVector2DPropertyType:{
				CVFloat* floatData = data;
				fprintf(networkFile,"#v \"%s\" v2\n",name);
				CVIndex i;
				for(i=0;i<theNetwork->verticesCount;i++){
					fprintf(networkFile,"%"CVFloatScan" %"CVFloatScan"\n",floatData[i*2],floatData[i*2+1]);
				}
				break;
			}
			case CVVector3DPropertyType:{
				CVFloat* floatData = data;
				fprintf(networkFile,"#v \"%s\" v3\n",name);
				CVIndex i;
				for(i=0;i<theNetwork->verticesCount;i++){
					fprintf(networkFile,"%"CVFloatScan" %"CVFloatScan" %"CVFloatScan"\n",floatData[i*3],floatData[i*3+1],floatData[i*3+2]);
				}
				break;
			}
			default:{
			}
		}
	}
}



CVNetwork* CVNewNetworkFromXNETFile(FILE* networkFile){
	CVSize verticesCount = 0;

	CVBool isReadingVertices = CVFalse;
	CVBool isReadingEdges = CVFalse;
	CVBool isReadingProperty = CVFalse;

	CVNetwork * theNetwork = NULL;
	CVSize* fromIndices = NULL;
	CVSize* toIndices = NULL;
	CVFloat* edgesWeights = NULL;
	CVSize edgesCount = 0;
	CVSize edgesCapacity = 0;
	CVIndex currentVertex = 0;
	CVString propertyName = NULL;
	CVPropertyType propertyType;
	CVIndex propertyVertexIndex = 0;
	void* propertyData;
	CVString parsingError = NULL;

	CVIndex currentLine = 0;
#define __CVMAX_LINE_BUFFER_SIZE 255
	CVString lineBuffer = malloc((__CVMAX_LINE_BUFFER_SIZE+1)*sizeof(CVChar));
	while(fgets(lineBuffer, __CVMAX_LINE_BUFFER_SIZE, networkFile)&&!feof(networkFile)){
		CVString lineSegment = lineBuffer;
		currentLine++;
		if(lineSegment&&CVStringScanCharacters(&lineSegment, '#')){
			//printf("Reading Line: %s\n",lineSegment);
			if(CVStringScan(&lineSegment, "vertices")){
				CVStringScanCharacters(&lineSegment, ' ');
				if(CVStringScanIndex(&lineSegment, &verticesCount)){
					CVStringScanCharacters(&lineSegment, ' ');
					//printf("VerticesCount: %ld\n", verticesCount);
					theNetwork = CV_NewAllocationNetwork(verticesCount);
				}
				if(CVStringScan(&lineSegment, "weighted")){
					theNetwork->vertexWeighted = CVTrue;
				}
				isReadingVertices=CVTrue;
				isReadingEdges=CVFalse;
				isReadingProperty = CVFalse;
				currentVertex = 0;
			}else if(CVStringScan(&lineSegment, "edges")){
				CVStringScanCharacters(&lineSegment, ' ');
				if(CVStringScan(&lineSegment, "weighted")){
					theNetwork->edgeWeighted = CVTrue;
				}
				CVStringScanCharacters(&lineSegment, ' ');
				if(CVStringScan(&lineSegment, "directed")){
					theNetwork->directed = CVTrue;
				}
				CVStringScanCharacters(&lineSegment, ' ');
				if(CVStringScan(&lineSegment, "weighted")){
					theNetwork->edgeWeighted = CVTrue;
				}
				isReadingVertices=CVFalse;
				isReadingEdges=CVTrue;
				isReadingProperty = CVFalse;
			}else if(CVStringScan(&lineSegment,"v")){
				CVStringScanCharacters(&lineSegment, ' ');
				CVStringScanCharacters(&lineSegment, '\"');
				free(propertyName);
				propertyName = CVNewStringScanningUpToCharacter(&lineSegment, '\"');
				propertyVertexIndex=0;
				CVStringScanCharacters(&lineSegment, '\"');
				CVStringScanCharacters(&lineSegment, ' ');
				if(CVStringScan(&lineSegment, "n")){
					propertyType=CVNumberPropertyType;
					isReadingProperty = CVTrue;
					propertyData = calloc(verticesCount,sizeof(CVFloat));
				}else if(CVStringScan(&lineSegment, "v2")){
					propertyType=CVVector2DPropertyType;
					isReadingProperty = CVTrue;
					propertyData = calloc(verticesCount*2,sizeof(CVFloat));
				}else if(CVStringScan(&lineSegment, "v3")){
					propertyType=CVVector3DPropertyType;
					isReadingProperty = CVTrue;
					propertyData = calloc(verticesCount*3,sizeof(CVFloat));
				}else if(CVStringScan(&lineSegment, "s")){
					propertyType=CVStringPropertyType;
					isReadingProperty = CVTrue;
					propertyData = calloc(verticesCount,sizeof(CVString));
				}else{
					isReadingProperty = CVFalse;
					//printf("Unnalocationg %s\n",propertyName);
					free(propertyName);
					parsingError = "Error reading header for property.";
					break;
				}
				isReadingVertices=CVFalse;
				isReadingEdges=CVFalse;
			}else{
				isReadingVertices=CVFalse;
				isReadingEdges=CVFalse;
				isReadingProperty = CVFalse;
			}
		}else{
			if(isReadingVertices){
				if(currentVertex<verticesCount){
					if(!theNetwork->vertexNames){
						theNetwork->vertexNames = calloc(verticesCount, sizeof(CVString));
					}
					CVStringTrim(lineSegment, "\"\n \t");
					CVSize lineLength = strlen(lineSegment);
					theNetwork->vertexNames[currentVertex] = calloc(lineLength+1, sizeof(CVChar));
					strncpy(theNetwork->vertexNames[currentVertex], lineSegment, lineLength);
					currentVertex++;
				}else{
					isReadingVertices=CVFalse;
				}
			}else if(isReadingEdges){
				unsigned long _longFromIndex = 0;
				unsigned long _longToIndex = 0;
				CVFloat _doubleWeight = 1.0;
				if(sscanf(lineSegment, "%ld %ld " "%"CVFloatScan,&_longFromIndex,&_longToIndex,&_doubleWeight)>=2){
					edgesCount++;
					if(CVUnlikely(edgesCapacity < edgesCount)){
						edgesCapacity = CVCapacityGrow(edgesCount);
						fromIndices = realloc(fromIndices, sizeof(CVIndex)*edgesCapacity);
						toIndices = realloc(toIndices, sizeof(CVIndex)*edgesCapacity);
						if(theNetwork->edgeWeighted){
							edgesWeights = realloc(edgesWeights, sizeof(CVFloat)*edgesCapacity);
						}
					}
					fromIndices[edgesCount-1]=_longFromIndex;
					toIndices[edgesCount-1]=_longToIndex;
					if(theNetwork->edgeWeighted){
						edgesWeights[edgesCount-1] = _doubleWeight;
					}
				}
			}else if(isReadingProperty){
				CVStringTrim(lineSegment, "\"\n \t");
				switch (propertyType) {
					case CVNumberPropertyType:{
						CVFloat currentValue = 0.0f;
						if(sscanf(lineSegment, "%"CVFloatScan,&currentValue)>0&&propertyVertexIndex<verticesCount){
							CVFloat* currentData = propertyData;
							currentData[propertyVertexIndex] = currentValue;
							propertyVertexIndex++;
							if(propertyVertexIndex==verticesCount){
								CVNetworkAppendProperty(theNetwork, propertyName, propertyType, currentData);
							}
						}
						break;
					}
					case CVStringPropertyType:{
						CVStringScanCharacters(&lineSegment, '\"');
						CVString currentString = CVNewStringScanningUpToCharacter(&lineSegment, '\"');
						
						CVString* currentData = propertyData;
						currentData[propertyVertexIndex] = currentString;
						propertyVertexIndex++;
						if(propertyVertexIndex==verticesCount){
							CVNetworkAppendProperty(theNetwork, propertyName, propertyType, currentData);
						}
						break;
					}
					case CVVector2DPropertyType:{
						CVFloat currentValue1 = 0.0f;
						CVFloat currentValue2 = 0.0f;
						if(sscanf(lineSegment, "%"CVFloatScan" %"CVFloatScan,&currentValue1,&currentValue2)>0&&propertyVertexIndex<verticesCount){
							CVFloat* currentData = propertyData;
							currentData[propertyVertexIndex*2] = currentValue1;
							currentData[propertyVertexIndex*2+1] = currentValue2;
							propertyVertexIndex++;
							if(propertyVertexIndex==verticesCount){
								CVNetworkAppendProperty(theNetwork, propertyName, propertyType, currentData);
							}
						}
						break;
					}
					case CVVector3DPropertyType:{
						CVFloat currentValue1 = 0.0f;
						CVFloat currentValue2 = 0.0f;
						CVFloat currentValue3 = 0.0f;
						if(sscanf(lineSegment, "%"CVFloatScan" %"CVFloatScan" %"CVFloatScan,&currentValue1,&currentValue2,&currentValue3)>0&&propertyVertexIndex<verticesCount){
							CVFloat* currentData = propertyData;
							currentData[propertyVertexIndex*3] = currentValue1;
							currentData[propertyVertexIndex*3+1] = currentValue2;
							currentData[propertyVertexIndex*3+2] = currentValue3;
							propertyVertexIndex++;
							if(propertyVertexIndex==verticesCount){
								CVNetworkAppendProperty(theNetwork, propertyName, propertyType, currentData);
							}
						}
						break;
					}
					default:
						break;
				}
			}

		}
	}
	
	if(parsingError==NULL){
		if(theNetwork && theNetwork->verticesCount>0){
			CVNetworkAddNewEdges(theNetwork, fromIndices,toIndices,edgesWeights, edgesCount);
		}
	}else{
		fprintf(stderr, "Parsing error occurred[at line %"CVIndexScan"]: %s\n", currentLine, parsingError);
		CVNetworkDestroy(theNetwork);
		theNetwork = NULL;
	}
	free(fromIndices);
	free(toIndices);
	free(edgesWeights);

	free(lineBuffer);
#undef __CVMAX_LINE_BUFFER_SIZE
	return theNetwork;
}

void CVNetworkDestroy(CVNetwork* network){
	CVIndex i;
	CV_NetworkDestroyProperties(network);
	
	
	
	
	
	for(i=0;i<network->verticesCount;i++){
		free(network->vertexEdgesLists[i]);
		free(network->vertexEdgesIndices[i]);

		//if(network->directed){
			free(network->vertexInEdgesLists[i]);
			free(network->vertexInEdgesIndices[i]);
			//}
	}

	free(network->vertexCapacityOfEdges);
	free(network->vertexNumOfEdges);
	free(network->vertexEdgesLists);
	free(network->vertexEdgesIndices);

	//if(network->directed){
		free(network->vertexNumOfInEdges);
		free(network->vertexInEdgesLists);
		free(network->vertexInEdgesIndices);
		free(network->vertexCapacityOfInEdges);
	//}

	free(network->edgeFromList);
	free(network->edgeToList);

	//if(network->edgeWeighted)
	free(network->edgesWeights);

	free(network->verticesWeights);
	free(network->verticesEnabled);
	
	free(network);
}

void CVNetworkGetDegree(const CVNetwork* network, CVIntegerArray* degrees){
	CVSize verticesCount = network->verticesCount;
	CVIntegerArrayReallocToCapacity(verticesCount, degrees);
	CVIntegerArraySetCount(verticesCount, degrees);
	
	CVIndex vertexIndex;
	for (vertexIndex=0; vertexIndex<verticesCount; vertexIndex++) {
		degrees->data[vertexIndex] = network->vertexNumOfEdges[vertexIndex];
	}
}

CVNetwork* CVNewRegularNetwork(CVSize* gridSize, CVSize dimension, CVFloat connectionRadius, CVBool toroidal){
	CVFloat maxGridSize = 0;
	CVIndex curDim;
	CVSize verticesCount = 1;
	for(curDim = 0; curDim<dimension;curDim++){
		maxGridSize = CVMAX(maxGridSize,gridSize[curDim]);
		verticesCount*=gridSize[curDim];
	}
	
	CVFloat maxRadius = CVMIN(connectionRadius, maxGridSize);
	
	CVInteger maxRadiusSize = ceil(maxRadius);
	CVSize inRadiusMaxCount = (CVSize)ipow(2*(int)maxRadiusSize+1,dimension);
	CVInteger* inRadiusIndices = calloc(inRadiusMaxCount*dimension, sizeof(CVInteger));
	CVFloat* inRadiusDistances = calloc(inRadiusMaxCount, sizeof(CVFloat));
	CVSize inRadiusCount = 0;
	
	CVIndex i;
	CVInteger* coordinateVector = calloc(dimension, sizeof(CVInteger));
	CVInteger* coordinateCenter = calloc(dimension, sizeof(CVInteger));
	CVSize* coordinateRadiusBox = calloc(dimension, sizeof(CVSize));
	
	for(i=0;i<dimension;i++){
		coordinateCenter[i] = maxRadiusSize;
		coordinateRadiusBox[i] = maxRadiusSize*2+1;
	}
	for (i=0; i<inRadiusMaxCount; i++) {
		CVGridGetCoordinatesFromLinearIndex(i, coordinateRadiusBox, dimension, coordinateVector);
		CVFloat distance = CVGridCalcIntegerDistance(coordinateCenter,coordinateVector, dimension);
		if(distance<=maxRadius){
			memcpy(inRadiusIndices+(inRadiusCount*dimension), coordinateVector, dimension*sizeof(CVInteger));
			inRadiusDistances[inRadiusCount] = distance;
			inRadiusCount++;
		}
	}
	CVNetwork* theNetwork = CVNewNetwork(verticesCount, CVFalse, CVFalse);
	
	CVInteger cellCount = verticesCount;
	
	CVInteger* newCoordinateVector = calloc(dimension, sizeof(CVInteger));
	CVIndex cellIndex;
	for(cellIndex=0;cellIndex<cellCount;cellIndex++){
		CVGridGetCoordinatesFromLinearIndex(cellIndex, gridSize, dimension, coordinateVector);
		CVIndex rIndex;
		for (rIndex=0; rIndex<inRadiusCount; rIndex++) {
			CVGridGetCoordinatesFromLinearIndex(cellIndex, gridSize, dimension, coordinateVector);
			//CVFloat distance = inRadiusDistances[rIndex];
			if(toroidal
			   ?CVGridGetDisplacedToroidalCoordinate(coordinateVector, inRadiusIndices+(rIndex*dimension),  coordinateCenter, gridSize, dimension,newCoordinateVector)
			   :CVGridGetDisplacedCoordinate(coordinateVector, inRadiusIndices+(rIndex*dimension),  coordinateCenter, gridSize, dimension,newCoordinateVector)){
				CVIndex displacedIndex = CVGridLinearIndexFromCoordinates(newCoordinateVector, gridSize, dimension);
				if(cellIndex<displacedIndex){
					CVNetworkAddNewEdge(theNetwork, cellIndex, displacedIndex, 1.0f);
				}
			}
			
		}
	}
	free(coordinateCenter);
	free(coordinateVector);
	free(inRadiusDistances);
	free(inRadiusIndices);
	return theNetwork;
}

CVNetwork* CVNewRegular2DNetwork(CVSize rows, CVSize columns, CVBool toroidal){
	CVSize verticesCount = rows*columns;
	CVSize maximumEdgesCount = verticesCount*2;
	CVSize edgesCount = 0;
	CVIndex* fromList = calloc(maximumEdgesCount, sizeof(CVIndex));
	CVIndex* toList = calloc(maximumEdgesCount, sizeof(CVIndex));
	CVFloat* positions = calloc(verticesCount*3, sizeof(CVFloat));
	CVIndex i,j;
	for (i=0; i<rows; i++) {
		for (j=0; j<columns; j++) {
			positions[(i*columns+j)*3+0]=(i-rows*0.5f)*200.0f/CVMAX(rows, columns);
			positions[(i*columns+j)*3+1]=(j-columns*0.5f)*200.0f/CVMAX(rows, columns);
			positions[(i*columns+j)*3+2]=0.0f;
			if(toroidal){
				fromList[edgesCount] = i*columns+j;
				toList[edgesCount] =(i)*columns+((j+1)%columns);
				edgesCount++;
				
				fromList[edgesCount] = i*columns+j;
				toList[edgesCount] =((i+1)%rows)*columns+(j);
				edgesCount++;
			}else{
				if(j+1<columns){
					fromList[edgesCount] = i*columns+j;
					toList[edgesCount] =(i)*columns+(j+1);
					edgesCount++;
				}
				if(i+1<rows){
					fromList[edgesCount] = i*columns+j;
					toList[edgesCount] =(i+1)*columns+(j);
					edgesCount++;
				}
			}
		}
	}
	CVNetwork* theNetwork = CVNewNetwork(verticesCount, CVFalse, CVFalse);
	CVNetworkAddNewEdges(theNetwork, fromList, toList, NULL, edgesCount);
	CVNetworkAppendProperty(theNetwork, "Position", CVVector3DPropertyType, positions);
	free(positions);
	free(fromList);
	free(toList);
	return theNetwork;
}


CVNetwork* CVNewRandomNetwork(CVSize verticesCount, CVFloat degree){
	CVSize averageEdgesCount = floorf(verticesCount*degree)+1;
	CVSize edgesCapacity = averageEdgesCount;
	CVIndex* fromList = calloc(edgesCapacity, sizeof(CVIndex));
	CVIndex* toList = calloc(edgesCapacity, sizeof(CVIndex));
	CVSize edgesCount = 0;
	
	double probability = degree/(double)verticesCount;
	
	CVIndex toIndex,fromIndex;
	for (fromIndex=0; fromIndex<verticesCount; fromIndex++) {
		for (toIndex=fromIndex+1; toIndex<verticesCount; toIndex++) {
			if (CVRandomFloat()<probability){
				if(edgesCapacity<edgesCount+1){
					edgesCapacity = CVCapacityGrow(edgesCount+1);
					fromList = realloc(fromList, sizeof(CVIndex)*edgesCapacity);
					toList = realloc(toList, sizeof(CVIndex)*edgesCapacity);
				}
				fromList[edgesCount] = fromIndex;
				toList[edgesCount] = toIndex;
				edgesCount++;
			}
		}
	}
	
	CVNetwork* theNetwork = CVNewNetwork(verticesCount, CVFalse, CVFalse);
	CVNetworkAddNewEdges(theNetwork, fromList, toList, NULL, edgesCount);
	free(fromList);
	free(toList);
	return theNetwork;
}


CVNetwork* CVNewWaxmanNetwork(CVSize verticesCount,CVFloat alpha, CVFloat beta, CVSize dimension){
	CVSize edgesCapacity = verticesCount*3;
	CVIndex* fromList = calloc(edgesCapacity, sizeof(CVIndex));
	CVIndex* toList = calloc(edgesCapacity, sizeof(CVIndex));
	CVSize edgesCount = 0;
	CVFloat* positions = calloc(verticesCount*dimension, sizeof(CVFloat));
	CVFloat* positions3D = calloc(verticesCount*3, sizeof(CVFloat));
	
	CVIndex vertexIndex;
	for (vertexIndex=0; vertexIndex<verticesCount; vertexIndex++) {
		CVIndex dimensionIndex;
		for (dimensionIndex=0; dimensionIndex<dimension; dimensionIndex++) {
			positions[dimension*vertexIndex+dimensionIndex] = CVRandomFloat();
			if(dimensionIndex<3){
				positions3D[3*vertexIndex+dimensionIndex] = (positions[dimension*vertexIndex+dimensionIndex] - 0.5)*200;
			}
		}
	}
	
	CVIndex toIndex,fromIndex;
	for (fromIndex=0; fromIndex<verticesCount; fromIndex++) {
		for (toIndex=fromIndex+1; toIndex<verticesCount; toIndex++) {
			double distanceSquared = 0.0;
			CVIndex dimensionIndex;
			for (dimensionIndex=0; dimensionIndex<dimension; dimensionIndex++) {
				double positionFrom = positions[dimension*fromIndex+dimensionIndex];
				double positionTo = positions[dimension*toIndex+dimensionIndex];
				distanceSquared += (positionFrom-positionTo)*(positionFrom-positionTo);
			}
			double probability = alpha*exp(-sqrt(distanceSquared)/(beta*sqrt(dimension)));
			if (CVRandomFloat()<probability){
				if(edgesCapacity<edgesCount+1){
					edgesCapacity = CVCapacityGrow(edgesCount+1);
					fromList = realloc(fromList, sizeof(CVIndex)*edgesCapacity);
					toList = realloc(toList, sizeof(CVIndex)*edgesCapacity);
				}
				fromList[edgesCount] = fromIndex;
				toList[edgesCount] = toIndex;
				edgesCount++;
			}
		}
	}
	
	CVNetwork* theNetwork = CVNewNetwork(verticesCount, CVFalse, CVFalse);
	CVNetworkAddNewEdges(theNetwork, fromList, toList, NULL, edgesCount);
	CVNetworkAppendProperty(theNetwork, "Position", CVVector3DPropertyType, positions3D);
	free(fromList);
	free(toList);
	free(positions);
	free(positions3D);
	return theNetwork;
}



CVNetwork* CVNewRandomGeographicNetwork(CVSize verticesCount, CVFloat maximumDistance, CVSize dimension){
	CVSize edgesCapacity = verticesCount*3;
	CVIndex* fromList = calloc(edgesCapacity, sizeof(CVIndex));
	CVIndex* toList = calloc(edgesCapacity, sizeof(CVIndex));
	CVSize edgesCount = 0;
	CVFloat* positions = calloc(verticesCount*dimension, sizeof(CVFloat));
	CVFloat* positions3D = calloc(verticesCount*3, sizeof(CVFloat));
	
	CVIndex vertexIndex;
	for (vertexIndex=0; vertexIndex<verticesCount; vertexIndex++) {
		CVIndex dimensionIndex;
		for (dimensionIndex=0; dimensionIndex<dimension; dimensionIndex++) {
			positions[dimension*vertexIndex+dimensionIndex] = CVRandomFloat();
			if(dimensionIndex<3){
				positions3D[3*vertexIndex+dimensionIndex] = (positions[dimension*vertexIndex+dimensionIndex] - 0.5)*200;
			}
		}
	}
	
	CVIndex toIndex,fromIndex;
	for (fromIndex=0; fromIndex<verticesCount; fromIndex++) {
		for (toIndex=fromIndex+1; toIndex<verticesCount; toIndex++) {
			double distanceSquared = 0.0;
			CVIndex dimensionIndex;
			for (dimensionIndex=0; dimensionIndex<dimension; dimensionIndex++) {
				double positionFrom = positions[dimension*fromIndex+dimensionIndex];
				double positionTo = positions[dimension*toIndex+dimensionIndex];
				distanceSquared += (positionFrom-positionTo)*(positionFrom-positionTo);
			}
			if (sqrt(distanceSquared)<maximumDistance){
				if(edgesCapacity<edgesCount+1){
					edgesCapacity = CVCapacityGrow(edgesCount+1);
					fromList = realloc(fromList, sizeof(CVIndex)*edgesCapacity);
					toList = realloc(toList, sizeof(CVIndex)*edgesCapacity);
				}
				fromList[edgesCount] = fromIndex;
				toList[edgesCount] = toIndex;
				edgesCount++;
			}
		}
	}
	
	CVNetwork* theNetwork = CVNewNetwork(verticesCount, CVFalse, CVFalse);
	CVNetworkAddNewEdges(theNetwork, fromList, toList, NULL, edgesCount);
	CVNetworkAppendProperty(theNetwork, "Position", CVVector3DPropertyType, positions3D);
	free(fromList);
	free(toList);
	free(positions);
	free(positions3D);
	return theNetwork;
}

CVNetwork* CVNewNetworkFromRandomRewiringEdgeList(CVIndex* fromList,CVIndex* toList, CVSize edgesCount, CVSize verticesCount, CVBool directed, CVFloat rewireProbability){
	CVIndex edgeIndex = 0;
	struct __cv_edge{
		CVIndex from;
		CVIndex to;
	};
	struct __cv_edge_element{
		struct __cv_edge edge;
		UT_hash_handle hh;
	};
	
	struct __cv_edge_element* edgesHash = NULL;
	
	//	HASH_ADD_KEYPTR(hh, edgesHash, edgesHash->edge, sizeof(struct __cv_edge), edgesHash);
	
	for (edgeIndex=0; edgeIndex<edgesCount; edgeIndex++) {
		struct __cv_edge_element* edgeElement = calloc(1, sizeof(struct __cv_edge_element));
		if(!directed){
			edgeElement->edge.from = CVMIN(fromList[edgeIndex],toList[edgeIndex]);
			edgeElement->edge.to = CVMAX(fromList[edgeIndex],toList[edgeIndex]);
		}else{
			edgeElement->edge.from = fromList[edgeIndex];
			edgeElement->edge.to = toList[edgeIndex];
		}
		HASH_ADD_KEYPTR(hh, edgesHash, (&(edgeElement->edge)), sizeof(struct __cv_edge), edgeElement);
	}
	
	for (edgeIndex=0; edgeIndex<edgesCount; edgeIndex++) {
		if(CVRandomFloat()<rewireProbability){
			CVBool edgeExists;
			do{
				struct __cv_edge_element* edgeElement = NULL;
				struct __cv_edge edgesKey;
				edgesKey.from = fromList[edgeIndex];
				edgesKey.to = CVRandomInRange(0, verticesCount);
				toList[edgeIndex] = edgesKey.to;
				HASH_FIND_PTR(edgesHash, &edgesKey, edgeElement);
				if(edgeElement || edgesKey.from==edgesKey.to){
					edgeExists=CVTrue;
				}else{
					edgeExists=CVFalse;
				}
			}while(edgeExists);
		}
	}
	
	CVNetwork* theNetwork = CVNewNetwork(verticesCount, CVFalse, directed);
	CVNetworkAddNewEdges(theNetwork, fromList, toList, NULL, edgesCount);
	
	struct __cv_edge_element* edgeElement, *tempElement;
	
	HASH_ITER(hh, edgesHash, edgeElement, tempElement) {
		HASH_DEL(edgesHash, edgeElement);
		free(edgeElement);
    }
	
	return theNetwork;
}

CVNetwork* CVNewNetworkFromRandomRewiring(const CVNetwork* originalNetwork, CVFloat rewiringProbability){
	CVIndex* fromList = calloc(originalNetwork->edgesCount, sizeof(CVIndex));
	CVIndex* toList = calloc(originalNetwork->edgesCount, sizeof(CVIndex));
	
	memcpy(fromList, originalNetwork->edgeFromList, sizeof(CVSize)*originalNetwork->edgesCount);
	memcpy(toList, originalNetwork->edgeToList, sizeof(CVSize)*originalNetwork->edgesCount);
	
	CVNetwork* theNetwork = CVNewNetworkFromRandomRewiringEdgeList(fromList, toList, originalNetwork->edgesCount, originalNetwork->verticesCount, originalNetwork->directed,rewiringProbability);
	
	free(fromList);
	free(toList);
	return theNetwork;
}


CVNetwork* CVNewWattsStrogatzNetwork(CVSize approximateNumberOfVertices, CVSize dimension, CVFloat connectionRadius, CVBool toroidal, CVFloat rewiringProbability){
	CVSize* gridSize = calloc(dimension, sizeof(CVSize));
	CVSize sizePerDimension = CVMAX(1,ceil(pow(approximateNumberOfVertices, 1.0/dimension)));
	CVIndex gridIndex;
	for (gridIndex=0; gridIndex<dimension; gridIndex++) {
		gridSize[gridIndex] = sizePerDimension;
	}
	CVNetwork* gridNetwork = CVNewRegularNetwork(gridSize, dimension, connectionRadius, toroidal);
	
	CVNetwork* wsNetwork = CVNewNetworkFromRandomRewiring(gridNetwork, rewiringProbability);
	free(gridSize);
	CVNetworkDestroy(gridNetwork);
	return wsNetwork;
}


CVNetwork* CVNewBarabasiAlbertNetwork(CVSize initialSize, CVSize degreeGrowth, CVSize iterations){
	CVSize edgesCount = iterations*degreeGrowth;
	CVIndex* fromList = calloc(edgesCount, sizeof(CVIndex));
	CVIndex* toList = calloc(edgesCount, sizeof(CVIndex));
	CVSize verticesCount = iterations+initialSize;
	CVIndex currentLink=0;
	CVIndex currentVertex;
	
	CVSize* distribGenerator = calloc(initialSize+(degreeGrowth*2)*iterations, sizeof(CVSize));
	CVIndex distribPointer=0;
	CVIndex i;
	for(i=0;i<initialSize;i++){
		distribGenerator[distribPointer]=i;
		distribPointer++;
	}
	
	for(currentVertex=initialSize;currentVertex < iterations+initialSize;currentVertex++){
		CVIndex m;
		for(m=0;m<degreeGrowth;m++){
			CVIndex connectTo = 0;
			CVBool linkExist = CVFalse;
			do{
				linkExist = CVFalse;
				connectTo = distribGenerator[CVRandomInRange(0, distribPointer-1)];
				CVIndex curLink;
				for(curLink=0;curLink<currentLink;curLink++){
					if((fromList[curLink]==currentVertex&&toList[curLink]==connectTo) ||
					   (fromList[curLink]==connectTo&&toList[curLink]==currentVertex)    ){
						linkExist=CVTrue;
					}
				}
			}while(linkExist);
			fromList[currentLink] = currentVertex;
			toList[currentLink] = connectTo;
			currentLink++;
			distribGenerator[distribPointer]=connectTo;
			distribPointer++;
		}
		for(m=0;m<degreeGrowth;m++){
			distribGenerator[distribPointer]=currentVertex;
			distribPointer++;
		}
	}
	
	CVNetwork* theNetwork = CVNewNetwork(verticesCount, CVFalse, CVFalse);
	CVNetworkAddNewEdges(theNetwork, fromList, toList, NULL, edgesCount);
	free(fromList);
	free(toList);
	free(distribGenerator);
	return theNetwork;
}

CVBool CVNetworkCouldBeIsomorphic(const CVNetwork* aNetwork,const CVNetwork* bNetwork){
	if(aNetwork->directed||bNetwork->directed){ //DIRECTED NOT SUPPORTED
		return CVFalse;
	}
	if(aNetwork->vertexWeighted||bNetwork->vertexWeighted){ //DIRECTED NOT SUPPORTED
		return CVFalse;
	}
	if(aNetwork->verticesCount!=bNetwork->verticesCount){
		return CVFalse;
	}
	
	if(aNetwork->edgesCount!=bNetwork->edgesCount){
		return CVFalse;
	}
	
	CVIntegerArray aDegrees;
	CVIntegerArray bDegrees;
	
	CVIntegerArrayInitWithCapacity(1, &aDegrees);
	CVIntegerArrayInitWithCapacity(1, &bDegrees);
	
	CVNetworkGetDegree(aNetwork, &aDegrees);
	CVNetworkGetDegree(bNetwork, &bDegrees);
	
	CVIntegerArrayQuickSort3(&aDegrees);
	CVIntegerArrayQuickSort3(&bDegrees);
	
	CVIndex i;
	CVBool degreeOk = CVTrue;;
	for (i=0; i<aNetwork->verticesCount; i++) {
		if(aDegrees.data[i]!=bDegrees.data[i]){
			degreeOk = CVFalse;
		}
	}
	
	
	if(!degreeOk){
		CVIntegerArrayDestroy(&aDegrees);
		CVIntegerArrayDestroy(&bDegrees);
		return CVFalse;
	}
	
	
	//Neighbors
	CVIntegerArray aEdgeDegrees;
	CVIntegerArray bEdgeDegrees;
	
	CVNetworkGetDegree(aNetwork, &aDegrees);
	CVNetworkGetDegree(bNetwork, &bDegrees);
	
	CVIntegerArrayInitWithCapacity(1, &aEdgeDegrees);
	CVIntegerArrayInitWithCapacity(1, &bEdgeDegrees);
	
	CVBool edgesDegreeOk = CVTrue;
	for (i=0; i<aNetwork->edgesCount; i++) {
		CVInteger aEdgeDegree = aDegrees.data[aNetwork->edgeFromList[i]]+aDegrees.data[aNetwork->edgeToList[i]];
		CVInteger bEdgeDegree = bDegrees.data[bNetwork->edgeFromList[i]]+bDegrees.data[bNetwork->edgeToList[i]];
		
		CVIntegerArrayAdd(aEdgeDegree, &aEdgeDegrees);
		CVIntegerArrayAdd(bEdgeDegree, &bEdgeDegrees);
	}
	
	CVIntegerArrayQuickSort3(&aEdgeDegrees);
	CVIntegerArrayQuickSort3(&bEdgeDegrees);
	
	for (i=0; i<aNetwork->edgesCount; i++) {
		if(aEdgeDegrees.data[i]!=bEdgeDegrees.data[i]){
			edgesDegreeOk = CVFalse;
		}
	}
	
	CVIntegerArrayDestroy(&aDegrees);
	CVIntegerArrayDestroy(&bDegrees);
	CVIntegerArrayDestroy(&aEdgeDegrees);
	CVIntegerArrayDestroy(&bEdgeDegrees);
	
	if(!edgesDegreeOk){
		return CVFalse;
	}
	
	
	
	return CVTrue;
}





CVBool CVNetworkAddNewEdge(CVNetwork* network, CVIndex fromIndex, CVIndex toIndex, CVFloat weight){
	CVIndex initialEdge = network->edgesCount;
	CVSize verticesCount = network->verticesCount;
	if (toIndex >= verticesCount||fromIndex >= verticesCount) {
		return CVFalse;
	}
	CVNetworkGrowEdgesCapacity(network,1);
	network->edgeFromList[initialEdge]=fromIndex;
	network->edgeToList[initialEdge]=toIndex;
	
	CVNetworkGrowVertexSetEdgeForVertex(network,initialEdge,fromIndex,toIndex);
	
	if(network->edgeWeighted&&weight>=0){
		network->edgesWeights[initialEdge]=weight;
	}else if(network->edgeWeighted){
		network->edgesWeights[initialEdge]=1.0f;
	}
	if(!network->directed){
		//printf("Index: %lu toIndex:%lu fromIndex:%lu\n",i+initialEdge,toIndex,fromIndex);
		CVNetworkGrowVertexSetEdgeForVertex(network,initialEdge,toIndex,fromIndex);
		//printf("OK\n");
	}else{
		CVNetworkGrowVertexSetInEdgeForVertex(network,initialEdge,toIndex,fromIndex);
	}
	network->edgesCount++;
	return CVTrue;
}

CVBool CVNetworkAddNewEdgeAndIntegrateWeight(CVNetwork* network, CVIndex fromIndex, CVIndex toIndex, CVFloat weight){
	CVIndex initialEdge = network->edgesCount;
	CVSize verticesCount = network->verticesCount;
	if (toIndex >= verticesCount||fromIndex >= verticesCount) {
		return CVFalse;
	}
	
	CVBool edgeFound = CVFalse;
	CVIndex i;
	CVSize toVerticesCount = network->vertexNumOfEdges[fromIndex];
	CVIndex* toVertices = network->vertexEdgesLists[fromIndex];
	for (i=0; i<toVerticesCount; i++) {
		if(toVertices[i]==toIndex){
			edgeFound = CVTrue;
			break;
		}
	}
	
	if(edgeFound){
		if(network->edgeWeighted&&weight>0){
			network->edgesWeights[network->vertexEdgesIndices[fromIndex][i]]+=weight;
		}
	}else{
		CVNetworkGrowEdgesCapacity(network,1);
		network->edgeFromList[initialEdge]=fromIndex;
		network->edgeToList[initialEdge]=toIndex;
		
		CVNetworkGrowVertexSetEdgeForVertex(network,initialEdge,fromIndex,toIndex);
		
		if(network->edgeWeighted&&weight>=0){
			network->edgesWeights[initialEdge]=weight;
		}else if(network->edgeWeighted){
			network->edgesWeights[initialEdge]=1.0f;
		}
		if(!network->directed){
			//printf("Index: %lu toIndex:%lu fromIndex:%lu\n",i+initialEdge,toIndex,fromIndex);
			//FIXME: Directed Networks
			CVNetworkGrowVertexSetEdgeForVertex(network,initialEdge,toIndex,fromIndex);
			//printf("OK\n");
		}else{
			CVNetworkGrowVertexSetInEdgeForVertex(network,initialEdge,toIndex,fromIndex);
		}
		network->edgesCount++;
	}
	return CVTrue;
}


void CVNetworkWriteToGMLFile(CVNetwork* theNetwork, FILE* networkFile){
	fprintf(networkFile, "graph [\n");
	if(theNetwork->directed){
		fprintf(networkFile, "  directed 1\n");
	}
	
	CVIndex i;
	for(i=0;i<theNetwork->verticesCount;i++){
		fprintf(networkFile,"  node [\n");
		
		fprintf(networkFile,"    id %"CVIndexScan"\n",i);
		if(theNetwork->vertexNames){
			fprintf(networkFile,"    label \"%s\"\n",theNetwork->vertexNames[i]);
		}
		
		CVIndex propertyIndex;
		for(propertyIndex=0;propertyIndex<theNetwork->propertiesCount;propertyIndex++){
			CVPropertyType type = theNetwork->propertiesTypes[propertyIndex];
			void* data = theNetwork->propertiesData[propertyIndex];
			CVString name = theNetwork->propertiesNames[propertyIndex];
			switch (type) {
				case CVNumberPropertyType:{
					CVFloat* floatData = data;
					fprintf(networkFile,"    ");
					CVBool nextUpper = CVFalse;
					while (*name) {
						if(isalnum(*name)){
							fputc(nextUpper?toupper(*name):*name, networkFile);
							nextUpper = CVFalse;
						}else if(isspace(*name)){
							nextUpper=CVTrue;
						}
						name++;
					}
					fprintf(networkFile," ");
					fprintf(networkFile,"%"CVFloatScan"\n",floatData[i]);
					break;
				}
				case CVStringPropertyType:{
					CVString* stringData = data;
					fprintf(networkFile,"    ");
					CVBool nextUpper = CVFalse;
					while (*name) {
						if(isalnum(*name)){
							fputc(nextUpper?toupper(*name):*name, networkFile);
							nextUpper = CVFalse;
						}else if(isspace(*name)){
							nextUpper=CVTrue;
						}
						name++;
					}
					fprintf(networkFile," ");
					fprintf(networkFile,"\"%s\"\n",stringData[i]);
					break;
				}
				case CVVector2DPropertyType:{
					CVFloat* floatData = data;
					fprintf(networkFile,"    ");
					CVBool nextUpper = CVFalse;
					while (*name) {
						if(isalnum(*name)){
							fputc(nextUpper?toupper(*name):*name, networkFile);
							nextUpper = CVFalse;
						}else if(isspace(*name)){
							nextUpper=CVTrue;
						}
						name++;
					}
					fprintf(networkFile," ");
					fprintf(networkFile,"[ x %"CVFloatScan" y %"CVFloatScan" ]\n",floatData[i*2],floatData[i*2+1]);
					break;
				}
				case CVVector3DPropertyType:{
					CVFloat* floatData = data;
					if(strcmp(name, "Position")==0 || strcmp(name, "position")==0 ){
						fprintf(networkFile,"    graphics");
					}else{
						fprintf(networkFile,"    ");
						CVBool nextUpper = CVFalse;
						while (*name) {
							if(isalnum(*name)){
								fputc(nextUpper?toupper(*name):*name, networkFile);
								nextUpper = CVFalse;
							}else if(isspace(*name)){
								nextUpper=CVTrue;
							}
							name++;
						}
					}
					fprintf(networkFile," ");
					fprintf(networkFile,"[ x %"CVFloatScan" y %"CVFloatScan" z %"CVFloatScan" ]\n",floatData[i*3],floatData[i*3+1],floatData[i*3+2]);
					break;
				}
				default:{
				}
			}
		}
		
		fprintf(networkFile,"  ]\n");
	}
	CVIndex edgeIndex;
	
	CVIndex* edgesFrom = theNetwork->edgeFromList;
	CVIndex* edgesTo = theNetwork->edgeToList;
	CVFloat* edgesWeights = theNetwork->edgesWeights;
	for (edgeIndex=0; edgeIndex<theNetwork->edgesCount; edgeIndex++) {
		fprintf(networkFile,"  edge [\n");
		CVIndex fromVertex = edgesFrom[edgeIndex];
		CVIndex toVertex = edgesTo[edgeIndex];
		if(theNetwork->edgeWeighted){
			CVFloat weight = edgesWeights[edgeIndex];
			fprintf(networkFile,"    source %"CVIndexScan"\n    target %"CVIndexScan" \n    weight %"CVFloatScan"\n",fromVertex,toVertex,weight);
		}else{
			fprintf(networkFile,"    source %"CVIndexScan"\n    target %"CVIndexScan"\n",fromVertex,toVertex);
		}
		fprintf(networkFile,"  ]\n");
	}
	
	fprintf(networkFile,"]\n");
}


CVFloat CVNetworkClusteringCoefficient(const CVNetwork* aNetwork, CVIndex nodeIndex){
	CVSize vertexEdgesCount = aNetwork->vertexNumOfEdges[nodeIndex];
	CVIndex* vertexEdgesList = aNetwork->vertexEdgesLists[nodeIndex];
	CVSize inLevelConnections = 0;
	CVIndex ni;
	CVBitArray isNeighbor = CVNewBitArray(aNetwork->verticesCount);
	for(ni=0;ni<vertexEdgesCount;ni++){
		CVBitArraySet(isNeighbor, vertexEdgesList[ni]);
	}
	for(ni=0;ni<vertexEdgesCount;ni++){
		CVIndex neighborVertex = vertexEdgesList[ni];
		CVSize neighborEdgesCount = aNetwork->vertexNumOfEdges[neighborVertex];
		CVIndex* neighborEdgesList = aNetwork->vertexEdgesLists[neighborVertex];
		CVIndex nni;
		for(nni=0;nni<neighborEdgesCount;nni++){
			if(CVBitArrayTest(isNeighbor,neighborEdgesList[nni])){
				inLevelConnections++;
			}
		}
	}
	CVBitArrayDestroy(isNeighbor);
	if((vertexEdgesCount-1.0) > 0.0){
		return (inLevelConnections)/(CVFloat)(vertexEdgesCount*(vertexEdgesCount-1.0f));
	}else{
		return 0.0f;
	}
}

