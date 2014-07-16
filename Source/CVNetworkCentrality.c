//
//  CVNetworkCentrality.c
//  CVNetwork
//
//  Created by Filipi Nascimento Silva on 8/27/13.
//  Copyright (c) 2013 Filipi Nascimento Silva. All rights reserved.
//

#include "CVNetworkCentrality.h"
#include "fib.h"
#include <float.h>


#if CV_ENABLE_PARALLELISM

CVBool CVNetworkCalculateCentrality_weighted_parallel_implementation(const CVNetwork* network,CVFloatArray* centrality, CVOperationControl* operationControl){
	CVSize verticesCount = network->verticesCount;
	
	CVSize unrolledLoops = kCVDefaultParallelBlocks;
	CVInteger* currentProgress = NULL;
	void (*updateCallback)(CVOperationControl*)  = NULL;
	if(operationControl){
		operationControl->maxProgress = verticesCount;
		operationControl->currentProgress = 0;
		currentProgress = &(operationControl->currentProgress);
		if(operationControl->maxParallelBlocks>0){
			unrolledLoops = operationControl->maxParallelBlocks;
		}
		updateCallback = operationControl->updateCallback;
	}
	
	
	CVFloatArrayReallocToCapacity(verticesCount, centrality);
	CVFloatArraySetCount(verticesCount, centrality);
		//memset(centrality->data, 0, verticesCount*sizeof(CVFloat));
	CVIndex i;
	for (i=0; i<verticesCount; i++) {
		if(network->verticesEnabled[i]){
			centrality->data[i]  = 0.0f;
		}
	}
	CVSize unrolledSize = 1 + ((verticesCount - 1) / unrolledLoops);
	
	CVParallelForStart(centralityLoop, blockIndex, unrolledLoops){
		
		const CVFloat* verticesWeights = network->verticesWeights;
		const CVBool* verticesEnabled = network->verticesEnabled;
		
		CVUIntegerArray* P = calloc(verticesCount, sizeof(CVUIntegerArray));
		double* sigma = calloc(verticesCount, sizeof(double));
		double* d = calloc(verticesCount, sizeof(double));
		double* delta = calloc(verticesCount, sizeof(double));
		
		double* localCentrality = calloc(verticesCount, sizeof(double));
		
		double* seen = calloc(verticesCount, sizeof(double));
		
		CVIndex i;
		for (i=0; i<verticesCount; i++) {
			CVUIntegerArray newArray;
			CVUIntegerArrayInitWithCapacity(1, &newArray);
			P[i] = newArray;
		}
		
		CVIndex s;
		CVIntegerStack S = CVIntegerStackMake();
		struct fibheap* Q = fh_makekeyheap();
		
		CVSize maxIndex = CVMIN((blockIndex+1)*unrolledSize, verticesCount);
		for(s=blockIndex*unrolledSize;s<maxIndex;s++){
			if(currentProgress){
				CVAtomicIncrementInteger(currentProgress);
				if(updateCallback){
					updateCallback(operationControl);
				}
			}
			if(CVUnlikely(!verticesEnabled[s])){
				continue;
			}
			CVFloat sWeight = verticesWeights[s];
			S.count = 0;
			for (i=0; i<verticesCount; i++) {
				P[i].count=0;
				d[i] = -1.0;
				sigma[i] = 0.0;
				delta[i] = 0.0;
				seen[i] = -1.0;
			}
			sigma[s] = 1.0f;
			seen[s] = 0.0f;
			CVInteger v;
			CVInteger prev;
			fh_data sdata;
			sdata.data = s;
			sdata.prev = s;
			fh_insertkey(Q,0.0,sdata);
			
			fh_data vdata;
			double dist;
			while (fh_dequeue(Q,&vdata,&dist)) {
				v = vdata.data;
				prev = vdata.prev;
				
				if(d[v]!=-1.0){
					continue;
				}
				sigma[v]+=sigma[prev];
				CVIntegerStackPush(v, &S);
				d[v] = dist;
				
				CVSize vEdgesCount = network->vertexNumOfEdges[v];
				CVIndex* vNeighbor = network->vertexEdgesLists[v];
				CVIndex* vEdge = network->vertexEdgesIndices[v];
				CVFloat* edgesWeights = network->edgesWeights;
				CVIndex w,e;
				for (e=0; e<vEdgesCount; e++) {
					w = vNeighbor[e];
					double weight = exp(-edgesWeights[vEdge[e]]);
					if(CVLikely(verticesEnabled[w])){
						double vwdist = d[v]+weight;
						if(d[w]<0.0 && ( seen[w]<0.0 || vwdist < seen[w])){
							seen[w]= vwdist;
							fh_data wdata;
							wdata.data = w;
							wdata.prev = v;
							fh_insertkey(Q, vwdist, wdata);
							sigma[w] = 0.0;
							P[w].count = 0;
							CVUIntegerArrayAdd(v, P+w);
						}else{
							if(vwdist==seen[w]){
								sigma[w]+= sigma[v];
								CVUIntegerArrayAdd(v, P+w);
							}
						}
					}
				}
			}
			while(S.count>0){
				CVIndex w = CVIntegerStackPop(&S);
				CVIndex v,p;
				CVSize PCount = P[w].count;
				CVUInteger* PData = P[w].data;
				for (p=0; p<PCount; p++) {
					v = PData[p];
					delta[v] += sigma[v]/sigma[w]*(1.0+delta[w]);
				}
				if(w!=s){
					localCentrality[w] += sWeight*delta[w];
				}
			}
		}
		for (i=0; i<verticesCount; i++){
			CVUIntegerArrayDestroy(P+i);
		}
		free(P);
		free(sigma);
		free(d);
		free(delta);
		free(seen);
		CVIntegerArrayDestroy(&S);
		fh_deleteheap(Q);
		
		CVParallelLoopCriticalRegionStart(centralityLoop){
			CVFloat* centralityData = centrality->data;
			for (i=0; i<verticesCount; i++) {
				centralityData[i] += localCentrality[i];
			}
		}
		CVParallelLoopCriticalRegionEnd(centralityLoop);
		free(localCentrality);
	}CVParallelForEnd(centralityLoop);
	
	return CVTrue;
}

#endif //CV_ENABLE_PARALLELISM



CVBool CVNetworkCalculateCentrality_weighted_implementation(const CVNetwork* network,CVFloatArray* centrality, CVOperationControl* operationControl){
	CVSize verticesCount = network->verticesCount;
	
	CVInteger* currentProgress = NULL;
	void (*updateCallback)(CVOperationControl*)  = NULL;
	if(operationControl){
		operationControl->maxProgress = verticesCount;
		operationControl->currentProgress = 0;
		currentProgress = &(operationControl->currentProgress);
		updateCallback = operationControl->updateCallback;
	}
	
	CVFloatArrayReallocToCapacity(verticesCount, centrality);
	CVFloatArraySetCount(verticesCount, centrality);
		//memset(centrality->data, 0, verticesCount*sizeof(CVFloat));
	CVIndex i;
	for (i=0; i<verticesCount; i++) {
		if(network->verticesEnabled[i]){
			centrality->data[i]  = 0.0;
		}
	}
	
	const CVFloat* verticesWeights = network->verticesWeights;
	const CVBool* verticesEnabled = network->verticesEnabled;
	
	CVFloat* centralityData = centrality->data;
	
	CVUIntegerArray* P = calloc(verticesCount, sizeof(CVUIntegerArray));
	double* sigma = calloc(verticesCount, sizeof(double));
	double* d = calloc(verticesCount, sizeof(double));
	double* delta = calloc(verticesCount, sizeof(double));
	double* seen = calloc(verticesCount, sizeof(double));
	
	for (i=0; i<verticesCount; i++) {
		CVUIntegerArray newArray;
		CVUIntegerArrayInitWithCapacity(1, &newArray);
		P[i] = newArray;
	}
	CVIndex s;
	CVIntegerStack S = CVIntegerStackMake();
	struct fibheap* Q = fh_makekeyheap();
	
	
	for (s=0; s<verticesCount; s++) {
		if(currentProgress){
			CVAtomicIncrementInteger(currentProgress);
			if(updateCallback){
				updateCallback(operationControl);
			}
		}
		if(CVUnlikely(!verticesEnabled[s])){
			continue;
		}
		CVFloat sWeight = verticesWeights[s];
			//printf("%ld/%ld\n",s,verticesCount);
		S.count = 0;
		for (i=0; i<verticesCount; i++) {
			P[i].count=0;
			d[i] = -1.0;
			sigma[i] = 0.0;
			delta[i] = 0.0;
			seen[i] = -1.0;
		}
		sigma[s] = 1.0f;
		seen[s] = 0.0f;
		CVInteger v;
		CVInteger prev;
		fh_data sdata;
		sdata.data = s;
		sdata.prev = s;
		fh_insertkey(Q,0.0,sdata);
		
		fh_data vdata;
		double dist;
		while (fh_dequeue(Q,&vdata,&dist)) {
			v = vdata.data;
			prev = vdata.prev;
			
			if(d[v]!=-1.0){
				continue;
			}
			sigma[v]+=sigma[prev];
			CVIntegerStackPush(v, &S);
			d[v] = dist;
			
			CVSize vEdgesCount = network->vertexNumOfEdges[v];
			CVIndex* vNeighbor = network->vertexEdgesLists[v];
			CVIndex* vEdge = network->vertexEdgesIndices[v];
			CVFloat* edgesWeights = network->edgesWeights;
			CVIndex w,e;
			for (e=0; e<vEdgesCount; e++) {
				w = vNeighbor[e];
				double weight = exp(-edgesWeights[vEdge[e]]);
				if(CVLikely(verticesEnabled[w])){
					double vwdist = d[v]+weight;
					if(d[w]<0.0 && ( seen[w]<0.0 || vwdist < seen[w])){
						seen[w]= vwdist;
						fh_data wdata;
						wdata.data = w;
						wdata.prev = v;
						fh_insertkey(Q, vwdist, wdata);
						sigma[w] = 0.0;
						P[w].count = 0;
						CVUIntegerArrayAdd(v, P+w);
					}else{
						if(vwdist==seen[w]){
							sigma[w]+= sigma[v];
							CVUIntegerArrayAdd(v, P+w);
						}
					}
				}
			}
		}
		while(S.count>0){
			CVIndex w = CVIntegerStackPop(&S);
			CVIndex v,p;
			CVSize PCount = P[w].count;
			CVUInteger* PData = P[w].data;
			for (p=0; p<PCount; p++) {
				v = PData[p];
				delta[v] += sigma[v]/sigma[w]*(1.0+delta[w]);
					//if(v==2538){
						//printf("delta[v]: %f  sigma[v]: %f  sigma[w]: %f delta[w]: %f  s:%"CVIndexScan"  v:%"CVIndexScan"  w:%"CVIndexScan"\n",delta[v],sigma[v],sigma[w],delta[w], s,v,w);
						//}
			}
			if(w!=s){
				centralityData[w] += sWeight*delta[w];
					//if(w==2538&&s<=338){
					//	printf("sWeight: %f  delta[w]: %f  centralityData[w]: %f  s:%"CVIndexScan"\n",sWeight,delta[w],centralityData[w], s);
					//}
			}
		}
	}
		//CVQueueDestroy(&Q);
	CVIntegerArrayDestroy(&S);
	free(P);
	free(sigma);
	free(d);
	free(delta);
	free(seen);
	
	fh_deleteheap(Q);
	
	/*
	 printf("----\nFinished\n-----\n");
	 CVSize N = network->verticesCount;
	 for (i=0; i<verticesCount; i++) {
	 printf("%g\n",centralityData[i]/(N-1)/(N-2));
	 }
	 */
	return CVTrue;
}







#if CV_ENABLE_PARALLELISM
CVBool CVNetworkCalculateCentrality_parallel_implementation(const CVNetwork* network,CVFloatArray* centrality, CVOperationControl* operationControl){
	CVSize verticesCount = network->verticesCount;
	
	CVSize unrolledLoops = kCVDefaultParallelBlocks;
	CVInteger* currentProgress = NULL;
	void (*updateCallback)(CVOperationControl*)  = NULL;
	if(operationControl){
		operationControl->maxProgress = verticesCount;
		operationControl->currentProgress = 0;
		currentProgress = &(operationControl->currentProgress);
		if(operationControl->maxParallelBlocks>0){
			unrolledLoops = operationControl->maxParallelBlocks;
		}
		updateCallback = operationControl->updateCallback;
	}
	
	
	CVFloatArrayReallocToCapacity(verticesCount, centrality);
	CVFloatArraySetCount(verticesCount, centrality);
		//memset(centrality->data, 0, verticesCount*sizeof(CVFloat));
	CVIndex i;
	for (i=0; i<verticesCount; i++) {
		if(network->verticesEnabled[i]){
			centrality->data[i]  = 0.0f;
		}
	}
	CVSize unrolledSize = 1 + ((verticesCount - 1) / unrolledLoops);
	
	CVParallelForStart(centralityLoop, blockIndex, unrolledLoops){
		
		
		const CVFloat* verticesWeights = network->verticesWeights;
		const CVBool* verticesEnabled = network->verticesEnabled;
		
		CVUIntegerArray* P = calloc(verticesCount, sizeof(CVUIntegerArray));
		CVInteger* sigma = calloc(verticesCount, sizeof(CVInteger));
		CVInteger* d = calloc(verticesCount, sizeof(CVInteger));
		double* delta = calloc(verticesCount, sizeof(double));
		
		double* localCentrality = calloc(verticesCount, sizeof(double));
		
		CVIndex i;
		for (i=0; i<verticesCount; i++) {
			CVUIntegerArray newArray;
			CVUIntegerArrayInitWithCapacity(1, &newArray);
			P[i] = newArray;
		}
		
		CVIndex s;
		CVIntegerStack S = CVIntegerStackMake();
		CVQueue Q = CVQueueCreate();
		
		CVSize maxIndex = CVMIN((blockIndex+1)*unrolledSize, verticesCount);
		for(s=blockIndex*unrolledSize;s<maxIndex;s++){
			if(currentProgress){
				CVAtomicIncrementInteger(currentProgress);
				if(updateCallback){
					updateCallback(operationControl);
				}
			}
			if(CVUnlikely(!verticesEnabled[s])){
				continue;
			}
			CVFloat sWeight = verticesWeights[s];
			S.count = 0;
			for (i=0; i<verticesCount; i++) {
				P[i].count=0;
				d[i] = CVIntegerMAX;
				sigma[i] = 0;
				delta[i] = 0;
			}
			sigma[s] = 1;
			d[s] = 0;
			CVInteger v;
			CVQueuePush(&Q, s);
			while (CVQueueDequeue(&Q,&v)) {
				CVIntegerStackPush(v, &S);
				CVSize vEdgesCount = network->vertexNumOfEdges[v];
				CVIndex* vNeighbor = network->vertexEdgesLists[v];
				CVIndex w,e;
				for (e=0; e<vEdgesCount; e++) {
					w = vNeighbor[e];
					if(CVLikely(verticesEnabled[w])){
						if(d[w]==CVIntegerMAX){
							d[w] = d[v] + 1;// FIXME: Change to the w-v weight;
								//fh_insertkey(Q, d[v], w);
								//if(d[w]<=2){
							CVQueuePush(&Q, w);
								//}
						}
						if(d[w]== d[v]+1){
							sigma[w]+= sigma[v];
							CVUIntegerArrayAdd(v, &(P[w]));
						}
					}
				}
			}
			while(S.count>0){
				CVIndex w = CVIntegerStackPop(&S);
				CVIndex v,p;
				CVSize PCount = P[w].count;
				CVUInteger* PData = P[w].data;
				for (p=0; p<PCount; p++) {
					v = PData[p];
					delta[v] += sigma[v]/(double)sigma[w]*(1.0+delta[w]);
				}
				if(w!=s){
					localCentrality[w] += sWeight*delta[w];
				}
			}
		}
		for (i=0; i<verticesCount; i++){
			CVUIntegerArrayDestroy(P+i);
		}
		free(P);
		free(sigma);
		free(d);
		free(delta);
		CVIntegerArrayDestroy(&S);
		CVQueueDestroy(&Q);
		
		CVParallelLoopCriticalRegionStart(centralityLoop){
			CVFloat* centralityData = centrality->data;
			for (i=0; i<verticesCount; i++) {
				centralityData[i] += localCentrality[i];
			}
		}
		CVParallelLoopCriticalRegionEnd(centralityLoop);
		free(localCentrality);
		
	}CVParallelForEnd(centralityLoop);
	
	return CVTrue;
}

#endif //CV_ENABLE_PARALLELISM



CVBool CVNetworkCalculateCentrality_implementation(const CVNetwork* network,CVFloatArray* centrality, CVOperationControl* operationControl){
	CVSize verticesCount = network->verticesCount;
	
	CVInteger* currentProgress = NULL;
	void (*updateCallback)(CVOperationControl*)  = NULL;
	if(operationControl){
		operationControl->maxProgress = verticesCount;
		operationControl->currentProgress = 0;
		currentProgress = &(operationControl->currentProgress);
		updateCallback = operationControl->updateCallback;
	}
	
	CVFloatArrayReallocToCapacity(verticesCount, centrality);
	CVFloatArraySetCount(verticesCount, centrality);
		//memset(centrality->data, 0, verticesCount*sizeof(CVFloat));
	CVIndex i;
	for (i=0; i<verticesCount; i++) {
		if(network->verticesEnabled[i]){
			centrality->data[i]  = 0.0;
		}
	}
	
	const CVFloat* verticesWeights = network->verticesWeights;
	const CVBool* verticesEnabled = network->verticesEnabled;
	
	CVFloat* centralityData = centrality->data;
	
	CVUIntegerArray* P = calloc(verticesCount, sizeof(CVUIntegerArray));
	CVInteger* sigma = calloc(verticesCount, sizeof(CVInteger));
	CVInteger* d = calloc(verticesCount, sizeof(CVInteger));
	double* delta = calloc(verticesCount, sizeof(double));
	for (i=0; i<verticesCount; i++) {
		CVUIntegerArray newArray;
		CVUIntegerArrayInitWithCapacity(1, &newArray);
		P[i] = newArray;
	}
	CVIndex s;
	CVIntegerStack S = CVIntegerStackMake();
	CVQueue Q = CVQueueCreate();
	
	for (s=0; s<verticesCount; s++) {
		if(currentProgress){
			CVAtomicIncrementInteger(currentProgress);
			if(updateCallback){
				updateCallback(operationControl);
			}
		}
		if(CVUnlikely(!verticesEnabled[s])){
			continue;
		}
		CVFloat sWeight = verticesWeights[s];
			//printf("%ld/%ld\n",s,verticesCount);
		S.count = 0;
		for (i=0; i<verticesCount; i++) {
			P[i].count=0;
			d[i] = CVIntegerMAX;
			sigma[i] = 0;
			delta[i] = 0;
		}
		sigma[s] = 1;
		d[s] = 0;
		CVInteger v;
		CVQueuePush(&Q, s);
		while (CVQueueDequeue(&Q,&v)) {
			CVIntegerStackPush(v, &S);
			CVSize vEdgesCount = network->vertexNumOfEdges[v];
			CVIndex* vNeighbor = network->vertexEdgesLists[v];
			CVIndex w,e;
			for (e=0; e<vEdgesCount; e++) {
				w = vNeighbor[e];
				if(CVLikely(verticesEnabled[w])){
					if(d[w]==CVIntegerMAX){
						d[w]= d[v] + 1;// FIXME: Change to the w-v weight;
						CVQueuePush(&Q, w);
					}
					if(d[w] == d[v] + 1){
						sigma[w]+= sigma[v];
						CVUIntegerArrayAdd(v, P+w);
						
					}
				}
			}
		}
		while(S.count>0){
			CVIndex w = CVIntegerStackPop(&S);
			CVIndex v,p;
			CVSize PCount = P[w].count;
			CVUInteger* PData = P[w].data;
			for (p=0; p<PCount; p++) {
				v = PData[p];
				delta[v] += sigma[v]/(double)sigma[w]*(1.0+delta[w]);
			}
			if(w!=s){
				centralityData[w] += sWeight*delta[w];
			}
		}
	}
	CVQueueDestroy(&Q);
	CVIntegerArrayDestroy(&S);
	free(P);
	free(sigma);
	free(d);
	free(delta);
	
	
	/*
	 printf("----\nFinished\n-----\n");
	 CVSize N = network->verticesCount;
	 for (i=0; i<verticesCount; i++) {
	 printf("%g\n",centralityData[i]/(N-1)/(N-2));
	 }
	 */
	return CVTrue;
}


CVBool CVNetworkCalculateCentrality(const CVNetwork* network,CVFloatArray* centrality, CVOperationControl* operationControl){
	CV_BenchmarkPrepare(CVNetworkCalculateCentrality);
	CV_BenchmarkStart(CVNetworkCalculateCentrality);
	
	CVBool returnValue;
	
#if CV_ENABLE_PARALLELISM
	CVInteger maxParallelBlocksCount = kCVDefaultParallelBlocks;
	CVSize problemSize = network->verticesCount;
	if(operationControl){
		maxParallelBlocksCount = operationControl->maxParallelBlocks;
	}
	
	if(network->edgeWeighted){
		if(network&&problemSize>=128&&maxParallelBlocksCount>1){
			returnValue = CVNetworkCalculateCentrality_weighted_parallel_implementation(network, centrality, operationControl);
		}else{
			returnValue = CVNetworkCalculateCentrality_weighted_implementation(network, centrality, operationControl);
		}
	}else{
		if(network&&problemSize>=128&&maxParallelBlocksCount>1){
			returnValue = CVNetworkCalculateCentrality_parallel_implementation(network, centrality, operationControl);
		}else{
			returnValue = CVNetworkCalculateCentrality_implementation(network, centrality, operationControl);
		}
	}
#else
	
	if(network->edgeWeighted){
		returnValue = CVNetworkCalculateCentrality_weighted_implementation(network, centrality, operationControl);
	}else{
		returnValue = CVNetworkCalculateCentrality_implementation(network, centrality, operationControl);
	}
#endif //CV_ENABLE_PARALLELISM
	
	CV_BenchmarkStop(CVNetworkCalculateCentrality);
	CV_BenchmarkPrint(CVNetworkCalculateCentrality);
	return returnValue;
}



























//////Stress




#if CV_ENABLE_PARALLELISM
CVBool CVNetworkCalculateStressCentrality_parallel_implementation(const CVNetwork* network,CVFloatArray* centrality, CVOperationControl* operationControl){
	CVSize verticesCount = network->verticesCount;
	
	CVSize unrolledLoops = kCVDefaultParallelBlocks;
	CVInteger* currentProgress = NULL;
	void (*updateCallback)(CVOperationControl*)  = NULL;
	if(operationControl){
		operationControl->maxProgress = verticesCount;
		operationControl->currentProgress = 0;
		currentProgress = &(operationControl->currentProgress);
		if(operationControl->maxParallelBlocks>0){
			unrolledLoops = operationControl->maxParallelBlocks;
		}
		updateCallback = operationControl->updateCallback;
	}
	
	
	CVFloatArrayReallocToCapacity(verticesCount, centrality);
	CVFloatArraySetCount(verticesCount, centrality);
	//memset(centrality->data, 0, verticesCount*sizeof(CVFloat));
	CVIndex i;
	for (i=0; i<verticesCount; i++) {
		if(network->verticesEnabled[i]){
			centrality->data[i]  = 0.0f;
		}
	}
	CVSize unrolledSize = 1 + ((verticesCount - 1) / unrolledLoops);
	
	CVParallelForStart(centralityLoop, blockIndex, unrolledLoops){
		
		
		const CVFloat* verticesWeights = network->verticesWeights;
		const CVBool* verticesEnabled = network->verticesEnabled;
		
		CVUIntegerArray* P = calloc(verticesCount, sizeof(CVUIntegerArray));
		CVInteger* sigma = calloc(verticesCount, sizeof(CVInteger));
		CVInteger* d = calloc(verticesCount, sizeof(CVInteger));
		double* delta = calloc(verticesCount, sizeof(double));
		
		double* localCentrality = calloc(verticesCount, sizeof(double));
		
		CVIndex i;
		for (i=0; i<verticesCount; i++) {
			CVUIntegerArray newArray;
			CVUIntegerArrayInitWithCapacity(1, &newArray);
			P[i] = newArray;
		}
		
		CVIndex s;
		CVIntegerStack S = CVIntegerStackMake();
		CVQueue Q = CVQueueCreate();
		
		CVSize maxIndex = CVMIN((blockIndex+1)*unrolledSize, verticesCount);
		for(s=blockIndex*unrolledSize;s<maxIndex;s++){
			if(currentProgress){
				CVAtomicIncrementInteger(currentProgress);
				if(updateCallback){
					updateCallback(operationControl);
				}
			}
			if(CVUnlikely(!verticesEnabled[s])){
				continue;
			}
			CVFloat sWeight = verticesWeights[s];
			S.count = 0;
			for (i=0; i<verticesCount; i++) {
				P[i].count=0;
				d[i] = CVIntegerMAX;
				sigma[i] = 0;
				delta[i] = 0;
			}
			sigma[s] = 1;
			d[s] = 0;
			CVInteger v;
			CVQueuePush(&Q, s);
			while (CVQueueDequeue(&Q,&v)) {
				CVIntegerStackPush(v, &S);
				CVSize vEdgesCount = network->vertexNumOfEdges[v];
				CVIndex* vNeighbor = network->vertexEdgesLists[v];
				CVIndex w,e;
				for (e=0; e<vEdgesCount; e++) {
					w = vNeighbor[e];
					if(CVLikely(verticesEnabled[w])){
						if(d[w]==CVIntegerMAX){
							d[w] = d[v] + 1;// FIXME: Change to the w-v weight;
							//fh_insertkey(Q, d[v], w);
							//if(d[w]<=2){
							CVQueuePush(&Q, w);
							//}
						}
						if(d[w]== d[v]+1){
							sigma[w]+= sigma[v];
							CVUIntegerArrayAdd(v, &(P[w]));
						}
					}
				}
			}
			while(S.count>0){
				CVIndex w = CVIntegerStackPop(&S);
				CVIndex v,p;
				CVSize PCount = P[w].count;
				CVUInteger* PData = P[w].data;
				for (p=0; p<PCount; p++) {
					v = PData[p];
					delta[v] += (1.0+delta[w]);
				}
				if(w!=s){
					localCentrality[w] += sigma[w]*sWeight*delta[w];
				}
			}
		}
		for (i=0; i<verticesCount; i++){
			CVUIntegerArrayDestroy(P+i);
		}
		free(P);
		free(sigma);
		free(d);
		free(delta);
		CVIntegerArrayDestroy(&S);
		CVQueueDestroy(&Q);
		
		CVParallelLoopCriticalRegionStart(centralityLoop){
			CVFloat* centralityData = centrality->data;
			for (i=0; i<verticesCount; i++) {
				centralityData[i] += localCentrality[i];
			}
		}
		CVParallelLoopCriticalRegionEnd(centralityLoop);
		free(localCentrality);
		
	}CVParallelForEnd(centralityLoop);
	
	return CVTrue;
}

#endif //CV_ENABLE_PARALLELISM



CVBool CVNetworkCalculateStressCentrality_implementation(const CVNetwork* network,CVFloatArray* centrality, CVOperationControl* operationControl){
	CVSize verticesCount = network->verticesCount;
	
	CVInteger* currentProgress = NULL;
	void (*updateCallback)(CVOperationControl*)  = NULL;
	if(operationControl){
		operationControl->maxProgress = verticesCount;
		operationControl->currentProgress = 0;
		currentProgress = &(operationControl->currentProgress);
		updateCallback = operationControl->updateCallback;
	}
	
	CVFloatArrayReallocToCapacity(verticesCount, centrality);
	CVFloatArraySetCount(verticesCount, centrality);
	//memset(centrality->data, 0, verticesCount*sizeof(CVFloat));
	CVIndex i;
	for (i=0; i<verticesCount; i++) {
		if(network->verticesEnabled[i]){
			centrality->data[i]  = 0.0;
		}
	}
	
	const CVFloat* verticesWeights = network->verticesWeights;
	const CVBool* verticesEnabled = network->verticesEnabled;
	
	CVFloat* centralityData = centrality->data;
	
	CVUIntegerArray* P = calloc(verticesCount, sizeof(CVUIntegerArray));
	CVInteger* sigma = calloc(verticesCount, sizeof(CVInteger));
	CVInteger* d = calloc(verticesCount, sizeof(CVInteger));
	double* delta = calloc(verticesCount, sizeof(double));
	for (i=0; i<verticesCount; i++) {
		CVUIntegerArray newArray;
		CVUIntegerArrayInitWithCapacity(1, &newArray);
		P[i] = newArray;
	}
	CVIndex s;
	CVIntegerStack S = CVIntegerStackMake();
	CVQueue Q = CVQueueCreate();
	
	for (s=0; s<verticesCount; s++) {
		if(currentProgress){
			CVAtomicIncrementInteger(currentProgress);
			if(updateCallback){
				updateCallback(operationControl);
			}
		}
		if(CVUnlikely(!verticesEnabled[s])){
			continue;
		}
		CVFloat sWeight = verticesWeights[s];
		//printf("%ld/%ld\n",s,verticesCount);
		S.count = 0;
		for (i=0; i<verticesCount; i++) {
			P[i].count=0;
			d[i] = CVIntegerMAX;
			sigma[i] = 0;
			delta[i] = 0;
		}
		sigma[s] = 1;
		d[s] = 0;
		CVInteger v;
		CVQueuePush(&Q, s);
		while (CVQueueDequeue(&Q,&v)) {
			CVIntegerStackPush(v, &S);
			CVSize vEdgesCount = network->vertexNumOfEdges[v];
			CVIndex* vNeighbor = network->vertexEdgesLists[v];
			CVIndex w,e;
			for (e=0; e<vEdgesCount; e++) {
				w = vNeighbor[e];
				if(CVLikely(verticesEnabled[w])){
					if(d[w]==CVIntegerMAX){
						d[w]= d[v] + 1;// FIXME: Change to the w-v weight;
						CVQueuePush(&Q, w);
					}
					if(d[w] == d[v] + 1){
						sigma[w]+= sigma[v];
						CVUIntegerArrayAdd(v, P+w);
						
					}
				}
			}
		}
		while(S.count>0){
			CVIndex w = CVIntegerStackPop(&S);
			CVIndex v,p;
			CVSize PCount = P[w].count;
			CVUInteger* PData = P[w].data;
			for (p=0; p<PCount; p++) {
				v = PData[p];
				delta[v] += (1.0+delta[w]);
			}
			if(w!=s){
				centralityData[w] += sigma[w]*sWeight*delta[w];
			}
		}
	}
	CVQueueDestroy(&Q);
	CVIntegerArrayDestroy(&S);
	free(P);
	free(sigma);
	free(d);
	free(delta);
	
	
	/*
	 printf("----\nFinished\n-----\n");
	 CVSize N = network->verticesCount;
	 for (i=0; i<verticesCount; i++) {
	 printf("%g\n",centralityData[i]/(N-1)/(N-2));
	 }
	 */
	return CVTrue;
}


CVBool CVNetworkCalculateStressCentrality(const CVNetwork* network,CVFloatArray* centrality, CVOperationControl* operationControl){
	CV_BenchmarkPrepare(CVNetworkCalculateStressCentrality);
	CV_BenchmarkStart(CVNetworkCalculateStressCentrality);
	
	CVBool returnValue;
	
#if CV_ENABLE_PARALLELISM
	CVInteger maxParallelBlocksCount = kCVDefaultParallelBlocks;
	CVSize problemSize = network->verticesCount;
	if(operationControl){
		maxParallelBlocksCount = operationControl->maxParallelBlocks;
	}
	
	if(network&&problemSize>=128&&maxParallelBlocksCount>1){
		returnValue = CVNetworkCalculateStressCentrality_parallel_implementation(network, centrality, operationControl);
	}else{
		returnValue = CVNetworkCalculateStressCentrality_implementation(network, centrality, operationControl);
	}
#else
	returnValue = CVNetworkCalculateStressCentrality_implementation(network, centrality, operationControl);
#endif //CV_ENABLE_PARALLELISM
	
	CV_BenchmarkStop(CVNetworkCalculateStressCentrality);
	CV_BenchmarkPrint(CVNetworkCalculateStressCentrality);
	return returnValue;
}

