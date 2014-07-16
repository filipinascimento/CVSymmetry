//
//  CVSimpleCVQueue.m
//  CVNetwork
//
//  Created by Filipi Nascimento Silva on 11/23/12.
//  Copyright (c) 2012 Filipi Nascimento Silva. All rights reserved.
//

#include "CVSimpleQueue.h"


CV_INLINE CVQueueNode* getNodeFromPool(CVQueue* queue){
	if(queue->nodesPoolCount>0){
		queue->nodesPoolCount--;
		return queue->nodesPool[queue->nodesPoolCount];
	}else{
		return (CVQueueNode*)malloc(sizeof(CVQueueNode));
	}
}

CV_INLINE void putNodeInPool(CVQueue* queue,CVQueueNode* node){
	queue->nodesPoolCount++;
	if(CVUnlikely(queue->nodesPoolCapacity < queue->nodesPoolCount)){
		queue->nodesPoolCapacity = CVCapacityGrow(queue->nodesPoolCount);
		queue->nodesPool = realloc(queue->nodesPool, sizeof(CVQueueNode*)*queue->nodesPoolCapacity);
	}
	queue->nodesPool[queue->nodesPoolCount-1] = node;
}


/**
 * Push an item CVIntegero queue, if this is the first item,
 * both queue->head and queue->tail will poCVInteger to it,
 * otherwise the oldtail->next and tail will poCVInteger to it.
 */
void CVQueuePush (CVQueue* queue, CVInteger item) {
	// Create a new node
	CVQueueNode* n = getNodeFromPool(queue);
	n->item = item;
	n->next = NULL;
	
	if (queue->head == NULL) { // no head
		queue->head = n;
	} else{
		queue->tail->next = n;
	}
	queue->tail = n;
	queue->size++;
}
/**
 * Return and remove the first item.
 */
CVInteger CVQueuePop (CVQueue* queue) {
    // get the first item
	CVQueueNode* head = queue->head;
	CVInteger item = head->item;
	// move head poCVIntegerer to next node, decrease size
	queue->head = head->next;
	queue->size--;
	// free the memory of original head
	free(head);
	return item;
}


/**
 * Return and remove the first item.
 */
CVBool CVQueueDequeue (CVQueue* queue,CVInteger *value) {
    // get the first item
	if (queue->size>0) {
		CVQueueNode* head = queue->head;
		CVInteger item = head->item;
		// move head poCVIntegerer to next node, decrease size
		queue->head = head->next;
		queue->size--;
		// free the memory of original head
		putNodeInPool(queue,head);
		*value = item;
		return CVTrue;
	}else{
		return CVFalse;
	}
}

/**
 * Return but not remove the first item.
 */
CVInteger CVQueuePeek (CVQueue* queue) {
	CVQueueNode* head = queue->head;
	return head->item;
}


/**
 * Create and initiate a CVQueue
 */
CVQueue CVQueueCreate () {
	CVQueue queue;
	queue.size = 0;
	queue.head = NULL;
	queue.tail = NULL;
	queue.push = &CVQueuePush;
	queue.pop = &CVQueuePop;
	queue.peek = &CVQueuePeek;
	queue.nodesPool = calloc(2, sizeof(CVQueueNode*));
	queue.nodesPoolCapacity = 2;
	queue.nodesPoolCount = 0;
	return queue;
}

void CVQueueDestroy (CVQueue* queue) {
	CVInteger value;
	while(CVQueueDequeue(queue,&value)){};
	CVIndex i;
	for (i=0; i<queue->nodesPoolCount; i++) {
		free(queue->nodesPool[i]);
	}
	free(queue->nodesPool);
}
