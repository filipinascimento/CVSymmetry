//
//  CVSimpleCVQueue.h
//  CVNetwork
//
//  Created by Filipi Nascimento Silva on 11/23/12.
//  Copyright (c) 2012 Filipi Nascimento Silva. All rights reserved.
//

#ifndef CVNetwork_CVSimpleCVQueue_h
#define CVNetwork_CVSimpleCVQueue_h
#include "CVCommons.h"
/**
 * This sample is about how to implement a queue in c
 *
 * Type of item is CVInteger
 * Add item to tail
 * Get item from head
 * Can get the size
 * Can display all content
 */
/**
 * The Node struct,
 * contains item and the poCVIntegerer that poCVInteger to next node.
 */
typedef struct CVQueueNode {
	CVInteger item;
	struct CVQueueNode* next;
} CVQueueNode;
/**
 * The CVQueue struct, contains the poCVIntegerers that
 * poCVInteger to first node and last node, the size of the CVQueue,
 * and the function poCVIntegerers.
 */
typedef struct CVQueue {
	CVQueueNode* head;
	CVQueueNode* tail;
	
	void (*push) (struct CVQueue*, CVInteger); // add item to tail
	// get item from head and remove it from queue
	CVInteger (*pop) (struct CVQueue*);
	// get item from head but keep it in queue
	CVInteger (*peek) (struct CVQueue*);
	// size of this queue
	CVInteger size;
	CVQueueNode** nodesPool;
	CVInteger nodesPoolCount;
	CVInteger nodesPoolCapacity;
	
} CVQueue;
/**
 * Push an item CVIntegero queue, if this is the first item,
 * both queue->head and queue->tail will poCVInteger to it,
 * otherwise the oldtail->next and tail will poCVInteger to it.
 */
void CVQueuePush (CVQueue* queue, CVInteger item);
/**
 * Return and remove the first item.
 */
CVInteger CVQueuePop (CVQueue* queue);
/**
 * Return but not remove the first item.
 */
CVInteger CVQueuePeek (CVQueue* queue);
/**
 * Create and initiate a CVQueue
 */
CVQueue CVQueueCreate ();
CVBool CVQueueDequeue (CVQueue* queue, CVInteger *value);
void CVQueueDestroy (CVQueue* queue);


#endif