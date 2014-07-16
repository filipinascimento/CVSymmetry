//
//  CVNetworkCentrality.h
//  CVNetwork
//
//  Created by Filipi Nascimento Silva on 8/27/13.
//  Copyright (c) 2013 Filipi Nascimento Silva. All rights reserved.
//

#ifndef CVNetwork_CVNetworkCentrality_h
#define CVNetwork_CVNetworkCentrality_h

#include "CVNetwork.h"
#include "CVBasicArrays.h"
#include "CVSimpleQueue.h"

CVBool CVNetworkCalculateCentrality(const CVNetwork* network,CVFloatArray* centrality, CVOperationControl* operationControl);
CVBool CVNetworkCalculateStressCentrality(const CVNetwork* network,CVFloatArray* centrality, CVOperationControl* operationControl);


#endif
