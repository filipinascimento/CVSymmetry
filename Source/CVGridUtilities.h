//
//  CVGridUtilities.h
//  CVNetwork
//
//  Created by Filipi Nascimento Silva on 28/05/14.
//  Copyright (c) 2014 Filipi Nascimento Silva. All rights reserved.
//

#ifndef CVNetwork_CVGridUtilities_h
#define CVNetwork_CVGridUtilities_h

#include "CVCommons.h"

CV_INLINE CVFloat CVGridCalcIntegerDistance(const CVInteger* origin, const CVInteger* coordVector,const CVUInteger dimensions){
	switch (dimensions) {
		case 1:{
			return CVABS(origin[0]-coordVector[0]);
		}
		case 2:{
			CVInteger dx = origin[0]-coordVector[0];
			CVInteger dy = origin[1]-coordVector[1];
			return sqrtf(dx*dx + dy*dy);
		}
		case 3:{
			CVInteger dx = origin[0]-coordVector[0];
			CVInteger dy = origin[1]-coordVector[1];
			CVInteger dz = origin[2]-coordVector[2];
			return sqrtf(dx*dx + dy*dy + dz*dz);
		}
		default:{
			CVInteger sqrSum = 0;
			CVIndex d;
			for(d=0;d<dimensions;d++){
				CVInteger dDiff = coordVector[d]-origin[d];
				sqrSum+= dDiff*dDiff;
			}
			return sqrtf(sqrSum);
		}
	}
}


CV_INLINE CVFloat CVGridCalcIntegerDistanceFromOrigin(const CVInteger* coordVector, CVInteger dimensions){
	switch (dimensions) {
		case 1:{
			return CVABS(coordVector[0]);
		}
		case 2:{
			CVInteger dx = coordVector[0];
			CVInteger dy = coordVector[1];
			return sqrtf(dx*dx + dy*dy);
		}
		case 3:{
			CVInteger dx = coordVector[0];
			CVInteger dy = coordVector[1];
			CVInteger dz = coordVector[2];
			return sqrtf(dx*dx + dy*dy + dz*dz);
		}
		default:{
			CVInteger sqrSum = 0;
			CVIndex d;
			for(d=0;d<dimensions;d++){
				CVInteger dDiff = coordVector[d];
				sqrSum+= dDiff*dDiff;
			}
			return sqrtf(sqrSum);
		}
	}
}
CV_INLINE CVIndex CVGridLinearIndexFromCoordinates(const CVInteger* coordinates, const CVSize* gridSize,CVUInteger dimensions){
	switch (dimensions) {
		case 1:{
			return coordinates[0];
		}
		case 2:{
			return coordinates[0] + gridSize[0]*coordinates[1];
		}
		case 3:{
			return coordinates[0] + gridSize[0]*(coordinates[1]+gridSize[1]*(coordinates[2]));
		}
		case 4:{
			return coordinates[0] + gridSize[0]*(coordinates[1]+gridSize[1]*(coordinates[2]+gridSize[3]*coordinates[3]));
		}
		default:{
			CVInteger linearIndex = 0;
			CVInteger coefficient = 1;
			CVIndex curDim;
			for(curDim=0;curDim<dimensions;curDim++){
				linearIndex += coefficient*coordinates[curDim];
				coefficient*=gridSize[curDim];
			}
			return linearIndex;
		}
	}
}

CV_INLINE void CVGridGetCoordinatesFromLinearIndex(CVInteger linearIndex, const CVSize* gridSize,CVUInteger dimensions,CVInteger* coordinates){
	switch (dimensions) {
		case 1:{
			coordinates[0] = linearIndex;
			return;
		}
		case 2:{
			coordinates[0] = linearIndex % gridSize[0];
			coordinates[1] = (linearIndex/gridSize[0]) % gridSize[1];
			return;
		}
		case 3:{
			coordinates[0] = linearIndex % gridSize[0];
			linearIndex /= (gridSize[0]);
			coordinates[1] = linearIndex % gridSize[1];
			linearIndex /= (gridSize[1]);
			coordinates[2] = linearIndex % gridSize[2];
			return;
		}
		case 4:{
			coordinates[0] = linearIndex % gridSize[0];
			linearIndex /= (gridSize[0]);
			coordinates[1] = linearIndex % gridSize[1];
			linearIndex /= (gridSize[1]);
			coordinates[2] = linearIndex % gridSize[2];
			linearIndex /= (gridSize[2]);
			coordinates[3] = linearIndex % gridSize[3];
			return;
		}
		default:{
			CVIndex curDim;
			for(curDim=0;curDim<dimensions;curDim++){
				coordinates[curDim] = linearIndex % gridSize[curDim];
				linearIndex /= (gridSize[curDim]);
			}
		}
	}
}


CV_INLINE CVBool CVGridCheckCoordinateBounds(CVInteger coordinateValue, CVSize bound){
	if(coordinateValue<0 || coordinateValue>=bound){
		return CVFalse;
	}else{
		return CVTrue;
	}
}


CV_INLINE CVBool CVGridGetDisplacedCoordinate(const CVInteger* originalCoordinates, const  CVInteger* displacement,  const  CVInteger* center,const  CVSize* gridSize, CVUInteger dimensions, CVInteger* dstCoordinates){
	switch (dimensions) {
		case 1:{
			CVInteger newCoordinate = originalCoordinates[0] + displacement[0] - center[0];
			if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[0])){
				dstCoordinates[0] = newCoordinate;
			}else{
				return CVFalse;
			}
			break;
		}
		case 2:{
			CVInteger newCoordinate = originalCoordinates[0] + displacement[0] - center[0];
			if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[0])){
				dstCoordinates[0] = newCoordinate;
			}else{
				return CVFalse;
			}
			newCoordinate = originalCoordinates[1] + displacement[1] - center[1];
			if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[1])){
				dstCoordinates[1] = newCoordinate;
			}else{
				return CVFalse;
			}
			break;
		}
		case 3:{
			CVInteger newCoordinate = originalCoordinates[0] + displacement[0] - center[0];
			if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[0])){
				dstCoordinates[0] = newCoordinate;
			}else{
				return CVFalse;
			}
			newCoordinate = originalCoordinates[1] + displacement[1] - center[1];
			if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[1])){
				dstCoordinates[1] = newCoordinate;
			}else{
				return CVFalse;
			}
			newCoordinate = originalCoordinates[2] + displacement[2] - center[2];
			if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[2])){
				dstCoordinates[2] = newCoordinate;
			}else{
				return CVFalse;
			}
			break;
		}
		case 4:{
			CVInteger newCoordinate = originalCoordinates[0] + displacement[0] - center[0];
			if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[0])){
				dstCoordinates[0] = newCoordinate;
			}else{
				return CVFalse;
			}
			newCoordinate = originalCoordinates[1] + displacement[1] - center[1];
			if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[1])){
				dstCoordinates[1] = newCoordinate;
			}else{
				return CVFalse;
			}
			newCoordinate = originalCoordinates[2] + displacement[2] - center[2];
			if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[2])){
				dstCoordinates[2] = newCoordinate;
			}else{
				return CVFalse;
			}
			newCoordinate = originalCoordinates[3] + displacement[3] - center[3];
			if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[3])){
				dstCoordinates[3] = newCoordinate;
			}else{
				return CVFalse;
			}
			break;
		}
		default:{
			CVIndex i;
			for(i=0;i<dimensions;i++){
				CVInteger newCoordinate =originalCoordinates[i]+displacement[i] - center[i];
				if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[i])){
					dstCoordinates[i] = newCoordinate;
				}else{
					return CVFalse;
				}
			}
			break;
		}
	}
	return CVTrue;
}



CV_INLINE CVBool CVGridGetDisplacedToroidalCoordinate(const CVInteger* originalCoordinates, const  CVInteger* displacement,  const  CVInteger* center,const  CVSize* gridSize, CVUInteger dimensions, CVInteger* dstCoordinates){
	switch (dimensions) {
		case 1:{
			CVInteger newCoordinate = (gridSize[0] + originalCoordinates[0] + displacement[0] - center[0]) % gridSize[0];
			if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[0])){
				dstCoordinates[0] = newCoordinate;
			}else{
				return CVFalse;
			}
			break;
		}
		case 2:{
			CVInteger newCoordinate = (gridSize[0] + originalCoordinates[0] + displacement[0] - center[0]) % gridSize[0];
			if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[0])){
				dstCoordinates[0] = newCoordinate;
			}else{
				return CVFalse;
			}
			newCoordinate = (gridSize[1] + originalCoordinates[1] + displacement[1] - center[1]) % gridSize[1];
			if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[1])){
				dstCoordinates[1] = newCoordinate;
			}else{
				return CVFalse;
			}
			break;
		}
		case 3:{
			CVInteger newCoordinate = (gridSize[0] + originalCoordinates[0] + displacement[0] - center[0]) % gridSize[0];
			if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[0])){
				dstCoordinates[0] = newCoordinate;
			}else{
				return CVFalse;
			}
			newCoordinate = (gridSize[1] + originalCoordinates[1] + displacement[1] - center[1]) % gridSize[1];
			if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[1])){
				dstCoordinates[1] = newCoordinate;
			}else{
				return CVFalse;
			}
			newCoordinate = (gridSize[2] + originalCoordinates[2] + displacement[2] - center[2]) % gridSize[2];
			if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[2])){
				dstCoordinates[2] = newCoordinate;
			}else{
				return CVFalse;
			}
			break;
		}
		case 4:{
			CVInteger newCoordinate = (gridSize[0] + originalCoordinates[0] + displacement[0] - center[0]) % gridSize[0];
			if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[0])){
				dstCoordinates[0] = newCoordinate;
			}else{
				return CVFalse;
			}
			newCoordinate = (gridSize[1] + originalCoordinates[1] + displacement[1] - center[1]) % gridSize[1];
			if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[1])){
				dstCoordinates[1] = newCoordinate;
			}else{
				return CVFalse;
			}
			newCoordinate = (gridSize[2] + originalCoordinates[2] + displacement[2] - center[2]) % gridSize[2];
			if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[2])){
				dstCoordinates[2] = newCoordinate;
			}else{
				return CVFalse;
			}
			newCoordinate = (gridSize[3] + originalCoordinates[3] + displacement[3] - center[3]) % gridSize[3];
			if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[3])){
				dstCoordinates[3] = newCoordinate;
			}else{
				return CVFalse;
			}
			break;
		}
		default:{
			CVIndex i;
			for(i=0;i<dimensions;i++){
				CVInteger newCoordinate = (gridSize[i] + originalCoordinates[i]+displacement[i] - center[i]) % gridSize[i];
				if(CVGridCheckCoordinateBounds(newCoordinate, gridSize[i])){
					dstCoordinates[i] = newCoordinate;
				}else{
					return CVFalse;
				}
			}
			break;
		}
	}
	return CVTrue;
}


#endif
