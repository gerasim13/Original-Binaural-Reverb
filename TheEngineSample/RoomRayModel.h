//
//  RoomRayModel.h
//  TheEngineSample
//
//  Created by Hans on 6/11/15.
//  Copyright © 2015 A Tasty Pixel. All rights reserved.
//

#ifndef RoomRayModel_h
#define RoomRayModel_h

#include <stdio.h>
#include "Vector3D.hpp"
#include "mactypes.h"

#define RRM_MAX_CORNERS 100

class RoomRayModel {
private:
    Vector3D corners[RRM_MAX_CORNERS];
    Vector3D wallOrientations[RRM_MAX_CORNERS];
    float wallLengths[RRM_MAX_CORNERS];
    size_t numCorners;
    float totalWallLength;
    void setBouncePoints(Vector3D* bouncePoints, Vector3D wallOrientation, Vector3D wallStart, float wallLength, size_t numPoints, float* outputGains2, float* inputGains2, Vector3D* BP);
    Vector3D getBP(float pointSpacing, Vector3D wallStart, size_t i, Vector3D wallOrientation, float randFlt);
        


    
    //original integration method
    float getGain(Vector3D start, Vector3D end, Vector3D loc);
    float integrate(Vector3D start, Vector3D end, float t, Vector3D loc, Vector3D vn);
    
    //Simple integration method
    float integrationSimple(Vector3D loc, float x, bool listLoc);
    Vector3D align(Vector3D point, Vector3D wallvector);
    float xAlignedIntegration(Vector3D loc, Vector3D ptStart, Vector3D ptEnd, bool listLoc);
    void setFloorBouncePointsGain(Vector3D* bouncePoints, float* inputGain, float* outputGain, size_t floorTaps);
    void gridBP(Vector3D* floorBouncePoints, size_t floorTaps);
    float pythagorasGain(Vector3D loc, Vector3D* bouncePoint, float height);
        Vector3D soundSourceLoc; Vector3D listenerLoc;   
    size_t floorTapsPerDimension;
    
    float getMaxGain(float xLower, float xUpper, float yLower, float yUpper);
    float calcMaxGain(float x, float y);
    float maxFloorGain = 1.55463f;
    float gridArea;
    
    float numWallPoints;
public:
    RoomRayModel();

    
    void setRoomGeometry(Vector3D* corners, size_t numCorners);
    
    void setLocation(float* rayLengths,size_t numTaps, Vector3D listenerLocation, Vector3D soundSourceLocation, Vector3D* bouncePoints, float* outputGains2, float* inputGains2, size_t floorTaps, Vector3D* BP);
    
    void setRayTracingPoints(Vector3D* bouncePoints, Vector3D ssLocation, float rheight, float rwidth, int numpoints, float* outputGains, float* inputGains,Vector3D listloc);
//    void setRayTracingPoints(Vector3D* bouncePoints, Vector3D ssLoc, float rheight, float rwidth, int numpoints,Vector3D listloc);

};

#endif /* RoomRayModel_h */
