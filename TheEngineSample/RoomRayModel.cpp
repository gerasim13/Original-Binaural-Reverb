//
//  RoomRayModel.c
//  TheEngineSample
//
//  Created by Hans on 6/11/15.
//  Copyright © 2015 A Tasty Pixel. All rights reserved.
//

#include "RoomRayModel.h"
#include "assert.h"
#include "string.h"
#include "math.h"
#include "FDN.h"

RoomRayModel::RoomRayModel(){
    numCorners = 0;
}


//Visibility check when setting a bounce point on the wall
void RoomRayModel::setBouncePoints(Point2d* bouncePoints, Point2d wallOrientation, Point2d wallStart, float wallLength, size_t numPoints, float* outputGains, float* inputGains){
    
    // average space between each pair of points
    float pointSpacing = wallLength / numPoints;
    
    Point2d prevStart = wallStart;

    // set the points at even but randomly jittered locations
    for (size_t i = 0; i < numPoints; i++) {
        float randFlt = (float)rand() / (float)RAND_MAX;
        

        bouncePoints[i] = getBP(pointSpacing, wallStart, i, wallOrientation, randFlt);
        if(i>0){
            Point2d start = prevStart;
            Point2d difference = (bouncePoints[i] - bouncePoints[i-1]).scalarMul(0.5f);
            Point2d end = bouncePoints[i-1] + difference;

            outputGains[i-1] =  sqrtf( xAlignedIntegration(listenerLoc, start, end, true));
            inputGains[i-1] = sqrtf(xAlignedIntegration(soundSourceLoc, start, end, false));
            prevStart = Point2d(end.x, end.y);
        }
        
    }
    
    //do the last gain
    Point2d end = wallStart + wallOrientation.scalarMul(wallLength);
    outputGains[numPoints-1] = sqrtf(xAlignedIntegration(listenerLoc, prevStart, end, true));
    inputGains[numPoints-1] =  sqrtf(xAlignedIntegration(soundSourceLoc, prevStart, end, false));
    
}

Point2d RoomRayModel::getBP(float pointSpacing, Point2d wallStart, size_t i, Point2d wallOrientation, float randFlt){
    float distance = (((float)i+randFlt) * pointSpacing);
    Point2d bp = wallStart + wallOrientation.scalarMul(distance);
    return bp;
}


void RoomRayModel::gridBP(Point2d* floorBouncePoints, size_t floorTaps){
    float xSpacing = wallLengths[0] / floorTaps;
    float ySpacing = wallLengths[1] / floorTaps;
    gridArea = xSpacing * ySpacing;
    
    for (size_t i = 0; i < floorTaps; i++){ //y
        for (size_t j = 0; j < floorTaps; j++) { //x
            float randFltX = (float)rand() / (float)RAND_MAX;
            float randFltY = (float)rand() / (float)RAND_MAX;
            //printf("index  %lu : ", i*floorTaps+j);
            floorBouncePoints[i*floorTaps + j] = Point2d(((float)j + randFltX) * xSpacing, ((float)i + randFltY) * ySpacing);
            //printf(" --- PT x %f y %f \n", floorBouncePoints[i*floorTaps + j].x, floorBouncePoints[i*floorTaps + j].y);
        }
    }
}

void RoomRayModel::setFloorBouncePointsGain(Point2d* bouncePoints, float* inputGain, float* outputGain, size_t floorTaps){
    for (size_t i = 0; i < floorTaps; i++){
      //  printf("GridArea is : %f \n", gridArea);
        inputGain[i] = gridArea * pythagorasGain(soundSourceLoc, &bouncePoints[i], HEIGHT);
        outputGain[i] = gridArea * pythagorasGain(listenerLoc, &bouncePoints[i], HEIGHT);
      //  printf("Floor input Gain : %f floor output Gain : %f \n", inputGain[i], outputGain[i]);
        if (inputGain[i] > maxFloorGain){
            inputGain[i] = maxFloorGain;
        }
        if (outputGain[i] > maxFloorGain){
            outputGain[i] = maxFloorGain;
        }
    }
}

float RoomRayModel::pythagorasGain(Point2d loc, Point2d* bouncePoint, float height){
    float zVal = ((float)rand()/float(RAND_MAX) * height);
    float distance = sqrtf( powf(loc.distance(*bouncePoint), 2.f) + powf(zVal, 2.f));
    bouncePoint->z = zVal;
    printf("z : %f", zVal);
    return 1.0f/distance;
}

float RoomRayModel::calcMaxGain(float x, float y){
    float a = y * logf(x + sqrtf(powf(x, 2.f)+powf(y, 2.f)));
    float b = x * (-1.f + logf(y + sqrtf(powf(x, 2.f)+powf(y, 2.f))));
    return a+b;
}
float RoomRayModel::getMaxGain(float xLower, float xUpper, float yLower, float yUpper){
    float a = calcMaxGain(yUpper, xUpper);
    float b = calcMaxGain(yUpper, xLower);
    float c = calcMaxGain(yLower, xUpper);
    float d = calcMaxGain(yLower, xLower);
    return ((a-b) - (c-d));
}
void RoomRayModel::setLocation(float* rayLengths, size_t numTaps, Point2d listenerLocation, Point2d soundSourceLocation, Point2d* bouncePoints, float* outputGains, float* inputGains, size_t floorTaps){
    srand (1);
    
//    numTaps -= floorTaps;
    floorTapsPerDimension = (size_t) sqrtf(floorTaps);
//
//    assert(numCorners > 0); // the geometry must be initialised before now
//    soundSourceLoc = soundSourceLocation;
//    listenerLoc = listenerLocation;
//    
//    // set the number of taps on each wall proportional to the
//    // length of the wall
//    size_t numTapsOnWall[RRM_MAX_CORNERS];
//    size_t totalTaps = 0;
//    for (size_t i = 0; i < RRM_MAX_CORNERS; i++) {
//        numTapsOnWall[i] = (size_t)floor(wallLengths[i]/totalWallLength * (float)numTaps);
//        totalTaps += numTapsOnWall[i];
//    }
//    
//    // if the number of taps now assigned isn't enough, add one tap to
//    // each wall until we have the desired number
//    size_t i = 0;
//    while (totalTaps < numTaps) {
//        numTapsOnWall[i]++;
//        i++;
//        totalTaps++;
//        if (i == RRM_MAX_CORNERS) i = 0;
//    }
//
//    // set bounce points for each wall
//    size_t j = 0;
//    for (size_t i = 0; i < numCorners; i++) {
//        //must be corner i-1 or shift the corner values firston
//        setBouncePoints(&bouncePoints[j], wallOrientations[i], corners[i], wallLengths[i], numTapsOnWall[i],&outputGains[j],&inputGains[j]);
//        j += numTapsOnWall[i];
//    }
//    
    
    // set bounce points for the floor
    size_t j = numTaps - floorTaps;
    gridBP(&bouncePoints[j], floorTapsPerDimension);

    
    setFloorBouncePointsGain(&bouncePoints[j], &inputGains[j], &outputGains[j], floorTaps);
    numTaps += floorTaps;
    
    
    
//    // normalize the total input gain to 1.0f
//    float totalSquaredInputGain = 0.0f;
//    for (size_t i = 0; i < numTaps; i++) {
//        inputGains[i] = fabsf(inputGains[i]*ROOMCEILING); //multiply by the room ceiling
//        totalSquaredInputGain += inputGains[i]*inputGains[i];
//    }
//
//    float inGainNormalize = 1.0f / sqrt(totalSquaredInputGain);
//    for (size_t i = 0; i < numTaps; i++) {
//        inputGains[i] *= inGainNormalize;
//       // printf("inputGains[%lu] : %f \n", i, inputGains[i]);
//    }
//    
//    //normalize the total out gain to 1.0f
//    float totalSquaredOutputGain = 0.0f;
//    for (size_t i = 0; i< numTaps; i++){
//        outputGains[i] = fabsf(outputGains[i]);
//     //   printf("Output gain: %f \n", outputGains[i]);
//        totalSquaredOutputGain += outputGains[i]*outputGains[i];
//    }
//
//    float outputGainNormalize = 1.0f / sqrtf(totalSquaredOutputGain);
//    for (size_t i = 0; i< numTaps; i++){
//        outputGains[i] *= outputGainNormalize;
//       // printf("OutputGain[%lu] : %f \n", i, outputGains[i]);
//    }
//    

}


void RoomRayModel::setRoomGeometry(Point2d* corners, size_t numCorners){
    assert(numCorners >= 3);
    
    this->numCorners = numCorners;
    
    // save the locations of the corners
    memcpy(this->corners,corners,sizeof(Point2d)*numCorners);
    
    // get normalized vectors to represent the orientation of each wall
    // and get length of each wall
    assert(numCorners < RRM_MAX_CORNERS);
    totalWallLength = 0.0f;
    for (size_t i = 1; i < numCorners; i++) {
        // get orientation vector
        wallOrientations[i] = corners[i] - corners[i-1];

        // get wall length
        wallLengths[i] = wallOrientations[i].length();
        totalWallLength += wallLengths[i];
        
        // normalize the orientation vector
        wallOrientations[i].normalize();
    }
    
    
    // get the values that wrap around from the end of the for loop above
    wallOrientations[0] = corners[0] - corners[numCorners-1];
    wallLengths[0] = wallOrientations[0].length();
    totalWallLength += wallLengths[0];
    wallOrientations[0].normalize();
    
    assert(totalWallLength > 0.0f);
    
    //change the corner indexes to match the wallOrientation indexes for setLocation method
    Point2d lastCorner = this->corners[numCorners-1];
    Point2d prevCorner = this->corners[0];
    Point2d currCorner;
    for (size_t i = 1; i<numCorners; i++){
        currCorner = this->corners[i];
        this->corners[i] = prevCorner;
        prevCorner = currCorner;
    }
    this->corners[0] = lastCorner;
}


//Simpler integration method with angle
float RoomRayModel::integrationSimple(Point2d loc, float x, bool listLoc){
    //With angle, for input gain, not listloc
    if (!listLoc){
    float a = -1.0f*loc.x + x;
    float b = sqrtf(powf(loc.x, 2.f) + pow(loc.y, 2.f) - 2.f * loc.x * x + pow(x, 2.f));
    return (a / (loc.y * b));
    }
    //Without angle, for output gain
    else{
        float a = (loc.x - x) / loc.y;
        return (- atan(a) / loc.y);
    }

}

Point2d  RoomRayModel::align(Point2d point, Point2d wallvector){
    //normalize wall vector
    wallvector.normalize();
    float x = wallvector.x * point.x + wallvector.y * point.y;
    float y = -1.0f*wallvector.y * point.x + wallvector.x * point.y;
    return Point2d(x,y);
}

//this returns the gain, can be used for both input and output
float  RoomRayModel::xAlignedIntegration(Point2d loc, Point2d ptStart, Point2d ptEnd, bool listLoc){
    Point2d wallVector = ptEnd - ptStart;
    
    Point2d alignedStart = align(ptStart, wallVector);
    Point2d alignedEnd = align(ptEnd, wallVector);
    Point2d alignedLoc = align(loc, wallVector);
    
    alignedEnd = alignedEnd - alignedStart;
    alignedLoc = alignedLoc - alignedStart;
    
    float endVal = integrationSimple(alignedLoc, alignedEnd.x, listLoc);
    float startVal = integrationSimple(alignedLoc, 0.0f, listLoc);
    
   // printf("endval %f startval %f \n", endVal, startVal);
    return fabs(endVal - startVal);
    
}


void RoomRayModel::setRayTracingPoints(Point2d* bouncePoints, Point2d ssLoc, float rheight, float rwidth, int numpoints, float* outputGains, float* inputGains, Point2d listloc){
    
    listenerLoc = listloc;
    soundSourceLoc = ssLoc;
    float yBot = 0.0f-ssLoc.y;
    float yTop = rheight - ssLoc.y;
    float xLeft = 0.0f - ssLoc.x;
    float xRight = rwidth - ssLoc.x;
    
    float w = ssLoc.x;
    float h = ssLoc.y;
    
    for (int i = 0; i < numpoints/2; i++){
        float angle = (360.f / float(numpoints)) * float(i);
       // printf("Angle : %f \n", angle);
        float m = 1.0f/tan(angle * M_PI / 180.f);
        //y = mx + 0
        Point2d pointArray[4] = {Point2d(yBot/m, yBot),
            Point2d(yTop/m, yTop),
            Point2d(xLeft, m*xLeft),
            Point2d(xRight, m*xRight)};
        
        for (int j = 0; j< 4; j++){
            float xO = pointArray[j].x + ssLoc.x;
            if (xO > rwidth or xO < 0.0f){
                pointArray[j].mark = false;
                continue;
            }
            float yO = pointArray[j].y + ssLoc.y;
            if (yO > rheight or yO < 0.0f){
                pointArray[j].mark = false;
                continue;
            }
            if (pointArray[j].mark == true){
                //check for x value
                if (pointArray[j].x >= 0){
                    bouncePoints[i].x = pointArray[j].x + w;
                    bouncePoints[i].y = pointArray[j].y + h;
                }
                else{
                    bouncePoints[i+numpoints/2].x = pointArray[j].x + w;
                    bouncePoints[i+numpoints/2].y = pointArray[j].y + h;
                }
            }
        }
    }
    
    Point2d prevStart = bouncePoints[0];
    
  //  printf("List loc %f %f ssloc %f %f \n", listenerLoc.x, listenerLoc.y, soundSourceLoc.x, soundSourceLoc.y);
    // set the points at even but randomly jittered locations
    for (size_t i = 1; i < numpoints; i++) {

            Point2d start = prevStart;
            Point2d difference = (bouncePoints[i] - bouncePoints[i-1]).scalarMul(0.5f);
            Point2d end = bouncePoints[i-1] + difference;
        
     //   printf("Start : %f %f , end : %f %f \n", start.x, start.y, end.x, end.y);
            outputGains[i-1] =  sqrtf( xAlignedIntegration(listloc, start, end, true));
            inputGains[i-1] = sqrtf(xAlignedIntegration(ssLoc, start, end, false));
//        printf("i : %d OutputGains : %f \n", i,outputGains[i-1]);
//        printf("i : %d InputGains : %f \n", i,inputGains[i-1]);

            prevStart = Point2d(end.x, end.y);
        
        
    }
    
    //do the last gain
    Point2d end = bouncePoints[numpoints-1];
    outputGains[numpoints-1] = sqrtf(xAlignedIntegration(listloc, prevStart, end, true));
    inputGains[numpoints-1] =  sqrtf(xAlignedIntegration(ssLoc, prevStart, end, false));
    
    printf("OutputGains : %f \n", outputGains[numpoints-1]);
    printf("InputGains : %f \n", inputGains[numpoints-1]);
//        printf("S loc x : %f , y : %f \n", ssLoc.x, ssLoc.y);
//        for (int i = 0; i<numpoints;i++){
//            printf("%f,", bouncePoints[i].x);
//        }
//    printf("\n\n");
//    for (int i = 0; i<numpoints;i++){
//        printf("%f,", bouncePoints[i].y);
//    }
}
