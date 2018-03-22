//
//  Parameter.hpp
//  TheEngineSample
//
//  Created by Natalie Agus on 28/10/15.
//  Copyright © 2015 A Tasty Pixel. All rights reserved.
//

#ifndef Parameter_hpp
#define Parameter_hpp


#include <stdio.h>
#include "Vector3D.hpp"
#include <math.h>


//
//
#define INITIALROOMSIZE 0.186f


#define INITIALWIDTH 1.0f

#define INITIALLENGTH 1.0f

#define INITIALSOUNDX 0.5f
#define INITIALSOUNDY 2.f/3.f

#define INITIALLISTENERX 0.5f
#define INITIALLISTENERY 1.f/3.f

#define INITIALRT60 2.0f
#define RADIUSOFHEAD 0.08f //8cm radius of head
#define ROOMSIZE 30.f //30 metres max

#define ROOMCEILING 3.0f
#define INITIALDIRECTGAIN 1.f ///(4.f*M_PI);
#define INITIALREVERBGAIN 1.f ///(4.f*M_PI);


typedef struct Parameter {
    
    //Initial values of parameter
    Parameter(){
        this->roomSize = INITIALROOMSIZE * ROOMSIZE;
        this->roomSizeRatio = INITIALROOMSIZE;

        this->roomWidthRatio = INITIALWIDTH;
        this->roomHeightRatio = INITIALLENGTH;
        this->roomWidth = this->roomWidthRatio * this->roomSize;
        this->roomHeight = this->roomHeightRatio * this->roomSize;
        this->roomCeiling = ROOMCEILING;
        
        Vector3D L = Vector3D(INITIALLISTENERX, INITIALLISTENERY);
        setListenerLocation(L);
        Vector3D S = Vector3D(INITIALSOUNDX, INITIALSOUNDY);
        setSoundLocation(S);
        this->listenerXYRatio = Vector3D(INITIALLISTENERX, INITIALLISTENERY);
        this->soundXYRatio = Vector3D(INITIALSOUNDX, INITIALSOUNDY);
    
        this->RT60 = INITIALRT60;
        
        this->directGain = INITIALDIRECTGAIN;
        this->reverbGain = INITIALREVERBGAIN;
        
        this->roomRayModelOn = true;
        this->reflection = 0.95f;

      //  printf("RT60: %f, room width : %f, room length : %f\n", this->RT60, this->roomWidth, this->roomHeight);
       // printf("Listenerloc : %f %f ssloc : %f %f \n", this->listenerLoc.x, this->listenerLoc.y, this->soundSourceLoc.x, this->soundSourceLoc.y);
        
        
    }
    
    void setListenerLocation(Vector3D Ratio);
    void setSoundLocation(Vector3D Ratio);
    void setRoomSize(float size);
    void setWidth(float ratio);
    void setLength(float ratio);
    
    Vector3D listenerXYRatio;
    Vector3D soundXYRatio;
    
    Vector3D listenerLoc;
    Vector3D listenerLocLeftEar;
    Vector3D listenerLocRightEar;
    Vector3D soundSourceLoc;
    float roomSizeRatio, RT60,roomSize;
    float roomWidth, roomHeight;
    float roomWidthRatio, roomHeightRatio;
    float roomCeiling;
    float directGain, reverbGain;
    float reflection;

    bool roomRayModelOn;
    
} Parameter;

#endif /* Parameter_hpp */
