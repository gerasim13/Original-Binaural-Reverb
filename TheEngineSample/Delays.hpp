//
//  Delays.hpp
//  TheEngineSample
//
//  Created by Natalie Agus on 13/11/15.
//  Copyright © 2015 A Tasty Pixel. All rights reserved.
//

#ifndef Delays_hpp
#define Delays_hpp

#include <stdio.h>
#include "Point2d.hpp"
typedef struct Delays {
    Delays(){
        delay = 0.0f;
        channel = 0;
        inputGainIndex = 0;
        outputGainIndex = 0;
        extraDelay = false;
    };
    Delays(float delay, size_t channel, float inputGain, float outputGain, bool extraDelay, Point2d bp){
        this->delay = delay;
        this->channel = channel;
        this->inputGainIndex = inputGain;
        this->outputGainIndex = outputGain;
        this->extraDelay = extraDelay;
        this->bp = bp;
        
    }
    float delay;
    size_t channel;
    float inputGainIndex;
    float outputGainIndex;
    bool extraDelay;
    Point2d bp;
    bool operator < (const Delays &d) const{
        return (this->delay < d.delay);
    }
    
} Delays;


#endif /* Delays_hpp */
