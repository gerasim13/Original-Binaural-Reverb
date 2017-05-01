//
//  Gains.hpp
//  surfaceIntegration
//
//  Created by Natalie Agus on 31/8/16.
//  Copyright © 2016 Hans. All rights reserved.
//

#ifndef Gains_hpp
#define Gains_hpp
#include "MonteCarlo.hpp"
#include "Cuboid.hpp"
#include "Vector3D.hpp"
#include "Plane3D.hpp"
#include "Cuboid.hpp"
#include <math.h>
#include <stdio.h>
#import "energyTerm.hpp"


#define DMIN 0.5f

//change them to make the first order reflections lose energy if you have to.
#define NUM_MONTECARLO 20
#define KD 1.0f
#define KS 0.0f
#define ENERGYINITIAL 1.0f
#define ALPHA 0.0f

typedef struct Gains{
    
    Gains(){
    
    }
    //Constructor with arguments
    Gains(float dmin, int numberDelays, float totalSurfaceArea, float* feedbackTapGains){
        this->dmin = dmin ;
        this->numberDelays = numberDelays;
        this->totalSurfaceArea = totalSurfaceArea;
        this->totalInputEnergy = 0.f;
        this->correctInputEnergy = 4.f * M_PI * powf(dmin, 2) * ENERGYINITIAL;
        this->feedbackTapGains = feedbackTapGains;
    };
    
    
    
    //methods
    float calculateGains(Plane3D* surfaces, Vector3D L, Vector3D S);
    

    void monteCarloUpsilon(Vector3D* points, Vector3D L, Vector3D S, Vector3D N, int numPoints, float* up, float area);
    void monteCarloBeta(Vector3D* points, Vector3D L, Vector3D S, Vector3D N, int numPoints, float* beta, float area);
    float calculateBeta(Vector3D x, Vector3D L, Vector3D S, Vector3D N, float visibilitySX, float visibilityXL);
    float integrateUpsilon(Vector3D x, Vector3D L, Vector3D N, float visibility);
    float getDBRDF();
    float pointCollectionFunction(Vector3D x, Vector3D L, Vector3D N, float visibility, float absorption);
    float reflectionKernel(Vector3D x, Vector3D L, Vector3D S, Vector3D N, float visibility);
    Vector3D Lambda(Vector3D u, Vector3D x);
    float phongBRDF(Vector3D S, Vector3D L, Vector3D N);
    void getGains(float* inputGains, float* outputGains);
    

    
    void cartesianToSpherical(Vector3D x, float* angles);
    Vector3D getDirectionVector(Vector3D S, Vector3D N);
    
    float lookUpNormalizationEnergyTerm(Vector3D S);
    
    //variables
    float dmin;
    float totalSurfaceArea;
    int numberDelays;
    
    float* upsilon;
    float* mu;
    float* beta;
    float* feedbackTapGains;

    
    float totalInputEnergy;
    float correctInputEnergy;
    
    energyTerm brdfPhong;
    
    
}Gains;

#endif /* Gains_hpp */
