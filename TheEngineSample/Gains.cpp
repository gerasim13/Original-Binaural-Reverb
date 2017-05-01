//
//  Gains.cpp
//  surfaceIntegration
//
//  Created by Natalie Agus on 31/8/16.
//  Copyright © 2016 Hans. All rights reserved.
//

#include "Gains.hpp"
#include <iostream>
#include <algorithm>

#import <Accelerate/Accelerate.h>
using namespace std;

//float Gains::calculateBeta(Vector3D x, Vector3D L, Vector3D S, Vector3D N, float visibilitySX, float visibilityXL){
//    float h = pointCollectionFunction(x, L, N, visibilityXL);
//    float R  = reflectionKernel(x, L, S, N, visibilitySX);
//    return sqrtf(powf(dmin, 2) * h * R);
//}
//
//
//float Gains::integrateUpsilon(Vector3D x, Vector3D L, Vector3D N, float visibility){
//    return pointCollectionFunction(x, L, N, visibility);
//    
//    //    printf("%d %f %f \n", numberDelays, M_PI, totalSurfaceArea);
//}

float Gains::pointCollectionFunction(Vector3D x, Vector3D L, Vector3D N, float visibility, float absorptionRate){
    Vector3D xLnormalized = Lambda(x, L);

    Vector3D xL = L.subtract(x);
//    printf("x %f %f %f, L %f %f %f, xL %f %f %f \t", x.x, x.y, x.z, L.x, L.y, L.z, xL.x, xL.y, xL.z);
    float g = xLnormalized.dotProduct(N) / powf(xL.magnitude(), 2);

    return visibility * g;
    
 //   printf("Point: %f %f %f \n", x.x, x.y, x.z);
//    printf("%f %f %f \n", N.x, N.y, N.z);
//    printf("%f %f %f %f\n", xLnorm.x, xLnorm.y, xLnorm.z, xLnorm.magnitude());
//    printf("%f %f %f \n", xL.x, xL.y, xL.z);
//    printf("%f xl mag \n", xL.magnitude());
//    printf("%f \n", xLnorm.dotProduct(N));
//    printf("%f \n", g);

}

float Gains::getDBRDF(){
    return 1.0f/M_PI;
}


float Gains::reflectionKernel(Vector3D x, Vector3D L, Vector3D S, Vector3D N, float visibility){
    Vector3D xSnormalized = Lambda(x, S);
    Vector3D xS = S.subtract(x);
    float g = xSnormalized.dotProduct(N) / powf(xS.magnitude(), 2);
    return visibility  * g  /M_PI;
 //   return visibility * getDBRDF() * g;
  //  float amt = phongBRDF(S, L, N);
//    printf("phongbrdf: %f \n", amt);
  
    
//    return visibility * phongBRDF(S, L, N) * g;
}

Vector3D Gains::Lambda(Vector3D u, Vector3D x){
    Vector3D ux = x.subtract(u);
//            printf("ux %f %f %f  x %f %f %f u %f %f %f \t", ux.x, ux.y, ux.z, x.x, x.y, x.z, u.x, u.y, u.z);
    ux = ux.normalize();
    return ux;
}



void Gains::monteCarloUpsilon(Vector3D *points, Vector3D L, Vector3D S, Vector3D N, int numPoints, float* up, float area){
    
    //Integrate h over the surface patch
    float hInt = 0.0f;
    for (int i = 0; i<numPoints; i++){
        hInt += (area * pointCollectionFunction(points[i], L, N, 1.0f, 0.0f))/(float)numPoints;
//        printf("PCF %f \t", pointCollectionFunction(points[i], L, N, 1.0f, 0.0f));
    }
    
    *up = sqrtf(((float) numberDelays / (M_PI * totalSurfaceArea)) * hInt  * ENERGYREDUCTIONDIFREB);
    
}

void Gains::monteCarloBeta(Vector3D *points, Vector3D L, Vector3D S, Vector3D N, int numPoints, float* beta, float area){
    
    //Integrate h * R over the surface patch
    float hRInt = 0.0f;
    for (int i = 0; i<numPoints; i++){
        float h = pointCollectionFunction(points[i], L, N, 1.0f, ALPHA);
        float R = reflectionKernel(points[i], L, S, N, 1.0f);
//        printf("Point : %f %f %f \t",points[i].x, points[i].y, points[i].z );
//             printf("beta R %f point %f %f %f \n ", h,points[i].x, points[i].y, points[i].z);
        hRInt += (area * h * R)/(float)numPoints;
    }
    
     *beta = sqrtf(dmin*dmin * hRInt);
    
}

float Gains::phongBRDF(Vector3D S, Vector3D L, Vector3D N){

    float cosAlpha = getDirectionVector(S, N).normalize().dotProduct(N);
    if(cosAlpha < 0.0f){
        cosAlpha = 0.0f;
    }
    
    float n = NSPECULAR;
    
    float normalisationEnergy = lookUpNormalizationEnergyTerm(S);
    
   // return (kD * 1.f/M_PI) + kS * (n+2.f)/(2*M_PI) * powf(cosAlpha, n);
        return (KD * 1.f/(M_PI))* ENERGYREDUCTIONDIF  + normalisationEnergy * powf(cosAlpha, n) * KS;
}


Vector3D Gains::getDirectionVector(Vector3D S, Vector3D N){
    return N.scalarMult(2*S.normalize().dotProduct(N)).subtract(S.normalize());
}



float Gains::calculateGains(Plane3D *surfaces, Vector3D L, Vector3D S){
        printf("SourceLoc : %f %f %f lLoc %f %f %f \n", S.x, S.y, S.z, L.x, L.y, L.z);
    printf("dmin %f \n", dmin);
    mu = new float[numberDelays];
    upsilon = new float[numberDelays];
    beta = new float[numberDelays];
    
    for (int i = 0; i<numberDelays; i++){
        beta[i] = 0.0f;
        upsilon[i] = 0.0f;
        mu[i] = 0.0f;
    }
    
    printf("Number delays : %d \n", numberDelays);
    for (int i = 0; i < numberDelays; i++){
        Vector3D points [NUM_MONTECARLO];
        Vector3D c = surfaces[i].corner;
        Vector3D s1 = surfaces[i].S1;
        Vector3D s2 = surfaces[i].S2;


//        printf("{{%f, %f, %f}, {%f, %f ,%f}, {%f, %f ,%f} },\n", c.x, c.y, c.z, s1.x, s1.y, s1.z, s2.x, s2.y, s2.z);
        
        randomPointsOnRectangle(c, s1, s2, points, NUM_MONTECARLO);
        
        monteCarloUpsilon(points, L, S, surfaces[i].normal, NUM_MONTECARLO, &upsilon[i], surfaces[i].getArea());
//        for (int i = 0; i<numberDelays; i++){
//            upsilon[i] = 1.f;
//        }
        
        monteCarloBeta(points, L, S, surfaces[i].normal, NUM_MONTECARLO, &beta[i], surfaces[i].getArea());
//          printf("Beta i %d : %f \n", i, beta[i]);
//        if (i==2){
//            printf("S: %f %f %f \n", S.x, S.y, S.z);
//            printf("{{%f, %f, %f}, {%f, %f ,%f}, {%f, %f ,%f} }, {%f, %f, %f}\n", c.x, c.y, c.z, s1.x, s1.y, s1.z, s2.x, s2.y, s2.z, surfaces[i].normal.x, surfaces[i].normal.y, surfaces[i].normal.z);
//            printf("beta: %f \n", beta[i]);
//            
//        }
//        printf("Upsilon i %d  is %f \n ",i, (upsilon[i] * upsilon[i]));
    }
    


    vDSP_vdiv(upsilon, 1, beta, 1, mu, 1, numberDelays);
    
    //Dont divide by feedbackTapGains?
    vDSP_vdiv(feedbackTapGains, 1, mu, 1, mu, 1, numberDelays);

    
    
//    std::sort(mu, mu+numberDelays);
//
    float sumbeta = 0.0;
    float sumup = 0.0;
    for (int i = 0; i<numberDelays; i++){
        sumbeta += beta[i];
        sumup += upsilon[i];
    }
    
    printf("Sumbeta: %f sumUp: %f \n", sumbeta, sumup);
//
    for (int i = 0; i< numberDelays; i++){
        totalInputEnergy += mu[i] * mu[i];
       // printf("%f \n", mu[i]);
    }
    
    printf("Total Input Energy: %f \n", totalInputEnergy);
    printf("Input energy should be : %f \n", correctInputEnergy);
    printf("We feed too much input energy by a factor of : %f (if < 1 then we put too little, if > 1 then we put too much)\n", totalInputEnergy/correctInputEnergy);

    
    return correctInputEnergy - totalInputEnergy;

    
}


void Gains::getGains(float *inputGains, float *outputGains){
    //Randomly -1 and +1
    
//     srand (time(NULL));
    

    for (int i = 0; i<numberDelays; i++){
//        int random_int = rand();
//        
//        if(random_int%2 == 0)
//        inputGains[i] = mu[i];
//        
//        else
//            inputGains[i] = mu[i] * -1.f;
//        outputGains[i] = upsilon[i];
        
        inputGains[i] = mu[i];
        outputGains[i] = upsilon[i];
    }
}



void Gains::cartesianToSpherical(Vector3D x, float *angles){
    
    angles[0] = acosf(x.normalize().z);
    angles[1] = acosf(x.normalize().x / sqrtf(powf(x.normalize().x, 2)+ powf(x.normalize().x, 2)));
}

float Gains::lookUpNormalizationEnergyTerm(Vector3D S){
    
    
    //get the theta
    float* angles = new float[2];
    cartesianToSpherical(S, angles);
    

    
    int numberbinsTheta = THETASEGMENTATION;

    

    //check the bin

    int quotientTheta = floorf(angles[0] / (M_PI/2.0f/float(numberbinsTheta))) ;

    
    
    return brdfPhong.SpecularBRDFnormalization[quotientTheta];

}

