//
//  Cuboid.cpp
//  surfaceIntegration
//
//  Created by Natalie Agus on 31/8/16.
//  Copyright © 2016 Hans. All rights reserved.
//

#include "Cuboid.hpp"
#include <algorithm>    // std::sort
using namespace std;

void Cuboid::segmentCube(int tilesPerSide){
    
    //Make sure tilesPerSide can be exactly square-rooted
    try
    {
        if(sqrtf(tilesPerSide) != floorf(sqrtf(tilesPerSide))){
            throw 1;
    }
    }
    catch (int e){
        printf("Please enter some number that can be square rooted for tilesPerSide \n");
        exit(1);
    }
    

    
    segmentedSides = new Plane3D[tilesPerSide*6];
    
    //now divide each side into tiles per side
    dimensions = static_cast<int>(sqrtf(tilesPerSide));

    
    int index = 0;
    for (int k = 0; k<6; k++){
        
        //get the plane
        Plane3D plane = sides[k];
        
        Vector3D s1 = plane.S1.scalarMult((float)1/dimensions);
        Vector3D s2 = plane.S2.scalarMult((float)1/dimensions);
        
        for(int i = 0; i<dimensions; i++){
            for (int j = 0; j<dimensions; j++){
                //Get the corner
                Vector3D c = plane.corner.add(plane.S1.scalarMult((float)j/dimensions));
                c = c.add(plane.S2.scalarMult((float)i/dimensions));
                
                segmentedSides[index] = Plane3D(c,s1,s2);
                segmentedSides[index].normal = plane.normal;
                index = index + 1;
            }
        }
    }
    
    elements = 6*tilesPerSide;
    
    //Printing the result
//    for (int i = 0; i< 6*tilesPerSide ; i++){
//        printf("c : {%f %f %f}, s1: {%f %f %f}, s2: {%f %f %f}, n: {%f %f %f} \n", segmentedSides[i].corner.x, segmentedSides[i].corner.y, segmentedSides[i].corner.z, segmentedSides[i].S1.x, segmentedSides[i].S1.y, segmentedSides[i].S1.z, segmentedSides[i].S2.x, segmentedSides[i].S2.y, segmentedSides[i].S2.z, segmentedSides[i].normal.x, segmentedSides[i].normal.y, segmentedSides[i].normal.z);
//    }
    
}

float Cuboid::segmentCubeOnce(int d){
    //segment floor and ceiling
    dimensions = d;
    
    int index = 0;
    for (int k = 0; k <2; k++){
        Plane3D plane = sides[k];
        
        Vector3D s1 = plane.S1.scalarMult((float)1/dimensions);
        Vector3D s2 = plane.S2.scalarMult((float)1/dimensions);
        
        for(int i = 0; i<dimensions; i++){
            for (int j = 0; j<dimensions; j++){
                //Get the corner
                Vector3D c = plane.corner.add(plane.S1.scalarMult((float)j/dimensions));
                c = c.add(plane.S2.scalarMult((float)i/dimensions));
                
                segmentedSides[index] = Plane3D(c,s1,s2);
                segmentedSides[index].normal = plane.normal;
                index = index + 1;
            }
        }
    }
    
    for (int k = 2; k<6; k++){
        
        //get the plane
        Plane3D plane = sides[k];
        
        segmentedSides[index] = plane;
        index ++;
    }
    
    return index;
    
    
}

float Cuboid::segmentCubeOnce(){
    //segment floor and ceiling
    int index = 0;

    for (int k = 0; k<6; k++){
        
        //get the plane
        Plane3D plane = sides[k];
        
        segmentedSides[index] = plane;
        index ++;
    }

    return index;
    
    
}


void Cuboid::segmentCubeBasedOnProjectedArea(int numDelays, Vector3D S, Vector3D L){
    segmentedSides = new Plane3D[numDelays];
    
    
    elements = numDelays;
    
    
  //  int index = segmentCubeOnce(3);
    int index = segmentCubeOnce();
    
   // int delaysCeilingFloor = (2*3*3);
    int delaysCeilingFloor = 0;
    
    while(index<numDelays){
        
        //get the side with the largest projected area
        int maxIndex = 0;
        float maxArea = 0;
        
        for (int i = delaysCeilingFloor; i<index; i++){

            //add is better
            float comparisonValue =   fabs(ProjectedArea_rectangleSubDiv(segmentedSides[i], S, 5)) + fabs(ProjectedArea_rectangleSubDiv(segmentedSides[i], L, 5)) ;

           // printf("%d subdivided ? \n ", segmentedSides[i].subdivided);
            if (comparisonValue > maxArea && segmentedSides[i].subdivided < numDelays){
                maxArea = comparisonValue;
                maxIndex = i;
            }
            
//            else if (segmentedSides[i].subdivided >= 5){
//                printf("SUBDIVIDED TOO MUCH... \n");
//            }
        }
        
        //split that patch
        index = dividePlane(segmentedSides[maxIndex], index, maxIndex, S, L);
        
     //   index = dividePlaneAlongS1(segmentedSides[maxIndex], index, maxIndex);
        
    }
    
//    printf("Index is finally: %d \n", index);
//    for (int i = 0; i<(index); i++){
//        printf("{{%f, %f, %f}, {%f, %f, %f}, {%f, %f, %f}},", segmentedSides[i].corner.x, segmentedSides[i].corner.y, segmentedSides[i].corner.z, segmentedSides[i].S1.x , segmentedSides[i].S1.y, segmentedSides[i].S1.z, segmentedSides[i].S2.x, segmentedSides[i].S2.y, segmentedSides[i].S2.z );
//    }
////
}

//
////NEED CHANGES
//float Cuboid::projectedAreaOfAPlane(Vector3D S, Vector3D L, Plane3D patch){
//    // (projectedArea(An,L) + projectedArea(An,S)) / delayLength^2
// //   float areaAnL = patch.getArea() * fabs(L.normalize().dotProduct(patch.normal)) / powf(L.distance(patch.getMidpoint()), 2);
//    //float areaAnS = patch.getArea() * fabs(S.normalize().dotProduct(patch.normal)) / powf(S.distance(patch.getMidpoint()), 2);
//    float delayLength = S.distance(patch.getMidpoint()) + L.distance(patch.getMidpoint());
//    return 1 / powf(delayLength, 2);
//  //  return patch.getArea() * fabs(S.normalize().dotProduct(patch.normal)) * (1.f/powf(S.distance(patch.getMidpoint()),2));
//}


int Cuboid::dividePlaneAlongS1(Plane3D divide, int index, int sourceIndex){
    
    Vector3D newCorner = divide.corner.add(divide.S2.scalarMult(0.5));
    Vector3D newS2 = divide.S2.scalarMult(0.5);
    segmentedSides[sourceIndex] = Plane3D(divide.corner, divide.S1, newS2);
    
    segmentedSides[index] = Plane3D(newCorner, divide.S1, newS2);
    index++;
    
    return index;
}

int Cuboid::dividePlane(Plane3D divide, int index, int sourceIndex, Vector3D S, Vector3D L ){
//    
    Vector3D omegaS1L= L.subtract(divide.corner.add(divide.S1.scalarMult(0.5))).normalize();
    Vector3D omegaS1S= S.subtract(divide.corner.add(divide.S1.scalarMult(0.5))).normalize();
    float projectedLengthS1 = divide.S1.magnitude() * omegaS1L.dotProduct(divide.normal) + divide.S1.magnitude() * omegaS1S.dotProduct(divide.normal);
    
    Vector3D omegaS2L = L.subtract(divide.corner.add(divide.S2.scalarMult(0.5))).normalize();
    Vector3D omegaS2S = S.subtract(divide.corner.add(divide.S2.scalarMult(0.5))).normalize();
    float projectedLengthS2 = divide.S2.magnitude() * omegaS2L.dotProduct(divide.normal) + divide.S2.magnitude() * omegaS2S.dotProduct(divide.normal);
//    
//    //project S1 into omega
//    float projectedAreaS1 = omega.normalize().scalarMult( divide.S1.dotProduct(omega)).magnitude();
//    //project S2 into omega
//    float projectedAreaS2 = omega.normalize().scalarMult( divide.S2.dotProduct(omega)).magnitude();
//    
//    
//    
//    
    Vector3D c2 = divide.corner.add(divide.S1.scalarMult(0.5));
    Vector3D s1 = divide.S1.scalarMult(0.5);
    Vector3D c22 = divide.corner.add(divide.S2.scalarMult(0.5));
    Vector3D s2 = divide.S2.scalarMult(0.5);

//    float projectedAreaS1 = fabs(ProjectedArea_rectangleSubDiv(Plane3D(divide.corner, s1, divide.S2), L, 2)) * fabs(ProjectedArea_rectangleSubDiv(Plane3D(divide.corner, s1, divide.S2), S, 2))  + fabs(ProjectedArea_rectangleSubDiv(Plane3D(c2, s1, divide.S2), L, 2)) * fabs(ProjectedArea_rectangleSubDiv(Plane3D(c2, s1, divide.S2), S, 2));
//
//
//    float projectedAreaS2 = fabs(ProjectedArea_rectangleSubDiv(Plane3D(divide.corner, divide.S1, s2), L, 2)) * fabs(ProjectedArea_rectangleSubDiv(Plane3D(divide.corner, divide.S1, s2), S, 2)) + fabs(ProjectedArea_rectangleSubDiv(Plane3D(c22, divide.S1, s2), L, 2)) * fabs(ProjectedArea_rectangleSubDiv(Plane3D(c22, divide.S1, s2), S, 2));
//    
//  

    if (projectedLengthS1 > projectedLengthS2){
        //divide along S1
        int subdivision = segmentedSides[sourceIndex].subdivided;
        segmentedSides[sourceIndex] = Plane3D(divide.corner, s1, divide.S2);
        segmentedSides[sourceIndex].subdivided = subdivision + 1;

        segmentedSides[index] = Plane3D(c2, s1, divide.S2);
        segmentedSides[index].subdivided = segmentedSides[sourceIndex].subdivided;
        
        index++;
    }
    else{
        //Divide along S2
        int subdivision = segmentedSides[sourceIndex].subdivided;
        segmentedSides[sourceIndex] = Plane3D(divide.corner, divide.S1, s2);
        segmentedSides[sourceIndex].subdivided = subdivision + 1;

        segmentedSides[index] = Plane3D(c22, divide.S1, s2);
        segmentedSides[index].subdivided = segmentedSides[sourceIndex].subdivided;
        index++;
        

    }
    return index;
}

//Either S1 or S2
bool Cuboid::longestDimension(Plane3D patch){
    if (patch.S1.magnitude() > patch.S2.magnitude()){
        return true;
    }
    else{
        return false;
    }
}


void Cuboid::getDelayValues(int *delayValues, Vector3D L, Vector3D S, int Hz){
    for (int i =0; i< elements; i++){
        Vector3D p = segmentedSides[i].getMidpoint();
        float d1 = S.subtract(p).magnitude();
        float d2 = L.subtract(p).magnitude();
        delayValues[i] = static_cast<int>((d1+d2)/SOUNDSPEED*Hz);
    }
}




float Cuboid::projAreaSubSec(Plane3D r, Vector3D L){
    // centre of r
    Vector3D c = r.getMidpoint();
    
    // area of r
    float A = r.getArea();
    
    // surface normal of r
    Vector3D n =   r.S1.crossProduct(r.S2).normalize();//  normalize(crossProduct(r.s1, r.s2));
    
    // L-c
    Vector3D cToL = L.subtract(c);
    
    // distance from L to centre of r
    float d = cToL.magnitude();
    
    // direction from c to L
    Vector3D Omega = cToL.normalize();
    
    // projected area = Area (n . Omega) / d^2
    return (A * fabs(n.dotProduct(Omega))) / (d*d);
}




/*!
 * Estimate the area of the rectangle r projected onto the unit sphere
 * centred at L.
 *
 * @param r         the rectangle
 * @param L         point at centre of the unit sphere of projection
 * @param divisions number of times to subdivide the rectangle in each
 *                  direction (l x w). higher values are more accurate.
 */
float Cuboid::ProjectedArea_rectangleSubDiv(Plane3D r, Vector3D L, size_t divisions){
    // how many sub patches?
    size_t numPatches = divisions*divisions;
    
    // allocate memory
    Plane3D* subMesh = new Plane3D[numPatches];
    
    // subdivide r
    Vector3D s1s = r.S1.scalarMult(1.0f/(float)divisions);//  scalarDivide(r.s1,divisions);
    Vector3D s2s = r.S2.scalarMult(1.0f/(float)divisions);;
    for(size_t i=0; i<divisions; i++){
        for(size_t j=0; j<divisions; j++){
            Vector3D cij =  r.corner.add(s1s.scalarMult(i)).add(s2s.scalarMult(j));// add3(r.c,
                             //  scalarProduct(i,s1s),
                             //  scalarProduct(j,s2s));
            Plane3D rij = Plane3D(cij, s1s, s2s);
            subMesh[divisions*i + j] = rij;
        }
    }
    
    // get projected area of each subdivision
    float projArea = 0.0f;
    for(size_t i=0; i<numPatches; i++)
        projArea += projAreaSubSec(subMesh[i],L);
    
    //
    free(subMesh);
    
    return projArea;
}


