#include "FDN.h"
#include <math.h>
#include <cstdlib>
#include <cstring>
#include <assert.h>
#include <random>
#include "Cuboid.hpp"
#include "Gains.hpp"
#include <stdlib.h>     /* srand, rand */
#include <random>


void FDN::setParameterSafe(Parameter params)
{
    
    clock_t begin = clock();
    
    
    parametersFDN = newParametersFDN;
    reverbOn = parametersFDN.roomSize > 0.05;
    
    printf("RT60: %f, room width : %f, room length : %f roomHeight: %f\n", parametersFDN.RT60, parametersFDN.roomWidth, parametersFDN.roomHeight, parametersFDN.roomCeiling);
    printf("Listenerloc : %f %f ssloc : %f %f \n", parametersFDN.listenerLoc.x, parametersFDN.listenerLoc.y, parametersFDN.soundSourceLoc.x, parametersFDN.soundSourceLoc.y);
    
    randomSeed = 0; // reset the seed every time for consistent results
    
    // this is required to keep the mixing matrices unitary
    matrixAttenuation = 1.0/sqrt((float)DELAYSPERUNIT);

    Room = Cuboid(parametersFDN.roomWidth, parametersFDN.roomHeight, parametersFDN.roomCeiling);
    
    
//    Room.segmentCube(4 * 4); //even segmentation
    
    //Approximately radial segmentation
//    Room.segmentCubeBasedOnProjectedArea(TOTALDELAYS-SMOOTHDELAY, parametersFDN.soundSourceLoc, parametersFDN.listenerLoc);
    Room.sliceCube(TOTALDELAYS-SMOOTHDELAY);
    //Set delay times
    Room.getDelayValues(delayTimes, parametersFDN.listenerLoc, parametersFDN.soundSourceLoc, SAMPLE_RATE_F);
    
    printf("Room elements %d \n", Room.elements);
    resetTapAttenuation(parametersFDN.RT60);
    GainValues = Gains(DMIN, Room.elements, Room.area, feedbackTapGains);
    
//    GainValues.brdfPhong = BRDFenergy;
    
    float insufficiency = GainValues.calculateGains(Room.segmentedSides, parametersFDN.listenerLoc, parametersFDN.soundSourceLoc);
    
    
    GainValues.getGains(inputGains, outputGains);

    

    //set the delay channel of the multitap too
    
    totalEnergyAfterAttenuation = GainValues.totalInputEnergy;
    
//    //Doing first order reflection Attenuation
//    for (int i = 0; i< TOTALDELAYS-SMOOTHDELAY; i++){
//        firstOrderReflectionAttenuation[i] = 1.0f;
//       // firstOrderReflectionAttenuation[i] = gain(parametersFDN.RT60, delayTimes[i]);
//        //printf("i: %d First order reflection attenuation from jot : %f  from eyring: %f \n", i, firstOrderReflectionAttenuation[i], gainEyring(delayTimes[i]));
//    }
//    
//    //Adjust the beta
//    float totalAttenuatedEnergy = 0;
//    float totalEnergy = 0;
//    for (int i = 0; i<TOTALDELAYS-SMOOTHDELAY; i++){
//        totalAttenuatedEnergy += firstOrderReflectionAttenuation[i] * GainValues.beta[i] * firstOrderReflectionAttenuation[i] * GainValues.beta[i];
//        
//        totalEnergy += GainValues.beta[i] * GainValues.beta[i];
//        
//    }
//    
//    float attenuationConstant = sqrtf(totalAttenuatedEnergy/totalEnergy);
//    
//    printf("Attenuation constant: %f \n",attenuationConstant);
//    
//    //adjust mu
//   totalEnergyAfterAttenuation = 0;
//    for (int i = 0; i<TOTALDELAYS-SMOOTHDELAY; i++){
//      //  printf("Gainvalues mu%f \n", GainValues.mu[i]);
//        inputGains[i] = inputGains[i] * attenuationConstant;
//        totalEnergyAfterAttenuation += inputGains[i] * inputGains[i];
//    }
//    
//    
//    //adjust correct energy
//    GainValues.correctInputEnergy = totalAttenuatedEnergy/totalEnergy * GainValues.correctInputEnergy;
//    
//    printf("New correct energy: %f \n", GainValues.correctInputEnergy);
   // GainValues.totalInputEnergy = totalEnergyAfterAttenuation;
    
    insufficiency = GainValues.correctInputEnergy - totalEnergyAfterAttenuation;
    
    printf("New input energy after attenuation before adjustment, %f \n", totalEnergyAfterAttenuation);
    
    setDelayChannels();
    setDelayNoOutputLength();

    
    
    //Handle energy excess or insufficiency here
    if (insufficiency > 0.0f){
        //set the mu of the last delay (delay no output) to top up to the energy
        float inputEnergyNoOutputDelayLine = insufficiency + ((float)TOTALDELAYS/float(TOTALDELAYS-SMOOTHDELAY) * GainValues.correctInputEnergy) - GainValues.correctInputEnergy;
        printf("Input energy to top up : %f \n", inputEnergyNoOutputDelayLine);
        float muOnNoOutputDelayLine = sqrtf(inputEnergyNoOutputDelayLine/SMOOTHDELAY);
        
        for (int l = TOTALDELAYS-SMOOTHDELAY; l<TOTALDELAYS; l++){
        inputGains[l] = muOnNoOutputDelayLine;
        }

        
    }
    

    
    else if (insufficiency == 0.0f){
        float inputEnergyOnThisDelayLine = ((float)TOTALDELAYS / float(TOTALDELAYS-SMOOTHDELAY) * GainValues.correctInputEnergy)  - GainValues.correctInputEnergy;
        float muOnThisDelayLine = sqrtf(inputEnergyOnThisDelayLine);
        inputGains[TOTALDELAYS-SMOOTHDELAY] = muOnThisDelayLine;
    }
    
    else if (insufficiency < 0.0f){
        excessEnergy = true;
        float currentInputEnergy = totalEnergyAfterAttenuation;
        
        int multiTapIndex = 0;
        
        while (currentInputEnergy > GainValues.correctInputEnergy){
            
            float maxInGain = 0.0f;
            float maxInGainSq = 0.0f;
            int indexMaxInGain = 0;
            
            for (int i = 0; i<TOTALDELAYS-SMOOTHDELAY; i++){
                if (inputGains[i] * inputGains[i] > maxInGainSq){
                    maxInGain = inputGains[i];
                    maxInGainSq = inputGains[i] * inputGains[i];
                    indexMaxInGain = i;
                }
            }
            
            printf("Index : %d caused excess, delay length: %d, removing this\n", indexMaxInGain, delayTimes[indexMaxInGain]);
            
            size_t dVal = delayTimes[indexMaxInGain];
            indicesL[multiTapIndex] = dVal;
            float betaVal = inputGains[indexMaxInGain] * outputGains[indexMaxInGain];
            gainsL[multiTapIndex] = betaVal;
            printf("Sides number %d with area %f \n", indexMaxInGain, Room.segmentedSides[indexMaxInGain].getArea());
            multiDelayLinePoints[multiTapIndex] = Room.segmentedSides[indexMaxInGain].getMidpoint();
            multiTapGains[multiTapIndex] =  gain(parametersFDN.RT60, delayTimes[indexMaxInGain]);
            
            //transform this into no inputGain
            inputGains[indexMaxInGain] = 0.0f;
            
            
            //recalculate total energy
            float newInputEnergy = 0.0f;
            for (int i = 0; i< TOTALDELAYS-SMOOTHDELAY; i++){
                newInputEnergy += inputGains[i] * inputGains[i];
            }
            
            

            currentInputEnergy = newInputEnergy;
            multiTapIndex++;

        }
        
        multiTapDelayTapsNumber = multiTapIndex;
        
        //recalculate total energy
        float newInputEnergy = 0.0f;
        for (int i = 0; i< TOTALDELAYS; i++){
            newInputEnergy += inputGains[i] * inputGains[i];
        }
        
        printf("New input Energy after take off the excess: %f \n" , newInputEnergy);
        //adjust the value of noOutputGain
        float insufficientEnergyDueToRemoval = (GainValues.correctInputEnergy - currentInputEnergy);
        float inputEnergyNoOutputDelayLine = insufficientEnergyDueToRemoval + (((float)(TOTALDELAYS)/(float)(TOTALDELAYS-SMOOTHDELAY)) * GainValues.correctInputEnergy) - GainValues.correctInputEnergy;
        
        
        float muOnNoOutputDelayLine = sqrtf(inputEnergyNoOutputDelayLine/SMOOTHDELAY);
        for (int l = TOTALDELAYS-SMOOTHDELAY; l<TOTALDELAYS; l++){
            inputGains[l] = muOnNoOutputDelayLine;
        }
//        
//        inputGains[TOTALDELAYS-SMOOTHDELAY] = muOnNoOutputDelayLine;
    }
    
    //recalculate total energy
    float newTotalInputEnergy = 0.0f;
    for (int i = 0; i< TOTALDELAYS; i++){
        newTotalInputEnergy += inputGains[i] * inputGains[i];
    }
    
    printf("Previous totalEnergy: %f, New total input energy: %f, correct energy: %f, effectiveEnergy: %f \n", GainValues.totalInputEnergy, newTotalInputEnergy, GainValues.correctInputEnergy, newTotalInputEnergy * (float)(TOTALDELAYS-SMOOTHDELAY)/(TOTALDELAYS));
    
//
//    for (int i = 1; i< TOTALDELAYS; i++){
//        printf("Delay times :%d \n", delayTimes[i]);
//    }
    
//    * @param stereo    is it stereo (true) or mono (false). Mono uses Left only
//    * @param indicesL  array of delay tap times in samples, left channel / mono
//    * @param indicesR  array of delay tap times in samples, right channel
//    * @param gainsL    array of output tap gains for left channel or mono
//        * @param gainsR    array of output tap gains for right channel
//            * @param numTapsL  length of gainsL and indicesL
//            * @param numTapsR  length of gainsR and indicesR
    
    //Initialize the multi-tap delay
    
    if (excessEnergy){
        BMMultiTapDelay_init(&multiTapDelay, false, indicesL, indicesR, gainsL, gainsR, multiTapDelayTapsNumber,0);
        
        setMultiTapDelayChannels();
    }
    
//    for (int i = 1; i< TOTALDELAYS; i++){
//        printf("Delay times :%d \n", delayTimes[i]);
//    }
    
    
//
//    float output = 0.0f;

//    
    //Init multitap delay
//    
//    void BMMultiTapDelay_init(BMMultiTapDelay* This,
//                              bool stereo,
//                              size_t* indicesL,
//                              size_t* indicesR,
//                              float* gainsL,
//                              float* gainsR,
//                              size_t numTapsL,
//                              size_t numTapsR);

    
    setDirectDelayTimes();
    setDirectGains();

    setDirectRayAngles();
    
    
//    for (int i = 1; i< TOTALDELAYS; i++){
//        printf("Delay times :%d \n", delayTimes[i]);
//    }
    
    setTempPoints();
    calculateAdditionalDelays();
    
    
    
    setDirectSingleTapDelay();
    
    
    setFilters();
    

    
    int totalDelayTime = 0;
    for(int i = 0; i < numDelays; i++) totalDelayTime += delayTimes[i];
    


    
    printf("totaldelay  %d \n ", totalDelayTime);
    resetDelay(totalDelayTime);
    resetReadIndices();
    
    setHFDecayMultiplier(1000,1.0f,parametersFDN.RT60);
    



    
    
//    printf("Printing Parameters: \n");
//    
//    for (int i = 0 ; i < numDelays ; i++){
//        printf("Index : %d, delay : %d, inputGain : %f, outputGain: %f \n",i, delayTimes[i], inputGains[i], outputGains[i]);
//        
//
//        
//    }
    
    clock_t end = clock();
    
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("Time elapsed: %f\n", elapsed_secs);
    time[timeindex] = elapsed_secs;
    timeindex++;
    
    if (timeindex > 30){
        timeindex = 0;
    }
    
    
    printf("\n\n======Setting End=======\n\n");
    //        printf("DA %f \n", directAttenuation);
}

// this is the function we call to process the audio from iFretless.
void  FDN::processIFretlessBuffer(float* input, size_t numFrames, float* outputL, float* outputR)
{
    // bypass reverb
    //for (int i = 0; i < numFrames; i++) {
    //    outputL[i] = outputR[i] = ioBuffer[i];
    //}
    // return;
    
    if (parameterNeedsUpdate){
        setParameterSafe(newParametersFDN);
        parameterNeedsUpdate = false;
        for (size_t i = 0; i < numFrames; i++){
            processReverb(&zero_f, outputL + i, outputR + i);
        }
//            printf("DA %f \n", directAttenuation);
    }
    
    else {
        if (reverbOn){
            for (UInt32 i = 0; i < numFrames; i++){
                float* inputPtr = input + i;
                processReverb(inputPtr, outputL + i, outputR + i);
            }
        }
        else{
            // spread the signal to the left and right channels
            memcpy(outputL, input, numFrames * sizeof(Float32));
            memcpy(outputR, input, numFrames * sizeof(Float32));
        }
    }
}


//The structure of the unit is different from the paper. In this case, there's input from previous unit + dry input which will be
//applied to a delay. Output taps are placed on this delay.
//Then, the output from this delay (both dry and previous unit) will be processed back to 4x4 mixing matrix.
//Apply high-shelf first order filter to this, somehow this filter depends on delays


inline void FDN::processReverb(float* pInput, float* pOutputL, float* pOutputR)
{
    
    
    // copy the filtered input so that we can process it without affecting the original value
    float xn = *pInput;
    float fdnTankOutsNew[CHANNELS] = {0};
    float directRaysOutput[2] = { xn *directAttenuation, xn *directAttenuation };
   // printf("%f \n", directAttenuation);
    
    
    // copy output taps to pre-filtered output buffer
    //rwIndices are output taps pointer
    
    vDSP_vgathra (const_cast<const float**>(rwIndices), 1, outputsPF, 1, numDelays);
    
    // scale the output taps according to the feedBackTapGains
    vDSP_vmul(feedbackTapGains, 1, outputsPF, 1, outputsPF, 1, numDelays);
    
    //Processing outputs
    vDSP_vmul(outputGains, 1, outputsPF, 1, outputsBinaural, 1, numDelays);

    if (excessEnergy){
    //get output from multiTaps
    BMMultiTapDelay_processMultiChannelOut(&multiTapDelay, xn, multiTapOutputs);
        
    vDSP_vmul(multiTapGains, 1, multiTapOutputs, 1, multiTapOutputs, 1, multiTapDelayTapsNumber);
    
    //add the output from multiTaps to tankout
    for (int i = 0; i<multiTapDelayTapsNumber; i++){
        size_t channel = multiDelayLinePointsChannel[i];
        fdnTankOutsNew[channel] += multiTapOutputs[i];
    }
    }
    
    processTankOut(fdnTankOutsNew); //convert to 8 channels from outputsBinaural
    
    float fdnTankOutLeft[CHANNELS] = {};
    float fdnTankOutRight[CHANNELS] = {};
    filterChannels(fdnTankOutsNew, directRaysOutput, fdnTankOutLeft, fdnTankOutRight); //HRTF
    
    float reverbOut[2] = {0.0f, 0.0f};
    processDirectRays(directRaysOutput, directRaysOutput); //delays the rays
    addReverbDelay(fdnTankOutLeft, fdnTankOutRight); //Temp point delays
    vDSP_sve(fdnTankOutLeft, 1, &reverbOut[0], CHANNELS);
    vDSP_sve(fdnTankOutRight, 1, &reverbOut[1], CHANNELS);

    *pOutputL = (directRaysOutput[0]*directPortionOn + reverbOut[0]*reverbPortionOn);
    *pOutputR = (directRaysOutput[1]*directPortionOn + reverbOut[1]*reverbPortionOn);
    
    //Continue processing reverb

    // vDSP_vsmul(outputsPF, 1, &parametersFDN.reflection, outputsPF, 1, numTaps);
    
    

    // apply a first order high-shelf filter to the feedback path
    //
    // This filter structure is direct form 2 from figure 14 in section 1.1.6
    // of Digital Filters for Everyone by Rusty Alred, second ed.
    //
    // t = outputsPF + (a1 * z1);
    
    vDSP_vma (a1, 1, z1, 1, outputsPF , 1, t, 1, numDelays);
    // outputsAF = b0*t + b1*z1;
    vDSP_vmma(b0, 1, t, 1, b1, 1, z1, 1, outputsAF , 1, numDelays);
    
    // z1 = t;
    memcpy(z1, t, numDelays * sizeof(float));
    
    // apply the matrix attenuation before mixing
    vDSP_vsmul (outputsAF, 1, &matrixAttenuation, outputsAF, 1, numDelays);
    
    BMFastHadamardTransform(outputsAF, inputs, temp1, temp2, TOTALDELAYS);
    
    // process mixing matrices
    //if (DELAYSPERUNIT == 4) {
    /*
     // The following block of code is replaced by vDSP vectorized code in the section below
     //
     for (i = 0; i < NUMDELAYS; i += DELAYSPERUNIT){
     ppxx = outputsAF[i + 0] + outputsAF[i + 1];
     xxpp = outputsAF[i + 2] + outputsAF[i + 3];
     pnxx = outputsAF[i + 0] - outputsAF[i + 1];
     xxpn = outputsAF[i + 2] - outputsAF[i + 3];
     fbIn[i + 0] = ppxx + xxpp;
     fbIn[i + 1] = ppxx - xxpp;
     fbIn[i + 2] = pnxx + xxpn;
     fbIn[i + 3] = pnxx - xxpn;
     }
     */
    
    // for (i = 0; i < NUMDELAYS; i += DELAYSPERUNIT){
    //    j = (i + 4) % NUMDELAYS;
    //    ppxx = outputsAF[i + 0] + outputsAF[i + 1];
    //    xxpp = outputsAF[i + 2] + outputsAF[i + 3];
    //    pnxx = outputsAF[i + 0] - outputsAF[i + 1];
    //    xxpn = outputsAF[i + 2] - outputsAF[i + 3];
    //
    
//    //replaced with hadamard
//    //  Here's a vectorised version of the code in the comment above:
//    vDSP_vadd(outputsAF + 0 , 4, outputsAF + 1 , 4, ppxxV, 1, delayUnits);
//    vDSP_vadd(outputsAF + 2 , 4, outputsAF + 3 , 4, xxppV, 1, delayUnits);
//    vDSP_vsub(outputsAF + 0 , 4, outputsAF + 1 , 4, pnxxV, 1, delayUnits);
//    vDSP_vsub(outputsAF + 2 , 4, outputsAF + 3 , 4, xxpnV, 1, delayUnits);
//    //
//    //  inputs[j + 0] += ppxx + xxpp;
//    //  inputs[j + 1] += ppxx - xxpp;
//    //  inputs[j + 2] += pnxx + xxpn;
//    //  inputs[j + 3] += pnxx - xxpn;
//    // }
//    // Here's the vectorised version:
//    vDSP_vadd(ppxxV, 1, xxppV, 1, inputs+0, 4, delayUnits);
//    vDSP_vsub(ppxxV, 1, xxppV, 1, inputs+1, 4, delayUnits);
//    vDSP_vadd(pnxxV, 1, xxpnV, 1, inputs+2, 4, delayUnits);
//    vDSP_vsub(pnxxV, 1, xxpnV, 1, inputs+3, 4, delayUnits);
//    //} else { // if DELAYSPERUNIT != 4, i.e. 2x2 mixing
//    

    // mixing in 2x2 matrices
    //for (i = 0; i < NUMDELAYS; i+= 2) {
    //    float a = outputsAF[i];
    //    float b = outputsAF[i+1];
    //    inputs[i] = a + b;
    //    inputs[i+1] = a - b;
    //}
    // vectorized mixing of 2x2 matrices
    //    vDSP_vadd(outputsAF, 2, outputsAF + 1, 2, inputs, 2, DELAYUNITS);
    //    vDSP_vsub(outputsAF, 2, outputsAF + 1, 2, inputs + 1, 2, DELAYUNITS);
    //}
    
    // scale and mix new input together with feedback input
    float scaledInput[TOTALDELAYS] = {};
    
    //Uncomment this  and comment the next two lines if you just want to check FOR
     //   vDSP_vsmul(inputGains, 1, &xn, inputs, 1, numDelays);
    vDSP_vsmul(inputGains, 1, &xn, scaledInput, 1, numDelays);
    
    vDSP_vadd(inputs, 1, scaledInput, 1, inputs, 1, numDelays);
    
    // feedback the mixed audio into the tank, shifting the feedback path over to the next unit
    size_t j;//,k;
    for (j = DELAYSPERUNIT; j < numDelays; j++) *(rwIndices[j]) = inputs[j-DELAYSPERUNIT];
    for (j = 0;j<DELAYSPERUNIT; j++) *(rwIndices[j]) = inputs[j+numDelays-DELAYSPERUNIT];
  //  for (j = 0; j < numDelays; j++) *(rwIndices[j]) = inputs[j];
    
    incrementIndices();
}


/*!
 * Fast Hadamard Transform
 *   works in place
 *
 * @param input   an array of floats with length = length
 * @param output  an array of floats with length = length
 * @param temp1   an array of floats with length = length (temp buffer)
 * @param temp2   an array of floats with length = length (temp buffer)
 * @param length  must be a power of 2, >= 16
 */
inline void FDN::BMFastHadamardTransform(float* input,
                                    float* output,
                                    float* temp1,
                                    float* temp2,
                                    size_t length){
    // length must be a power of 2
    assert(BMPowerOfTwoQ(length));
    // length must be at least 16 because we hard-code the
    // base case at order 16.
    assert(length >= 16);
    
    size_t partitionSize = length;
    float* tmpIn = input;
    float* tmpOut = temp1;
    bool tmpState = true;
    
    // iteratively calculated recursion
    while (partitionSize > 16) {
        size_t halfPartitionSize = partitionSize >> 1;
        for (size_t i = 0; i < length; i += partitionSize){
            // copy all the lower terms into place
            memcpy(tmpOut+i,tmpIn+i,sizeof(float)*halfPartitionSize);
            memcpy(tmpOut+i+halfPartitionSize,tmpIn+i,sizeof(float)*halfPartitionSize);
            // sum all the higher terms into place
            for (size_t j=i; j<halfPartitionSize+i; j++) {
                size_t idx2 = j+halfPartitionSize;
                tmpOut[j] += tmpIn[idx2];
                tmpOut[idx2] -= tmpIn[idx2];
            }
        }
        
        // swap temp buffers to avoid using the same memory for reading and writing.
        tmpIn = tmpOut;
        tmpState = !tmpState;
        if (tmpState) tmpOut = temp1;
        else tmpOut = temp2;
        partitionSize >>= 1;
    }
    
    // base case
    for (size_t i = 0; i < length; i += 16)
        BMFastHadamard16(tmpIn + i, output + i, tmpOut);
}




/*
 * Fast Hadamard transform of length 16
 *   works in place.
 *
 * @param input   an array of 16 floats
 * @param output  an array of 16 floats
 * @param temp16  an array of 16 floats that will be over-written
 *
 */
inline void FDN::BMFastHadamard16(const float* input, float* output, float* temp16){
    // level 1
    // +
    temp16[0] = input[0] + input[8];
    temp16[1] = input[1] + input[9];
    temp16[2] = input[2] + input[10];
    temp16[3] = input[3] + input[11];
    temp16[4] = input[4] + input[12];
    temp16[5] = input[5] + input[13];
    temp16[6] = input[6] + input[14];
    temp16[7] = input[7] + input[15];
    // -
    temp16[8] = input[0] - input[8];
    temp16[9] = input[1] - input[9];
    temp16[10] = input[2] - input[10];
    temp16[11] = input[3] - input[11];
    temp16[12] = input[4] - input[12];
    temp16[13] = input[5] - input[13];
    temp16[14] = input[6] - input[14];
    temp16[15] = input[7] - input[15];
    
    
    // level 2
    //+
    output[0] = temp16[0] + temp16[4];
    output[1] = temp16[1] + temp16[5];
    output[2] = temp16[2] + temp16[6];
    output[3] = temp16[3] + temp16[7];
    //-
    output[4] = temp16[0] - temp16[4];
    output[5] = temp16[1] - temp16[5];
    output[6] = temp16[2] - temp16[6];
    output[7] = temp16[3] - temp16[7];
    
    //+
    output[8] = temp16[8] + temp16[12];
    output[9] = temp16[9] + temp16[13];
    output[10] = temp16[10] + temp16[14];
    output[11] = temp16[11] + temp16[15];
    //-
    output[12] = temp16[8] - temp16[12];
    output[13] = temp16[9] - temp16[13];
    output[14] = temp16[10] - temp16[14];
    output[15] = temp16[11] - temp16[15];
    
    
    // level 3
    // +
    temp16[0] = output[0] + output[2];
    temp16[1] = output[1] + output[3];
    // -
    temp16[2] = output[0] - output[2];
    temp16[3] = output[1] - output[3];
    
    // +
    temp16[4] = output[4] + output[6];
    temp16[5] = output[5] + output[7];
    // -
    temp16[6] = output[4] - output[6];
    temp16[7] = output[5] - output[7];
    
    // +
    temp16[8] = output[8] + output[10];
    temp16[9] = output[9] + output[11];
    // -
    temp16[10] = output[8] - output[10];
    temp16[11] = output[9] - output[11];
    
    // +
    temp16[12] = output[12] + output[14];
    temp16[13] = output[13] + output[15];
    // -
    temp16[14] = output[12] - output[14];
    temp16[15] = output[13] - output[15];
    
    
    // level 4
    output[0] = temp16[0] + temp16[1];
    output[1] = temp16[0] - temp16[1];
    output[2] = temp16[2] + temp16[3];
    output[3] = temp16[2] - temp16[3];
    
    output[4] = temp16[4] + temp16[5];
    output[5] = temp16[4] - temp16[5];
    output[6] = temp16[6] + temp16[7];
    output[7] = temp16[6] - temp16[7];
    
    output[8] = temp16[8] + temp16[9];
    output[9] = temp16[8] - temp16[9];
    output[10] = temp16[10] + temp16[11];
    output[11] = temp16[10] - temp16[11];
    
    output[12] = temp16[12] + temp16[13];
    output[13] = temp16[12] - temp16[13];
    output[14] = temp16[14] + temp16[15];
    output[15] = temp16[14] - temp16[15];
}



// is x a power of 2?
// reference: http://www.exploringbinary.com/ten-ways-to-check-if-an-integer-is-a-power-of-two-in-c/
inline bool FDN::BMPowerOfTwoQ (size_t x)
{
    return ((x != 0) && ((x & (~x + 1)) == x));
}







//ITD for reverb
inline void FDN::addReverbDelay(float* fdnLeft, float*fdnRight){
    for(int i = 0; i < CHANNELS/2; i++){
        //left, delay channels 0-CHANNELS/2
        fdnLeft[i] = reverbDelays[i].process(fdnLeft[i]);
    }
    for (int i = CHANNELS/2 ; i < CHANNELS; i++){
        //right, delay channels CHANNELS/2 - CHANNELS
        fdnRight[i] = reverbDelays[i].process(fdnRight[i]);
    }
}

//ITD for direct rays
inline void FDN::processDirectRays(float* input, float* directRaysOutput){
    directRaysOutput[0] = directRays[0].process(input[0]);
    directRaysOutput[1] = directRays[1].process(input[1]);
}

//Multiplexer
inline void FDN::processTankOut(float fdnTankOut[CHANNELS]){
    for (size_t i = 0; i < numDelays; i++){
        size_t channel = delayTimesChannel[i];
        fdnTankOut[channel] += outputsBinaural[i];
    }
}

//HRTF
inline void FDN::filterChannels(float fdnTankOut[8], float directRay[2], float fdnTankOutLeft[8], float fdnTankOutRight[8]){
    //Filter left & right ear
    for (size_t i = 0; i < CHANNELS; i++){
        fdnTankOutLeft[i] = leftEarFilter[i].process(fdnTankOut[i]);
        fdnTankOutRight[i] = rightEarFilter[i].process(fdnTankOut[i]);
    }
    //Filter direct rays
    directRay[0] = directRayFilter[0].process(directRay[0]);
    directRay[1] = directRayFilter[1].process(directRay[1]);
}

FDN::FDN(void)
{
    bool powerSaveByDefault = false;
    initialise(powerSaveByDefault);
}

FDN::FDN(bool powerSaveMode)
{
    initialise(powerSaveMode);
}

void FDN::initialise(bool powerSaveMode){
    
    numDelays = TOTALDELAYS;
    delayUnits = DELAYUNITSSTD;
    
    delayBuffers = NULL;
    
    parametersFDN =  Parameter();
    
//    BRDFenergy = energyTerm();
    
    setParameterSafe(parametersFDN);

    
}

void FDN::resetReadIndices(){
    // set start, end indices for the first delay in the feedback tank
    rwIndices[0] = startIndices[0] = (float*)delayBuffers;
    endIndices[0] = startIndices[0] + delayTimes[0];
    
    //print delay times
    for (int i = 0; i < numDelays;i++){
        if (delayTimes[i] < 1){
            delayTimes[i] = 10;
        }
    };
    
    // set start / end indices for the second feedback delay tap onwards
    for (int i = 0 + 1; i < numDelays; i++){
        rwIndices[i] = startIndices[i] = startIndices[i-1] + delayTimes[i-1];
        endIndices[i] = startIndices[i] + delayTimes[i];
    }
    
    samplesUntilNextWrap = 0; // Initially, we trigger an unncecssary wrap operation to calculate the number of samples until the next actual wrap
}

// see: http://software.intel.com/en-us/articles/fast-random-number-generator-on-the-intel-pentiumr-4-processor/
// rather than returning an output, this function updates a class variable so that we only have to generate 1 random number for each sample.
inline void FDN::updateRand()
{
    randomSeed = (214013 * randomSeed + 2531011);
    rand = (randomSeed >> 16) & 0x7FFF;
}


void FDN::resetTapAttenuation(float rt60){
    for (int i = 0; i < numDelays; i++){
        updateRand();
        feedbackTapGains[i] = gain(rt60, delayTimes[i]);
      //  printf("Delay time: %d, gain %f \n", delayTimes[i], feedbackTapGains[i]);
    }
}

FDN::~FDN(){
    if (delayBuffers) delete[] delayBuffers;
}

void FDN::resetDelay(int totalDelayTime)
{
    if (delayBuffers) delete[] delayBuffers;
    delayBuffers = new float[totalDelayTime];
    this->totalDelayTime = totalDelayTime;
    
    // clear the buffers
    memset(delayBuffers, 0, sizeof(float) * totalDelayTime);
    memset(inputs, 0, sizeof(float) * numDelays);
    memset(outputsPF, 0, sizeof(float) * numDelays);
    memset(outputsAF, 0, sizeof(float) * numDelays);
    memset(outputsBinaural, 0, sizeof(float) * numDelays);
    
    memset(z1, 0, sizeof(float) * numDelays); //DF2
    memset(z1a, 0, sizeof(float) * numDelays); //DF1
    memset(z1b, 0, sizeof(float) * numDelays); //DF1
    
    randomSeed = 0;
}

inline void FDN::incrementIndices(){
    for (int i = 0; i < numDelays; i++) {
        float** rwIndex = rwIndices + i;
        // increment
        (*rwIndex)++;
        // wrap around back to the beginning if we're off the end
        //if ((*rwIndex) >= endIndices[i]) (*rwIndex) = startIndices[i];
    }
    samplesUntilNextWrap--;
    
    // if any pointer is at the end of the buffer, check all pointers and wrap back to the beginning where necessary
    if (samplesUntilNextWrap <= 0) {
        samplesUntilNextWrap = LONG_MAX;
        for (long i = 0; i< numDelays; i++){
            float** rwIndex = rwIndices + i;
            float** endIndex = endIndices + i;
            
            // wrap all pointers that need to be wrapped
            if ((*rwIndex) >= (*endIndex)) (*rwIndex) = startIndices[i];
            
            // find the distance until the next pointer needs to wrap
            long iThDistanceToEnd = (*endIndex) - (*rwIndex);
            if (iThDistanceToEnd < samplesUntilNextWrap) samplesUntilNextWrap = iThDistanceToEnd;
        }
    }
}

// Multiplies the decay rate for high Frequencies
// hfMult >= 1
// 0 <= fc < ~18000
//
// 1 means HF decays at the same rate as mid frequencies, 2 means HF decay twice as fast,
// 3 is three times as fast, etc.
//
//
// the formulae for calculating the filter coefficients are from Digital Filters for Everyone,
// second ed., by Rusty Alred.  Section 2.3.10: Shelf Filters

void FDN::setHFDecayMultiplier(float fc, float hfMult, float rt60){
    double g[numDelays];
    double D[numDelays];
    double gamma;
//    
//    fc = 8000.0f;
    
    gamma = tan((M_PI * fc) / double(SAMPLINGRATEF));
    
    for (int i = 0; i < numDelays; i++){
        // set the filter gains
        //
        // here we find the gain already applied to each delay line
        double broadbandGain = gain(rt60, delayTimes[i]);
        double desiredHFGain = gain(rt60 / hfMult, delayTimes[i]);
        // and now find how much additional filter gain to apply for the given increase in HF decay speed
        g[i] = desiredHFGain / broadbandGain;
        
        // just a temp variable
        D[i] = 1.0 / ((g[i] * gamma) + 1.0);
        
        // set the filter coefficients
        b0[i] = g[i] * (gamma + 1.0) * D[i];
        b1[i] = g[i] * (gamma - 1.0) * D[i];
        // Rusty Allred omits the negative sign in the next line.  We use it to avoid a subtraction in the optimized filter code.
        a1[i] = -1.0 * ((g[i] * gamma) - 1.0) * D[i];
    }
}


// computes the appropriate feedback gain attenuation
// to get a decay envelope with the specified RT60 time (in seconds)
// from a delay line of the specified length.
//
// This formula comes from solving EQ 11.33 in DESIGNING AUDIO EFFECT PLUG-INS IN C++ by Will Pirkle
// which is attributed to Jot, originally.
double FDN::gain(double rt60, double delayLengthInSamples) {
 //  printf("Gain : %f \n", pow(10.f, (-3.0 *  avgDelay) / (rt60 * SAMPLINGRATEF) ));
   //   printf("Gain 2 : %f \n", pow(M_E, (-3.0 *  avgDelay) / (rt60 * SAMPLINGRATEF) ));
    return pow(10.f, (-3.f *  delayLengthInSamples) / (rt60 * SAMPLINGRATEF) );
    
   // return  pow (0.1 , delayLengthInSamples/SAMPLINGRATEF);
}

void FDN::setParameter(Parameter params){
    newParametersFDN = params;
    parameterNeedsUpdate = true;
}



//Configure channel of each point
void FDN::setDelayChannels(){
    

    for (size_t i = 0; i < TOTALDELAYS-SMOOTHDELAY ; i++){
        delayTimesChannel[i] = determineChannel(Room.segmentedSides[i].getMidpoint().x, Room.segmentedSides[i].getMidpoint().y);
    }
    delayTimes[TOTALDELAYS-SMOOTHDELAY] = 0;
    

    
}

void  FDN::setMultiTapDelayChannels(){
    for (int i = 0; i<multiTapDelayTapsNumber; i++){
        multiDelayLinePointsChannel[i] = determineChannel(multiDelayLinePoints[i].x, multiDelayLinePoints[i].y);
    }
}

//Calculate delay time for each delay tap in the FDN
inline void FDN::setDelayTimes(){
    

    for (int i = 0; i < numDelays ; i++){
        delayTimes[i] = 100+i * (i/2);
    }
    
}

inline void FDN::setDirectGains(){
    
    float SL = parametersFDN.soundSourceLoc.distance(parametersFDN.listenerLoc);
    
//    if(SL<0.5f){
//        SL = 0.5f;
//    }
    
   // float averageDirectDelay = 0.5*((float)directDelayTimes[0] + (float)directDelayTimes[1]);
   // float attenuationEyring = expf(-SOUNDSPEED * averageDirectDelay * 0.0162f / 4.0f / parametersFDN.RT60);
    
    
    directAttenuation = DMIN/(SL) ;

}


//Setting direct delay value for two direct rays
inline void FDN::setDirectDelayTimes()
{
    
    //Calculate delay from source to receiver
    float directDelayLeft = parametersFDN.soundSourceLoc.distance(parametersFDN.listenerLocLeftEar) / SOUNDSPEED;
    float directDelayRight = parametersFDN.soundSourceLoc.distance(parametersFDN.listenerLocRightEar) / SOUNDSPEED;
    directDelayTimes[0] = directDelayLeft;

    directDelayTimes[1] = directDelayRight;
}

//Setting direct ray angle with respect to listener location
void FDN::setDirectRayAngles()
{
    float yDiff, xDiff;
    yDiff = parametersFDN.soundSourceLoc.y - parametersFDN.listenerLocLeftEar.y;
    xDiff = parametersFDN.soundSourceLoc.x - parametersFDN.listenerLocLeftEar.x;
    float angle = atan2(xDiff, yDiff) * 180.0f / M_PI;
    directRayAngles[0] = angle;
    
    float yDiff2, xDiff2;
    yDiff2 = parametersFDN.soundSourceLoc.y - parametersFDN.listenerLocRightEar.y;
    xDiff2 = parametersFDN.soundSourceLoc.x - parametersFDN.listenerLocRightEar.x;
    float angle2 = atan2(xDiff2, yDiff2) * 180.0f / M_PI;
    directRayAngles[1] = angle2;
      printf("directRayAngle L %f R %f \n", directRayAngles[0], directRayAngles[1]);
}




//Set points for additional delay to enter the ears, 8 points, 1 per channel
void FDN::setTempPoints(){
    
    float yBot = 0.0f-parametersFDN.listenerLoc.y;
    float yTop = parametersFDN.roomHeight - parametersFDN.listenerLoc.y;
    float xLeft = 0.0f - parametersFDN.listenerLoc.x;
    float xRight = parametersFDN.roomWidth - parametersFDN.listenerLoc.x;
    
    float w =parametersFDN.listenerLoc.x;
    float h = parametersFDN.listenerLoc.y;
    
    for (int i = 0; i < CHANNELS/2; i++){
        float angle = channeltoangleNormal(i);
        float m = 1.0f/tan(angle * M_PI / 180.f);
        //y = mx + 0
        Vector3D pointArray[4] = {Vector3D(yBot/m, yBot),
            Vector3D(yTop/m, yTop),
            Vector3D(xLeft, m*xLeft),
            Vector3D(xRight, m*xRight)};
        
        for (int j = 0; j< 4; j++){
            float xO = pointArray[j].x + parametersFDN.listenerLoc.x;
            if (xO > parametersFDN.roomWidth or xO < 0.0f){
                pointArray[j].mark = false;
                continue;
            }
            float yO = pointArray[j].y + parametersFDN.listenerLoc.y;
            if (yO > parametersFDN.roomHeight or yO < 0.0f){
                pointArray[j].mark = false;
                continue;
            }
            if (pointArray[j].mark == true){
                //check for x value
                if (pointArray[j].x >= 0){
                    tempPoints[i].x = pointArray[j].x + w;
                    tempPoints[i].y = pointArray[j].y + h;
                }
                else{
                    tempPoints[i+CHANNELS/2].x = pointArray[j].x + w;
                    tempPoints[i+CHANNELS/2].y = pointArray[j].y + h;
                }
            }
        }
    }
    
    
    
}

void FDN::calculateAdditionalDelays(){
    for (int i = 0; i < CHANNELS; i++)
    {
        if (i < CHANNELS/2){ //CHANNEL 0 -3 FOR RIGHT EAR
            //near to right use temp0,leftear- temp0,rightear, 0 to 3, FOR LEFT EAR
            float d = (tempPoints[i].distance(parametersFDN.listenerLocLeftEar) - tempPoints[i].distance(parametersFDN.listenerLocRightEar))  / SOUNDSPEED;
            
            if(d<0.00000001){
                printf("FALSE...\n");
                d = 0.0001;
            }
            additionalDelays[i] = d;
        }
        
        else{ //CHANNEL 4-7 FOR LEFT EAR
            //near to left use temp4,rightear - temp4,leftear, 4 to 7, FOR RIGHT EAR
            float d = (tempPoints[i].distance(parametersFDN.listenerLocRightEar) - tempPoints[i].distance(parametersFDN.listenerLocLeftEar))  / SOUNDSPEED;
            
            if(d<0.00000001){
                printf("FALSE...\n");
                d = 0.0001;
            }
            additionalDelays[i] = d;
        }
//         printf("DELAY ADDITIONAL : %f \n", additionalDelays[i]);
        reverbDelays[i].setTimeSafe(additionalDelays[i]);
    }
    
    }

void FDN::setDirectSingleTapDelay(){
    directRays[0].setTimeSafe(directDelayTimes[0]);
    directRays[1].setTimeSafe(directDelayTimes[1]);
}

void FDN::setFilters(){
    //Right is negative, p.151 DAFX
    for (size_t i = 0; i < CHANNELS; i++){
        float angleL = channelToAngleLeftEar(i);
        float angleR = channelToAngleRightEar(i);
                printf("i %d Channel to Angle : left %f right %f\n", i,angleL, angleR);
        
        leftEarFilter[i].setAngle(channelToAngleLeftEar(i), SAMPLINGRATEF, false);
        rightEarFilter[i].setAngle(channelToAngleRightEar(i), SAMPLINGRATEF, true);
    }
    
    directRayFilter[0].setAngle(directRayAngles[0]+90.f, SAMPLINGRATEF,false);
    directRayFilter[1].setAngle(directRayAngles[1]-90.f, SAMPLINGRATEF, true);
}

    //Helper functions
float FDN::channeltoangleNormal(size_t channel){
        float channelWidth = 360.0f / (float)CHANNELS;
     return (float)channel * channelWidth + (channelWidth / 2.0f);
}


float FDN::channelToAngleLeftEar(size_t channel){
    
    //azimuth between 0 to 360
    float azimuth = channeltoangleNormal(channel);
    
    if (azimuth >= 360.f){
        azimuth -= 360.f;
    }
    
    //change azimuth to be between 0 to 180 and 0 to -180
    if (azimuth > 180.f){
        azimuth = azimuth - 360.f;
    }
    
    if(azimuth > 90.f and azimuth < 180.f){
        float leftEarAngle = azimuth + 90.f - 360.f;
        return leftEarAngle;
    }
    
    return azimuth + 90.f;
    
    
}

float FDN::channelToAngleRightEar(size_t channel){
    //azimuth between 0 to 360
    float azimuth = channeltoangleNormal(channel);
    
    if (azimuth >= 360.f){
        azimuth -= 360.f;
    }
    
    //change azimuth to be between 0 to 180 and 0 to -180
    if (azimuth > 180.f){
        azimuth = azimuth - 360.f;
    }
    
    if(azimuth > -180.f and azimuth < -90.f){
        float leftEarAngle = azimuth - 90.f + 360.f;
        return leftEarAngle;
    }
    
    return azimuth - 90.f;
    
}

//float FDN::channelToAngle(size_t channel){
//    float channelWidth = 360.0f / (float)CHANNELS;
//        size_t midChannel = CHANNELS / 2 - 1;
//    
//    return 45.f;
////    return (float)channel * channelWidth  + (channelWidth / 2.0f);
////        if (channel <= midChannel){
////             printf("Channel : %lu Angle %f \n", channel, (float)channel * channelWidth + (channelWidth / 2.0f));
////            return (float)channel * channelWidth + (channelWidth / 2.0f);
////        }
////        else{
////              printf("Channel : %lu Angle %f \n", channel,(-1.f  * (float) (CHANNELS - channel) * channelWidth )+ (channelWidth / 2.0f));
////            return (-1.f  * (float) (CHANNELS - channel) * channelWidth )+ (channelWidth / 2.0f);
////        }
//////
//}

size_t FDN::angleToChannel(float angleInDegrees){
    
        printf(" original angle %f \n", angleInDegrees);
//    if (angleInDegrees<0){
//        angleInDegrees = fabs(angleInDegrees) + 90.f;
//    }
//    else if(0 <= angleInDegrees and angleInDegrees<= 90){
//        angleInDegrees = 90.f-angleInDegrees;
//    }
//    else if( 90 < angleInDegrees and angleInDegrees <= 180){
//        float difference = 180.f-angleInDegrees;
//        angleInDegrees = 270.f + difference;
//    }
//
//    if(angleInDegrees >= 360.f){
//        angleInDegrees = 360.f - angleInDegrees;
//    }
    if (fabs(angleInDegrees - 0.0f)<0.0001){
        angleInDegrees = 0.0f;
    }
    
    if (angleInDegrees < 0.0f ){
        angleInDegrees = 360.f + angleInDegrees;
    }
    
    float channelWidth = 360.0f / CHANNELS;
    
    //Channel 0 to 7 from 12 o'clock, clockwise according from DFX
    size_t channel =(size_t)floorf(angleInDegrees / channelWidth);
        printf("From angle: %f to channel %lu \n", angleInDegrees, channel);
    return channel;
}

//azimuth is clockwise 0 to 360 degree, p. 149 DAFX
size_t FDN::determineChannel(float x,float y){
    float xL = parametersFDN.listenerLoc.x;
    float yL = parametersFDN.listenerLoc.y;
    
    float xDistance = x - xL;
    float yDistance = y - yL;
    
    float angle2 = atan2f(xDistance, yDistance) * 180.0f / M_PI;
    printf("x %f y %f lLocx %f lLocy %f \n",x,y,xL, yL );
        printf("Angle %f \n", angle2);
    
    return angleToChannel(angle2);
    
}


void FDN::setDelayNoOutputLength(){
    
    //Set delay length for delay-without-output
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<int> uni(50, ROOMCEILING/SOUNDSPEED * SAMPLE_RATE_F); // guaranteed unbiased
    
    
    
    for (int i = 0; i<SMOOTHDELAY; i++){
    int random_integer = uni(rng);
        printf("Randomdelay length: %d \n", random_integer);
    int delayTimeNoOutput = random_integer;
    delayTimes[TOTALDELAYS-SMOOTHDELAY+i] = delayTimeNoOutput;
    }

    
}

float FDN::gainEyring(int delayTime){
    
    float power = -1.0f * SOUNDSPEED * delayTime/SAMPLE_RATE_F * 0.0405 / parametersFDN.RT60;
    
    return expf(power);
    
}
