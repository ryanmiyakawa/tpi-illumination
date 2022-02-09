/**
 * TPI illumination generator v0.0.1
 * Ryan Miyakawa
 * rhmiyakawa@lbl.gov
 * 
 */


#include <stdio.h>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string.h>

using namespace std;



// Squares a number, we do this a lot
double square(double x){
    return x*x;
}

// Writes coordinates to a text file 
void writeCoords(double * x, double * y, int numCoords){
    FILE * outputTextFile = NULL;

    if ((outputTextFile = fopen("illumination_coords.txt", "wb")) == NULL)
                printf("Cannot open file.\n");

    for (int k = 0; k < numCoords; k++){
        fprintf(outputTextFile, "%0.6f, %0.6f\n", x[k], y[k]);
    }
}

// Specify custom boolean masks for "primative" shapes.  
// All pupil fills are designed be composed of primative shapes
bool isInPrimativeShape(double r, double th, int primativeIdx, double * args){
    double sigmaLow, sigmaHigh, wedgeWidth, poleOffset, poleRadius, deltaSigma;
    switch (primativeIdx){
        case 0: // Annulus/disk
            sigmaLow = args[0];
            sigmaHigh = args[1];
            return r > sigmaLow && r < sigmaHigh;
            
        case 1: // Wedge pole
            sigmaLow = args[0];
            sigmaHigh = args[1];
            wedgeWidth = args[3] * M_PI/180;
            return r > sigmaLow && r < sigmaHigh &&
                    th > -wedgeWidth/2 && th < wedgeWidth/2;
            
        case 2: // Circular Pole
            poleOffset = args[0];
            poleRadius = args[1];

            return square(r*cos(th) - poleOffset) + square(r*sin(th)) < square(poleRadius);
            
        case 3: // Leaf Dipole
            poleOffset = 2 - args[0];
            poleRadius = 1;

            return square(r*cos(th) - poleOffset) + square(r*sin(th)) < square(poleRadius);


        case 4: // Leaf Quasar Pole
            deltaSigma = args[0];
            poleRadius = 1;

            return square(r*cos(th) - (2 - deltaSigma)) + square(r*sin(th)) < square(poleRadius) &&
            square(r*cos(th) - (1 - deltaSigma/2)) + square(r*sin(th) - (1 - deltaSigma/2)) < square(poleRadius) &&
            square(r*cos(th) - (1 - deltaSigma/2)) + square(r*sin(th) + (1 - deltaSigma/2)) < square(poleRadius);
            
        case 5: // Leaf Hexapole Pole
            deltaSigma = args[0];
            poleRadius = 1;

            return  square(r*cos(th + M_PI/6) - (2 - deltaSigma)) + square(r*sin(th + M_PI/6)) < square(poleRadius) &&
                    square(r*cos(th - M_PI/6) - (2 - deltaSigma)) + square(r*sin(th - M_PI/6)) < square(poleRadius);

    }
    return true;
}



/**
 * Populates x and y with coordinates of the primative shape given by primativeIdx
 * 
 * xBuff and yBuff are helper arrays for storing coordinates and should be larger than x and y
 * 
 * Returns the number of points used to create this primative
 */
int generatePrimative(double * x, double * y, double * xBuff, double * yBuff, int numberPoints, int targetPoints, 
                    int primativeIdx, double dl, double * radii, int numRadii, bool useRaster, double * params){

    int numTheta;
    double radius;
    double dTh = 0;
    double theta = 0;

    int count = 0;

    for (int k = 0; k < numRadii; k++){
        // Theta stepsize:
        dTh = dl/(2*M_PI * radii[k]);
        numTheta = 2*M_PI / dTh;

        // printf("NumTheta in ring %d: %d\n", k, numTheta);

        for (int m = 0; m < numTheta; m++){

            // Reverse scan direction for odd k:
            if (!useRaster || k % 2 == 0 ){
                theta = -M_PI + m*dTh;
            } else {
                theta =  M_PI - m*dTh;
            }

            if (isInPrimativeShape(radii[k], theta, primativeIdx, params)){
                xBuff[count] = radii[k]*cos(theta);
                yBuff[count] = radii[k]*sin(theta);
                count++;
            }
        }
    }

    for (int k = 0; k < min(count, numberPoints); k++){
        x[k] = xBuff[k];
        y[k] = yBuff[k];
    }

    return count - targetPoints;
}


// Discretized Newton's method for finding optimal theta sampling
int secantSolvePrimative(double * x, double * y, int numberPoints, int targetPoints, 
            int primativeIdx, double radiusWidth, double &drOpt, bool useRaster, double * params){


    // Set up R,Th grid:
    int numRadii = (int) ((1 - radiusWidth)/radiusWidth);
    double radii[numRadii];

    for (int k = 0; k < numRadii; k++){
        radii[k] = radiusWidth*(k + 1);
    }

    double dlGuess = M_PI/ ((double) numberPoints * radiusWidth) * 10;

    double tolX = 0.001;
    int maxIter = 50;
    
    //Make guesses:
    double x1, x2, fxm1, fxm2, R0, S0;

    double xBuff[numberPoints*3];
    double yBuff[numberPoints*3];

    x1 = dlGuess;
    x2 = x1*1.2;
    for (int currentIter = 0; currentIter < maxIter; currentIter++){
        fxm1 = (double) generatePrimative(x, y, xBuff, yBuff, numberPoints, targetPoints, primativeIdx, x1, radii, numRadii, useRaster, params);
        fxm2 = (double) generatePrimative(x, y, xBuff, yBuff, numberPoints, targetPoints, primativeIdx, x2, radii, numRadii, useRaster, params);

        if (fxm1 != fxm2){
            R0 = x1 - fxm1*(x1 - x2) / (fxm1 - fxm2);
            S0 = 1/x1 - fxm1*(1/x1 - 1/x2)/(fxm1 - fxm2);       
        } else {
            // Force convergence:
            S0 = 1/x1;
        }

        if (abs(S0 - 1/x1) < tolX){
            printf("Secant solve converged in %d iterations \n", currentIter);
            printf("optimal dl = %0.5f\n", x1);

            fxm1 = generatePrimative(x, y, xBuff, yBuff, numberPoints, targetPoints, primativeIdx, x1, radii, numRadii, useRaster, params);

            printf("Number of points in primative: %d\n", ((int) fxm1 ) + targetPoints);

            // Set dA:
            drOpt = sqrt(radiusWidth * x1);

            return ((int) fxm1 ) + targetPoints;
        }
    
        //Set new guesses
        x2 = x1;
        x1 = 1/S0;
    }
    printf("MAXIMUM ITERATIONS REACHED\n");
    return x1;
}




// Makes multiple copies of the primative shape, rotating around the unit circle.
// Can offset shapes
// Reconciles total number of points by either adding or removing points
void tileAndReconcile(double * x, double * y, int numberPoints, 
        int numberPointsPrim, int numTiles, double tileAngleSep, double tileAngleOffset){

    double ptR, ptTh;
    int count = 0;

    // Copy the primative into N tiles and applying angle offset.  Offset the parent last 
    // to not double count the offest;
    for (int n = numTiles - 1; n >= 0; n--){
        for (int k = 0; k < numberPointsPrim; k++){

            // Check if we are about to overrun array
            if (k + numberPointsPrim * n > numberPoints){
                // we have overrun the array
                break;
            }

            // Define rotated coordinates
            ptR = sqrt(square(x[k]) + square(y[k]));
            ptTh = atan2(y[k], x[k]) + tileAngleOffset + n*tileAngleSep;

            x[k + numberPointsPrim * n] = ptR*cos(ptTh);
            y[k + numberPointsPrim * n] = ptR*sin(ptTh);
            count++;
        }
    }

    // Reconcile count:
    int pointDifference = 0;
    if (count < numberPoints){

        // Need to add extra points to match total with desired number points
        pointDifference = numberPoints - count;

        for (int m = numberPoints - pointDifference; m < numberPoints; m++){
            x[m] = x[numberPoints - pointDifference - 1];
            y[m] = y[numberPoints - pointDifference - 1];
        }

    } else if (numberPointsPrim > numberPoints){
        // Need to get rid of some points but should never get here 

    }
}

// Generates illumination coordinates by building pupil fills from primative shapes that are subsequently tiled and offset/rotated
void getIlluminationCoordinates(double * x, double * y, int numberPoints, int patternNumber, double * params){



    int numberPointsPrim = 0;
    int targetPoints;
    int numTiles;
    int primativeNumber;
    double tileAngleSep;
    double tileAngleOffset;
    bool useRaster = true;

    switch(patternNumber){
        case 0: // Circular monopole
            targetPoints    = numberPoints;
            primativeNumber = 2;
            numTiles        = 1;
            tileAngleSep    = 0;
            tileAngleOffset = params[2] * M_PI/180;

        break;
        case 1: // Disk/annulus
            targetPoints    = numberPoints;
            primativeNumber = 0;
            numTiles        = 1;
            tileAngleSep    = 0;
            tileAngleOffset = M_PI_2;
            useRaster       = false;

        break;
        case 2: // Circular Dipole
            targetPoints    = numberPoints/2;
            primativeNumber = 2;
            numTiles        = 2;
            tileAngleSep    = M_PI;
            tileAngleOffset = params[2] * M_PI/180;

        break;
        case 3: // Crosspole (wedge dipole)
            targetPoints    = numberPoints/2;
            primativeNumber = 1;
            numTiles        = 2;
            tileAngleSep    = M_PI;
            tileAngleOffset = params[2] * M_PI/180;

        break;
        case 4: // Leaf Dipole 
            targetPoints    = numberPoints/2;
            primativeNumber = 3;
            numTiles        = 2;
            tileAngleSep    = M_PI;
            tileAngleOffset = params[1] * M_PI/180;

        break;
        case 5: // Circular Quadrupole
            targetPoints    = numberPoints/4;
            primativeNumber = 2;
            numTiles        = 4;
            tileAngleSep    = M_PI/2;
            tileAngleOffset = params[2] * M_PI/180;
        break;
        case 6: // Quasar (wedge quadrupole)
            targetPoints    = numberPoints/4;
            primativeNumber = 1;
            numTiles        = 4;
            tileAngleSep    = M_PI/2;
            tileAngleOffset = params[2] * M_PI/180;

        break;
        case 7: // Leaf Quasar 
            targetPoints    = numberPoints/4;
            primativeNumber = 4;
            numTiles        = 4;
            tileAngleSep    = M_PI/2;
            tileAngleOffset = params[1] * M_PI/180;
        break;
        case 8: // Circular N-pole
            targetPoints    = numberPoints/params[3];
            primativeNumber = 2;
            numTiles        = params[3];
            tileAngleSep    = 2*M_PI/params[3];
            tileAngleOffset = params[2] * M_PI/180;
        break;
        case 9: // Wedge N-pole
            targetPoints    = numberPoints/params[4];
            primativeNumber = 1;
            numTiles        = params[4];
            tileAngleSep    = 2*M_PI/params[4];
            tileAngleOffset = params[2] * M_PI/180;
        break;
        case 10: // Leaf Hexapole
            targetPoints    = numberPoints/6;
            primativeNumber = 5;
            numTiles        = 6;
            tileAngleSep    = M_PI/3;
            tileAngleOffset = params[1] * M_PI/180;
        break;
    }

    
    const double radiusWidth = 0.02;
    double drOpt = 0;

    // Numerically solve for optimal dl:
    numberPointsPrim = secantSolvePrimative(x, y, numberPoints, targetPoints, primativeNumber, radiusWidth, drOpt, useRaster, params);


    printf("Optimal dr = %0.5f\n", drOpt);
    // printf("Recomputing grid\n");
    // numberPointsPrim = secantSolvePrimative(x, y, numberPoints, targetPoints, primativeNumber, drOpt, drOpt, useRaster, params);
    
    printf("Tiling primitives\n");
    tileAndReconcile(x, y, numberPoints, numberPointsPrim, numTiles, tileAngleSep, tileAngleOffset);

    
}




int main(int argc, char** argv)
{

    char** argv_test;
    argv_test = argv + 1;

    int patternNumber, numberPoints;
    double param0, param1, param2, param3, param4;

    // Get commandline parameters
    if (argc == 1) {
        printf("No input parameters, setting default params \n");
        
        numberPoints    = 10000;
        patternNumber   = 10;
        param0       = 0.2;
        param1       = 0.3;
        param2       = 0;
        param3       = 0;
        param4       = 0;

        // Define some defaults for each pattern number for testing purposes only
        switch(patternNumber){
            case 0: // Offset monopole:
                param0 = 0.4;
                param1 = 0.15;
                break;
            case 1:
                param0 = 0.3;
                param1 = 0.9;
                break;
            case 2:
                param0 = 0.8;
                param1 = 0.1;
                param2 = 60;
                break;
            case 3:
                param0 = 0.3;
                param1 = 0.8;
                param2 = 0;
                param3 = 22.5;
                break;
            case 4:
                param0 = 0.4;
                param1 = 90;
                break;
            case 5:
                param0 = 0.6;
                param1 = 0.2;
                param2 = 22.5;
                break;
            case 6:
                param0 = 0.3;
                param1 = 0.8;
                param2 = 0;
                param3 = 22.5;
                break;
            case 7:
                param0 = 0.4;
                param1 = 90;
                break;
            case 8:
                param0 = 0.6;
                param1 = 0.15;
                param2 = 22.5;
                param3 = 7;
                break;
            case 9:
                param0 = 0.3;
                param1 = 0.8;
                param2 = 0;
                param3 = 25;
                param4 = 7;
                break;
            case 10:
                param0 = 0.6;
                param1 = 90;
                break;
        }
       

    } else {
        patternNumber   = atoi(*(argv_test++));
        numberPoints    = atoi(*(argv_test++));

        param0       = atof(*(argv_test++));
        param1       = atof(*(argv_test++));
        param2       = atof(*(argv_test++));
        param3       = atof(*(argv_test++));
        param4       = atof(*(argv_test++));
    }


    
    

    
    printf("Pattern number \t\t= %d \n", patternNumber);
    printf("numberPoints \t\t %d \n", numberPoints);
    printf("Arguments: \t\t [%0.2f, %0.2f, %0.2f, %0.2f, %0.2f]\n", param0, param1, param2, param3, param4); 

    double x[numberPoints];
    double y[numberPoints];

    double params[5] = {param0, param1, param2, param3, param4};


    // Populate arrays x, y with coordinates
    getIlluminationCoordinates(x, y, numberPoints, patternNumber, params);

    // Output coordinates to a text file
    writeCoords(x, y, numberPoints);

    return 0;
}




