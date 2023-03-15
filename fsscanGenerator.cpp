/**
 * TPI illumination generator v0.0.1
 * Ryan Miyakawa
 * rhmiyakawa@lbl.gov
 * 
 * Generates coordinates describing illumination patterns for TPI. Documentation 
 * on the pattern definitions and parameters can be found here:
 * https://cxro.atlassian.net/l/c/x0PcjVrA
 * 
 * Builds illuminations using primative shapes, and then uses a discretized
 * Newton's method (secant solver) to make sure that returned coordiante array
 * lengths match the desired number (default 10,000 points).  Coordinate array
 * lengths are controlled by adjusting samping along the tangential direction keeping
 * radial sampling fixed.  Coordinate sampling is chosen to have fixed point density
 * across the active illumination area in angle space.
 * 
 * Illumination coordinates are written to a text file.
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

/*
 * Populates the arrays x and y with coordinates as specified by pattern number and
 * paramters.  See https://cxro.atlassian.net/l/c/x0PcjVrA
 * for documentation on illuminations and paramters
 */
void getFSCoordinates(double * x, double * y, int numberPoints, int patternNumber, double * params){

    double dummy[5] = {0, 0, 0, 0, 0};
    double a = params[0];
    double b = params[1];   
    double rotation = 0; 
    double lineSpacing, spaceCorrection;

    int numLines, numberSubpoints, lineStart, lineEnd;
    int offsetByNumPoints = 0, offset = 0;
    double P, th, alpha, beta, N, phi, t;

    double * xr;
    double * yr;

    switch(patternNumber){
        case 0: // Ellipse
            rotation = params[2];
            // Make equal arc lengths, although i'm not sure what circumference of ellipse is:
            P = 2*M_PI * (a + b)/2;

            for (int k = 0; k < numberPoints; k++){
                th = k * 2*M_PI/numberPoints;
                x[k] = a*cos(th);
                y[k] = b*sin(th);
            }

        break;
        case 1: // MB spiral
            alpha = params[2];
            N = params[3];
            rotation = params[4];

            for (int k = 0; k < numberPoints; k++){
                phi = k * 2*M_PI * (N + 1)/numberPoints * N;
                
                x[k] = a/(1 + 1/alpha) * 
                        (cos(phi) + 1/alpha * cos(phi/(1 + N)));
                y[k] = b/(1 + 1/alpha) * 
                        (sin(phi) + 1/alpha * sin(phi/(1 + N)));
            }

        break;
        case 2: // Raster

            /**
             * @brief 03/14/23 - Updated raster to interlace lines on the half-integer line so that start and
             * end points are in proximity.  Also added the ability to cascade rasters to make the 2D raster
             * out of 2 X 1D rasters.  This is done by passing in a 5th parameter to the function. 
             */

            rotation = params[3];
            numLines = (int) params[2];

            // If params[4] == 1, then we are rastering the second half of the raster, so offset the array indices
            offsetByNumPoints = (int) params[4];
            if (offsetByNumPoints == 1){
                offset = numberPoints;
            }

            // Require rasters to have at least 2 lines to avoid divide by 0
            if (numLines <= 1){
                numLines = 2;
            }

            numberSubpoints = numberPoints/2;

            spaceCorrection =  ((double)numLines - 1)/((double)numLines - 0.5);
            lineSpacing = 2*b/((double)numLines - 1) * spaceCorrection;
           
            for (int k = 0; k < numLines; k++){
                lineStart = (int) ((k * numberSubpoints)/numLines); 
                lineEnd = (int) (((k + 1) * numberSubpoints)/numLines);

                for (int m = lineStart; m < lineEnd; m++){
                    x[m + offset] = ( 2*(k % 2) - 1) * (-a + (2*a/(lineEnd - lineStart)) * (m - lineStart));
                    if (!offsetByNumPoints){
                        // Reverse order of y coordinates
                        y[m + offset] = -b + lineSpacing * (numLines - k - 1);
                    } else {
                        y[m + offset] = -b + lineSpacing * k;
                    }
                }
            }

            // Now raster backward filling in the other half of the raster
            for (int k = numLines - 1; k >= 0; k--){
                lineStart = (int) ((k * numberSubpoints)/numLines);
                lineEnd = (int) (((k + 1) * numberSubpoints)/numLines);

                for (int m = lineStart; m < lineEnd; m++){
                    x[numberPoints - m - 1 + offset] = ( 2*(k % 2) - 1) * (-a + (2*a/(lineEnd - lineStart)) * (m - lineStart));

                     if (!offsetByNumPoints){
                        y[numberPoints - m - 1 + offset] = -b + lineSpacing * (numLines - k - 1)+ lineSpacing/2; // add a half line
                     }else {
                        y[numberPoints - m - 1 + offset] = -b + lineSpacing * k + lineSpacing/2; // add a half line
                     }
                }
            }


        break;

        case 3: // Interlaced raster

            // 3/14/23 - Now making 2D interlaced raster from 2 X 1D rasters
            params[4] = 0;
            getFSCoordinates(x, y, numberPoints/2, 2, params);

            // Swap params 0 and 1 to make the second raster
            params[4] = params[0];
            params[0] = params[1];
            params[1] = params[4];

            // Flag that we are now rastering the second half of the raster
            params[4] = 1;
            params[3] += 90;

            getFSCoordinates(x, y, numberPoints/2, 2, params);
            return;

            // rotation = params[3];
            // numberSubpoints = numberPoints/2;
            // numLines = (int) params[2];

            // // Require rasters to have at least 2 lines to avoid divide by 0
            // if (numLines <= 1){
            //     numLines = 2;
            // }

            // for (int k = 0; k < numLines; k++){
            //     int lineStart = (int) ((k * numberSubpoints)/numLines);
            //     int lineEnd = (int) (((k + 1) * numberSubpoints)/numLines);

            //     for (int m = lineStart; m < lineEnd; m++){
            //         x[m] = ( 2*(k % 2) - 1) * (-a + (2*a/(lineEnd - lineStart)) * (m - lineStart));
            //         y[m] = -b + (2*b/(numLines - 1)) * k;
            //     }
            // }

            //  for (int k = numLines; k < 2*numLines; k++){
            //     int lineStart = (int) ((k * numberSubpoints)/numLines);
            //     int lineEnd = (int) (((k + 1) * numberSubpoints)/numLines);

            //     for (int m = lineStart; m < lineEnd; m++){
            //         y[m] = ( 2*(k % 2) - 1) * (-b + (2*b/(lineEnd - lineStart)) * (m - lineStart));
            //         x[m] = -3*a + (2*a/(numLines - 1)) * (k - 1);
            //     }
            // }

        break;
        case 4: // Lissajous?
            alpha = params[2];
            beta = params[3];
            rotation = params[4];


            for (int k = 0; k < numberPoints; k++){
                t =  k * 2*M_PI/numberPoints;
                x[k] = a * sin(alpha * t + M_PI_2);
                y[k] = b * sin(beta * t);
            }

        break;
    }

    // Now handle rotation:
    double RX, RY;
    if (rotation != 0){
        th = rotation * M_PI/180; // rotation in radians


        if (offsetByNumPoints){
            for (int k = numberPoints; k < 2*numberPoints; k++){
                RX = cos(th) * x[k] - sin(th) * y[k];
                RY = sin(th) * x[k] + cos(th) * y[k];
                x[k] = RX;
                y[k] = RY;
            }
        } else {
            for (int k = 0; k < numberPoints; k++){
                RX = cos(th) * x[k] - sin(th) * y[k];
                RY = sin(th) * x[k] + cos(th) * y[k];
                x[k] = RX;
                y[k] = RY;
            }\
        }

    }
    
}


int main(int argc, char** argv)
{
    char** argv_test;
    argv_test = argv + 1;

    int patternNumber, numberPoints;
    double param0, param1, param2, param3, param4;


    // If not all inputs are given, then use default paramters
    if (argc < 8) {
        printf("No input parameters, setting default params \n");
        
        // Default parameters
        numberPoints    = 2500;
        patternNumber   = 3;

        // Define some defaults for each pattern number for testing purposes only
        switch(patternNumber){
            case 0: // Ellipse
                param0 = 0.8;
                param1 = 0.2;
                param2 = 22.5;
                break;
            case 1: // Rosetta
                param0 = 0.75;
                param1 = 1;
                param2 = 1.7;
                param3 = 4;
                param4 = 0;
                break;
            case 2: // Raster
                param0 = 0.5;
                param1 = 0.5;
                param2 = 4;
                param3 = 0;
                param4 = 1;
                break;
            case 3: // Interlaced raster
                param0 = 0.6;
                param1 = 0.4;
                param2 = 10;
                param2 = 5;
                param3 = 30;
                break;
            case 4: // Lissajous
                param0 = 0.75;
                param1 = 1;
                param2 = 6;
                param3 = 8;
                param4 = 22.5;
                break;
            case 5: // 2draster
                param0 = 0.5;
                param1 = 0.5;
                param2 = 10;
                param3 = 0;
                break;
        }

       

    } else {

    // Get commandline parameters
        patternNumber   = atoi(*(argv_test++));
        numberPoints    = atoi(*(argv_test++));

        param0          = atof(*(argv_test++));
        param1          = atof(*(argv_test++));
        param2          = atof(*(argv_test++));
        param3          = atof(*(argv_test++));
        param4          = atof(*(argv_test++));
    }

    printf("Pattern number \t\t= %d \n", patternNumber);
    printf("numberPoints \t\t %d \n", numberPoints);
    printf("Arguments: \t\t [%0.2f, %0.2f, %0.2f, %0.2f, %0.2f]\n", param0, param1, param2, param3, param4); 

    // Arrays to hold coordinates
    double x[numberPoints];
    double y[numberPoints];

    double params[5] = {param0, param1, param2, param3, param4};

    // Populate arrays x, y with coordinates
    getFSCoordinates(x, y, numberPoints, patternNumber, params);

    // Output coordinates to a text file
    writeCoords(x, y, numberPoints);

    return 0;
}




