#ifndef __IMAGEPROCESSING_LIB__
#define __IMAGEPROCESSING_LIB__

#include <math.h>
#include <iostream>
#include <iomanip>      // std::setprecision
#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cstdlib>
#include "msaLibs.h"




// ***********************************
//  
// 
template <typename T>
void blurMaskEdgeOnImage(T* I, T* maskI, T*imageOut, double sigma, int kernelLength , const unsigned long int nx, const unsigned long int ny, double thresholdMask = 0.6){

//std::cerr<<"sigma="<< sigma << "   kernelLengthint="<< kernelLength<<"\n";

const int nxy = nx * ny;
const int length=kernelLength;
const int nxK =  (2*length+1);
const int nyK = (2*length+1);
const int kernelSize = nxK*nyK;
float * kernel = new float [ kernelSize ];
float * tmpI = new float [ nxy ];
double angpix = 1;


double sum = 0;
double counter=0;

for (int jjT=-length; jjT<length+1; jjT++){
  for(int iiT=-length; iiT<length+1; iiT++){
          counter++;
        double distanceSquared =  pow(jjT,2.0)+pow(iiT,2.0);
        int X=iiT+length;
        int Y=jjT+length;
        kernel [  X + (Y * nxK) ] = 1 / (sigma * sigma * 2*3.1415927) 
              * exp (- (distanceSquared) / (2.0 * sigma * sigma));
        sum+=kernel [  X + (Y * nxK) ];
  }
}

//correct compensate for discretization (sum needs to be 1):
double normSum = (double((1.0-sum)))/double(kernelSize);
for (int ii=0; ii<kernelSize; ii++){
        kernel [ii] += normSum;
}


sum = 0;
for (int jjT=-length; jjT<length+1; jjT++){
  for(int iiT=-length; iiT<length+1; iiT++){
        //double distanceSquared =  pow(jjT,2.0)+pow(iiT,2.0);
        int X=iiT+length;
        int Y=jjT+length;
        sum+=kernel [ X +Y * nxK ];
 //       std::cerr << "kernel ["<< X +Y * nxK << "]=";
 //       std::cerr << kernel [ X +Y * nxK ] << " ";
  }
 // std::cerr<< " \n";
}
//std::cerr << " \n";
//std::cerr << "sum=" << sum <<" \n";



for (unsigned long int YY=length; YY<ny-length; YY++){
 for(unsigned long int XX=length; XX<nx-length; XX++){
   int numBlack=0;
   int numWhite=0;

   for (int jjT=-length; jjT<length; jjT++){
     for(int iiT=-length; iiT<length; iiT++){
        unsigned long int idx = (XX+iiT)+(YY+jjT)*nx;
        if ( maskI[idx] > 0.5f){
                numWhite++;  
        }else{
                numBlack++; 
        }
     }
   }


   //
   if (numBlack>0 && numWhite>0){
        double sum = 0;
        for (int jjT=-length; jjT<length; jjT++){
          for(int iiT=-length; iiT<length; iiT++){
                unsigned long int idx = (XX+iiT)+(YY+jjT)*nx;
                unsigned long int kernelIdx=iiT+length+(jjT+length)*nxK;
                sum+= kernel[kernelIdx]*I[idx];
          }
        }
        tmpI [ XX + YY * nx ] = sum;
   }else{
     tmpI [ XX + YY * nx ] = I [ XX + YY * nx ];
   }
 }
 

}

 for (unsigned long int ii=0;ii<nxy;ii++){
   imageOut[ii]= tmpI[ii];
 }

delete [] kernel;
delete [] tmpI;
}




template <typename T>
void fillGaussianKernel(T* kernel, int kernelLength, double sigmaSquared){
        const int length=kernelLength;
        const int nxK =  (2*length+1);
        const int nyK = (2*length+1);
        const int kernelSize = nxK*nyK;
        double sum = 0;
        double counter=0;
        double sigmaSquared_double = 2.0*sigmaSquared;
        double sigmaSquared_doublePI = sigmaSquared_double * M_PI;
        for (int jjT=-length; jjT<length+1; jjT++){
           int Y= jjT + length;
           int YnxK = Y * nxK;
           double jjT_squared=pow(jjT,2.0);
           for(int iiT=-length; iiT<length+1; iiT++){
                counter++;
                double distanceSquared =  jjT_squared+pow(iiT,2.0);
                int X=iiT+length;
                kernel [  X + YnxK ] = (1.0 / sigmaSquared_doublePI) * exp (- (distanceSquared) / (sigmaSquared_double));
                sum+=kernel [  X + YnxK ];
           }
        }
        //correct compensate for discretization (sum needs to be 1):
        double normSum = (double((1.0-sum)))/double(kernelSize);
        for (int ii=0; ii<kernelSize; ii++){
                kernel [ii] += normSum;
        }
}


template <typename T>
void gaussianBlur2D(T* I, T* maskI, T*imageOut, double sigma, int kernelLength , const unsigned long int nx, const unsigned long int ny, double maskThreshold=0.9){

const int nxy = nx * ny;
const int length=kernelLength;
const int nxK =  (2*length+1);
const int nyK = (2*length+1);
const int kernelSize = nxK*nyK;
float * kernel = new float [ kernelSize ];
float * tmpI = new float [ nxy ];
double angpix = 1;

fillGaussianKernel(kernel, kernelLength, sigma*sigma);


for ( unsigned long int YY=kernelLength; YY<ny-kernelLength; YY++){
 unsigned long int nxYY=nx*YY;
 for ( unsigned long int XX=kernelLength; XX<nx-kernelLength; XX++){
   unsigned long int idxC=XX+nxYY;
   if ( maskI[idxC] > maskThreshold ){
    double sum = 0;
    for ( int jjT=-kernelLength; jjT<kernelLength; jjT++){
      unsigned long int yImage =(YY+jjT)*nx;
      unsigned long int yKernel = kernelLength + (jjT+kernelLength)*nxK;
      for(int iiT=-kernelLength; iiT<kernelLength; iiT++){
          unsigned long int idx = (XX+iiT)+yImage;
          unsigned long int kernelIdx=iiT+yKernel;
          sum+= I[idx] * kernel[kernelIdx];
          //std::cerr<<kernel[kernelIdx]<<" ";
      }
    }
    //std::cerr<< " ->  " <<sum<<"\n";
    tmpI[idxC]=sum;
   }
 }
}


 for (unsigned long int ij=0; ij<nxy; ij++){
   if ( maskI[ij] > maskThreshold ){
      imageOut[ij]= tmpI[ij];
   }else{
      imageOut[ij]= I[ij];
   }
 }

delete [] kernel;
delete [] tmpI;
}



template <typename T>
void stackGaussianBlur2D(T* I, T* maskI, T* stackOut, std::vector<double> sigmaVector,  const unsigned long int nx, const unsigned long int ny, double maskThreshold=0.9){
  unsigned long int nxy = nx * ny;
  for (unsigned long int ii=0; ii<sigmaVector.size(); ii++ ){
    //std::cerr<<"sigmaVector[ii]="<<sigmaVector[ii]<<"\n";
    if (sigmaVector[ii]<0.01){
      for (unsigned long int ij=0; ij<nxy; ij++ ){
        stackOut[ ij + ii * nxy ] = I[ij];
      }
    }else{
       //gaussianBlur2D(I, maskI, &stackOut[ ii * nxy ], sigmaVector[ii], ceil(sigmaVector[ii]) , nx, ny, maskThreshold);
       GaussianBlurringFFT(I, &stackOut[ ii * nxy ], sigmaVector[ii], nx,  ny, 1);
    }
  }
}


// ***********************************
//  
// 
template <typename T>
void centerOfMassBasedRegularization(T* I, T* maskI, T*imageOut, double sigma, int kernelLength , const unsigned long int nx, const unsigned long int ny){

        const int nxy = nx * ny;
        const int length=kernelLength;
        const int sigmaSquared = sigma*sigma;
        const int nxK =  (2*length+1);
        const int nyK = (2*length+1);
        const int kernelSize = nxK*nyK;
        float * kernel = new float [ kernelSize ];
        float * tmpI = new float [ nxy ];
        double angpix = 1;

        double cmX = 0;
        double cmY = 0;

        unsigned long int counter = 0;
        unsigned long int minX = nx-length-1;
        unsigned long int maxX = length+1;
        unsigned long int minY = ny-length-1;
        unsigned long int maxY = length+1;
        for (unsigned long int YY=length+1; YY<ny-length-1; YY++){
          unsigned long int YYnx = YY * nx;
          for(unsigned long int XX=length+1; XX<nx-length-1; XX++){
                unsigned long int idx = XX + YYnx;
                if ( maskI[idx] > 0.5 ){
                        counter++;
                        cmX += XX;
                        cmY += YY;
                        if (XX < minX) minX = XX;
                        if (XX > maxX) maxX = XX;
                        if (YY < minY) minY = YY;
                        if (YY > maxY) maxY = YY;

                }
          }
        }
        if ( counter<1 || maxX <= minX || maxY <= minY ){
           return;
        }
        cmX /= counter;
        cmY /= counter;
        //std::cerr<<"CM=["<<cmX<<","<<cmY<<"]\n";
        double maskSquaredDistance = 0.2*(pow(maxX-minX,2.0)+pow(maxY-minY,2.0));


        for (unsigned long int YY=length+1; YY<ny-length-1; YY++){
           unsigned long int YYnx = YY * nx;
           for(unsigned long int XX=length+1; XX<nx-length-1; XX++){     
                unsigned long int idx = XX + YYnx;
                if ( maskI[idx] > 0.5 ){
                  //distance from CM
                  double distanceCM=pow(XX-cmX,2.0)+pow(YY-cmY,2.0)+0.0001;
                  distanceCM /= maskSquaredDistance;
                  if ( distanceCM > 0.05 ){
                        fillGaussianKernel ( kernel, length, sigmaSquared*distanceCM  );
                        double sum = 0;
                        for (int jjT=-length; jjT<length; jjT++){
                                unsigned long int jjT_length_nxK = (jjT+length)*nxK;
                                unsigned long int jjT_YY_nx = (YY+jjT)*nx;
                                for(int iiT=-length; iiT<length; iiT++){
                                        unsigned long int idx = (XX+iiT)+jjT_YY_nx;
                                        unsigned long int kernelIdx=iiT+length+jjT_length_nxK;
                                        sum+= kernel[kernelIdx]*I[idx];
                                }
                        }
                        tmpI [ XX + YYnx ] = sum;
                        //tmpI [ XX + YYnx ] = distanceCM;
                        //tmpI [ XX + YYnx ] = 1;
                        }else{
                        tmpI [ XX + YY * nx ] = I [ XX + YYnx ];
                        //tmpI [ XX + YYnx ] = 0;
                        }
                } else {
                        tmpI [ XX + YY * nx ] = I [ XX + YY * nx ];
                        //tmpI [ XX + YY * nx ] = 0;
                }
            }
        }

 for (unsigned long int ii=0;ii<nxy;ii++){
   imageOut[ii]= tmpI[ii];
 }

delete [] kernel;
delete [] tmpI;
}


// ***********************************
//  
// 
template<typename T>
void shiftImage2d( T* I,
                          T* outI,
                          const int xShift,
                          const int yShift,
                          const unsigned long int nx,
                          const unsigned long int ny
){
  unsigned long int nxy=nx*ny;
 T * tmpI = new T[nxy];
 for (unsigned long int ii=0;ii<nxy;ii++){
   tmpI[ii]= I[ii];
   outI[ii]= I[ii];
 }
 int absShiftX=abs(xShift);
 int absShiftY=abs(yShift);
 for (unsigned long int jj=absShiftY; jj<ny-absShiftY; jj++){
  long int nxjjShift=nx*(jj+yShift);
  long int nxjj=nx*(jj);
  for (unsigned long int ii=absShiftX; ii<nx-absShiftX; ii++){
    outI[ii+nxjj]= tmpI[ii-xShift+nxjjShift];
  }
 }

 delete [] tmpI;
 
}




std::vector<double> msaBlur1D (std::vector<double> inputMSA, double sigmaBlur, int padding = 0){

    if ( padding > inputMSA.size() ){
      padding=inputMSA.size();
    }
    unsigned long int unpaddedNx=inputMSA.size();
    unsigned long int nx=2*padding+inputMSA.size();
    std::vector<double> paddedInputMSA;
    for (unsigned long int ii=0;ii<padding;ii++){
      paddedInputMSA.push_back(inputMSA[padding-ii-1]);
    }
    for (unsigned long int ii=0;ii<unpaddedNx;ii++){
      paddedInputMSA.push_back(inputMSA[ii]);
    }
    for (unsigned long int ii=0;ii<padding;ii++){
      paddedInputMSA.push_back(inputMSA[padding-ii-1]);
    }
    
std::cerr<<"\n\n   AA \n\n";

    
    double meanIn=0;
    for (unsigned long int i = 0; i < unpaddedNx; i++){
    					meanIn += inputMSA[i];
    }
    meanIn/=unpaddedNx;

    fftw_complex *I, *F;
		I = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nx);
    F = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nx);
    //inputMSA[0]=1;
    for (unsigned long int i = 0; i < nx; i++){
					I[i][0] = (paddedInputMSA[i]-meanIn);
					I[i][1] = 0.0; //no imaginary part
    }
		fftw_plan planForward;
    planForward = fftw_plan_dft_1d(nx, I, F, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(planForward);
		fftw_destroy_plan(planForward);

    //perform gaussian blur
    for (unsigned long int i=0; i < nx; i++){
      F[ i ][0] /= (nx);
      F[ i ][1] /= (nx);
    }
    //std::cerr<<"Gaussian=";
    for (unsigned long int i=0; i < ceil(nx/2); i++){
      double gaussianValue = exp(-(i*i*sigmaBlur*sigmaBlur)/(nx*4));
    //   std::cerr<<gaussianValue<<",";
      F[ i ][0] *= gaussianValue;
      F[ nx-i ][1] *= gaussianValue;
    }
    //std::cerr<<"\n";

    //put back the transformed vector
    fftw_plan planBackward;
    planBackward = fftw_plan_dft_1d(nx, F, I, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(planBackward);
		fftw_destroy_plan(planBackward);

    std::vector<double> outputMSA;
    std::cerr<<"unpaddedNx="<<unpaddedNx<<"\n";
    for (unsigned long int i = padding; i < padding+unpaddedNx; i++){
          double value=meanIn+(I[i][0]);
          outputMSA.push_back(value);
    }

		fftw_free(I);
		fftw_free(F);

    std::cerr<<" outputMSA.size()="<<outputMSA.size()<<"\n";
    return outputMSA;

}




template<typename T, typename U>
void NormalizeVarianceImage (T * input, //!< input image
                     U * output, //!< output image
                     double sigma, //!< sigma
                     unsigned long int nx, //!< nx
                     unsigned long int ny, //!< nys
                     unsigned long int nz, //!< nz
                     T * maskI = (T*)NULL
 ){
 
//    if (maskI!=NULL){
//    }
  
  unsigned long int nxyz=nx*ny*nz;
  double meanIn = 0;
  double variance = 0;
  
  if (maskI==NULL){
          for (unsigned long int ii = 0; ii<nxyz; ii++){
           meanIn+=input[ii];
          }
          meanIn/=nxyz;
          for (unsigned long int ii = 0; ii<nxyz; ii++){
           variance+=pow(input[ii]-meanIn,2.0);
          }
          variance/=(nxyz-1);
          double SD = pow (variance, 0.5);
          if ( variance != 0 ){
           for (unsigned long int ii = 0; ii<nxyz; ii++){
            output[ii]=meanIn+sigma*(input[ii]-meanIn)/SD;
           }
          }
  }else{
          unsigned long int counter = 0;
          for (unsigned long int ii = 0; ii<nxyz; ii++){
           if (maskI[ii]>0){
             meanIn+=input[ii];
             counter++;
           }
          }
          if (counter > 0 ){
            meanIn/=counter;
          }

          for (unsigned long int ii = 0; ii<nxyz; ii++){
           if (maskI[ii]>0){
            variance+=pow(input[ii]-meanIn,2.0);
           }
          }
          if (counter > 1 ){
            variance/=(counter-1);
          }
          double SD = pow (variance, 0.5);
          if ( variance != 0 ){
           for (unsigned long int ii = 0; ii<nxyz; ii++){
            if (maskI[ii]>0){
              output[ii]=meanIn+sigma*(input[ii]-meanIn)/SD;
            }else{
             output[ii]=meanIn;
            }
           }
          }
  }
  std::cerr<<"MeanIn="<<meanIn<<"\n";
  std::cerr<<"varianceIn="<<variance<<"\n";

}


template<typename T>
void medianFilter (  const T * input, //!< input image
                     T * output, //!< output image
                     unsigned long int nx, //!< nx
                     unsigned long int ny, //!< ny
                     unsigned long int nz, //!< nz
                     int radiusX,
                     int radiusY,
                     int radiusZ,
                     T * mask = NULL
 ){

     unsigned long int nxyz = nx *ny *nz;
     unsigned long int nxy = nx *ny;
     T * tmpI = new T [nxyz];
     for (unsigned long int kk = 0, ijk=0; kk < nz; kk++ ){
        for (unsigned long int jj = 0; jj < ny; jj++ ){
           for (unsigned long int ii = 0; ii < nx; ii++, ijk++){
             long int startX = ii-radiusX; if (startX < 0) startX=0;
             long int startY = jj-radiusY; if (startY < 0) startY=0;
             long int startZ = kk-radiusZ; if (startZ < 0) startZ=0;
             long int endX = ii+radiusX;  if (endX > nx-1) endX=nx-1;
             long int endY = jj+radiusY;  if (endY > ny-1) endY=ny-1;
             long int endZ = kk+radiusZ;  if (endZ > nz-1) endZ=nz-1;
             std::vector<T> arrayValues;
             //std::cerr<<"("<<ii<<","<<jj<<","<<kk<<")   "<<"("<<startX<<"-"<<endX<<",  "<<startY<<"-"<<endY<<",  "<<startZ<<"-"<<endZ<<")   \n";
             for (unsigned long int kkR = startZ; kkR <= endZ; kkR++ ){
                unsigned long int indexK=nxy*kkR;
                for (unsigned long int jjR = startY; jjR <= endY; jjR++ ){
                   unsigned long int indexJK= nx*jjR+indexK;
                   for (unsigned long int iiR = startX; iiR <= endX; iiR++){
                      arrayValues.push_back(  input[iiR+indexJK]  );
                   }
                }
             }
             const int medianIndex = floor( arrayValues.size() / 2.0 );
             std::sort(arrayValues.begin(), arrayValues.end());
             tmpI[ijk]=arrayValues[medianIndex];
             if ( tmpI[ijk] != tmpI[ijk] ){
               tmpI[ijk]=input[ijk];
             }
           }
        }
     }
     
     for (unsigned long int ii = 0; ii < nxyz; ii++){
      output[ii]=tmpI[ii];
     }
     if (tmpI) delete [] tmpI;
}


template<typename T>
std::vector<double> bandpassFilter (  T * input, //!< input image
                     T * output, //!< output image
                     unsigned long int nx, //!< nx
                     unsigned long int ny, //!< ny
                     unsigned long int nz, //!< nz
                     double lowerBand,
                     double higherBand,
                     double sigmaBandBlur,
                     T * mask
 ){

     const unsigned long int nxyz = nx *ny *nz;
     //unsigned long int nxy = nx *ny;
     //T * tmpI = new T [nxyz];
     std::vector<double> msaOriginal=computeMsa (input, (T *) NULL,   nx,   ny,  nz);
     std::vector<double> processingCurve;
     for (unsigned long int ii=0; ii<msaOriginal.size();ii++){
        double value=0;
        double currentV = ((double)ii)/((double)msaOriginal.size()-1.0f);
        //std::cerr<<currentV<<"->";
        if (currentV<lowerBand){
          value=0.0f;
        }else if (currentV > higherBand){
           value=0.0f;
        }else{
         value=1.0f;
        }
        //std::cerr<<value<<"   ";
        processingCurve.push_back(value);
     }
     std::cerr<<"\nlowerBand="<<lowerBand<<"    higherBand="<<higherBand<<"\n    -- \n";
     //std::vector<double> msaBlur;
     //msaBlur=computeGaussianFunctionMSA(msaOriginal, msaBlur, sigma, nx,  ny, nz);
          std::cerr<<"AA  0001 \n";
     
     if (sigmaBandBlur>0){
        processingCurve=msaBlur1D (processingCurve, sigmaBandBlur, msaOriginal.size());
     }
          std::cerr<<"AA  0002 \n";
std::cerr<<"msaOriginal.size()="<<msaOriginal.size()<<" ";
std::cerr<<"processingCurve.size()="<<processingCurve.size()<<"\n"; 
     for (unsigned long int ii=0; ii<msaOriginal.size();ii++){
      msaOriginal[ii]*=processingCurve[ii];
      //std::cerr<<msaOriginal[ii]<<" ";
     }
          std::cerr<<"AA  0003 \n";
     applylMsa ( input, output, (T *) NULL, msaOriginal, nx,  ny, nz );

     return msaOriginal;
     //for (unsigned long int ii = 0; ii < nxyz; ii++){
      //output[ii]=input[ii];
     //}
     //if (tmpI) delete [] tmpI;
}



#endif
