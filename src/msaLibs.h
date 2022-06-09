/** @file src/msaLibs.h
*  @author   Mauro Maiorca, 2017, The Francis Crick Institute
*  @copyright Mauro Maiorca, 2017, The Francis Crick Institute.
*     This program is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*    You should have received a copy of the GNU General Public License
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
*/

/*
    (C) Mauro Maiorca, 2017, The Francis Crick Institute

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __MSA_LIBS__H___
#define __MSA_LIBS__H___



#include <math.h>
#include <iostream>
#include <iomanip>      // std::setprecision
#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cstdlib>


#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileReader.h"
#include "itkMRCImageIO.h"
#include "itkImageFileWriter.h"
#include "itkRecursiveGaussianImageFilter.h"
#define PI 3.14159265358979323846
#define MASK_THRESHOLD_VALUE 0.1f
#include "fftw3.h"
#include "msaLibs.h"
//#include "plot2d.h"

#define SHELL_WIDTH ((double)1.3)
#define EPSILON_ZERO ((double)0.0000000001)
#define __sizeMAV__  floor(  pow( floor( floor(pow(nz/2.0-1, 2))+ floor(pow(ny/2.0-1, 2))+ floor(pow(nx/2.0-1, 2)) ), 0.5 ) / SHELL_WIDTH ) ; //OK
#define __sizeMAV__2D__  floor(  pow( floor( floor(pow(ny/2.0-1, 2))+ floor(pow(nx/2.0-1, 2)) ), 0.5 )  ) ; //OK



// **********************
// **********************
//  I/O
//
// WRITE A VECTOR on a image
//
// **********************
// **********************
template<typename T>
void writeVectorMrcImage(const char * filenameMRC, T * Array, const unsigned long int nx, const unsigned long int ny, const unsigned long int nz){

 typedef float WorkingPixelType;
 typedef float OutputPixelType;
 const unsigned int Dimension = 3;
 typedef itk::Image< WorkingPixelType, Dimension > workingImageType;
 typedef itk::Image< OutputPixelType, Dimension > outputImageType;
 typedef itk::ImageRegionIterator< workingImageType> IteratorType;



      workingImageType::Pointer workingImage = workingImageType::New();
      workingImageType::IndexType workingImageIndex;
      workingImageIndex[0] = 0;
      workingImageIndex[1] = 0;
      workingImageIndex[2] = 0;
      workingImageType::SizeType workingImageSize;
      workingImageSize[0]=nx;
      workingImageSize[1]=ny;
      workingImageSize[2]=nz;
      workingImageType::SpacingType spacing;
        spacing[0] = 1; // spacing
      spacing[1] = 1; // spacing
      spacing[2] = 1; // spacing
      workingImage -> SetSpacing( spacing );
      workingImageType::RegionType workingImageRegion;
      workingImageRegion.SetIndex( workingImageIndex );
      workingImageRegion.SetSize( workingImageSize );
      workingImage->SetRegions( workingImageRegion );
      workingImage->Allocate();

    IteratorType  underProcessingIt( workingImage , workingImageRegion );
    unsigned long int kk=0;
    for (underProcessingIt.GoToBegin(); !underProcessingIt.IsAtEnd(); ++underProcessingIt){
        underProcessingIt.Set( (float)Array[kk]);
        kk++;
    }
    workingImage->Update();

    typedef itk::ImageFileWriter< outputImageType > WriterType;
    WriterType::Pointer Writer = WriterType::New();
    Writer->SetFileName(filenameMRC);
    Writer->SetImageIO( itk::MRCImageIO::New() );
    Writer->SetInput( workingImage );
    Writer->Update();
}


// *************************************
/** Recursive Gaussian Blurring (Uses ITK).
 *
 */
// **************************************
template<typename T, typename U>
void GaussianBlurring(T * input, //!< input image
                      U * output, //!< output image
                      double sigma, //!< sigma
                      const unsigned long int nx, //!< nx
                      const unsigned long int ny, //!< ny
                       const  unsigned long int nz //!< nz
 ){

  typedef float WorkingPixelType;
  typedef float OutputPixelType;
  const unsigned int Dimension = 3;
  typedef itk::Image< WorkingPixelType, Dimension > workingImageType;
  typedef itk::Image< OutputPixelType, Dimension > outputImageType;
  typedef itk::ImageRegionIterator< workingImageType> IteratorType;

	  workingImageType::Pointer workingImage = workingImageType::New();
	  workingImageType::IndexType workingImageIndex;
	  workingImageIndex[0] = 0;
	  workingImageIndex[1] = 0;
	  workingImageIndex[2] = 0;
	  workingImageType::SizeType workingImageSize;
	  workingImageSize[0]=nx;
	  workingImageSize[1]=ny;
	  workingImageSize[2]=nz;
	  workingImageType::SpacingType spacing;
	   spacing[0] = 1; // spacing
		spacing[1] = 1; // spacing
		spacing[2] = 1; // spacing
		workingImage -> SetSpacing( spacing );
	  workingImageType::RegionType workingImageRegion;
	  workingImageRegion.SetIndex( workingImageIndex );
	  workingImageRegion.SetSize( workingImageSize );
	  workingImage->SetRegions( workingImageRegion );
	  workingImage->Allocate();

	 IteratorType workingImageIt(workingImage, workingImageRegion );
	 unsigned long int kk = 0;
     for (workingImageIt.GoToBegin(); !workingImageIt.IsAtEnd(); ++workingImageIt){
		if ( !(input[kk]!=input[kk])){
        	workingImageIt.Set(input[kk]);
       	}else{
       	        	workingImageIt.Set(0);
       	}
        kk++;
    }
	workingImage->Update();


	workingImageType::Pointer blurredImage = workingImageType::New();
	if (sigma <= EPSILON_ZERO){
		 blurredImage = workingImage;
	}else{
		typedef itk::RecursiveGaussianImageFilter <workingImageType, workingImageType >  RecursiveGaussianFilterType;
		RecursiveGaussianFilterType::Pointer filterX = RecursiveGaussianFilterType::New();
		RecursiveGaussianFilterType::Pointer filterY = RecursiveGaussianFilterType::New();
		RecursiveGaussianFilterType::Pointer filterZ = RecursiveGaussianFilterType::New();
        blurredImage = workingImage;
        if ( nx > 4 ){
            filterX->SetDirection( 0 );   // 0 --> X direction; 1 --> Y direction
            filterX->SetOrder( RecursiveGaussianFilterType::ZeroOrder );
            filterX->SetNormalizeAcrossScale( true );
            filterX->SetSigma( sigma );
            filterX->SetInput( workingImage );
            filterX->Update();
            blurredImage = filterX->GetOutput();
        }
        if ( ny > 4 ){
            filterY->SetDirection( 1 );   // 0 --> X direction; 1 --> Y direction
            filterY->SetOrder( RecursiveGaussianFilterType::ZeroOrder );
            filterY->SetNormalizeAcrossScale( true );
            filterY->SetSigma( sigma );
            filterY->SetInput( blurredImage );
            filterY->Update();
            blurredImage = filterY -> GetOutput();
        }
		if ( nz > 4 ){
			filterZ->SetDirection( 2 );   // 0 --> X direction; 1 --> Y direction
			filterZ->SetOrder( RecursiveGaussianFilterType::ZeroOrder );
			filterZ->SetNormalizeAcrossScale( true );
			filterZ->SetSigma( sigma );
			filterZ->SetInput( blurredImage );
			filterZ->Update();
			blurredImage = filterZ -> GetOutput();
		}
	}

	 IteratorType outputImageIt(blurredImage , workingImageRegion );
	 kk=0;
     for (outputImageIt.GoToBegin(); !outputImageIt.IsAtEnd(); ++outputImageIt){
		output[kk]=outputImageIt.Get();
		kk++;
    }


}






// *************************************
/** Recursive Gaussian Blurring (Uses ITK).
 *
 */
// **************************************
template<typename T, typename U>
void GaussianFirstDerivative2D(T * input, //!< input image
                      U * output, //!< output image
                      double sigma, //!< sigma
                      int direction, //!< sigma
                      const unsigned long int nx, //!< nx
                      const unsigned long int ny //!< ny
 ){

  typedef float WorkingPixelType;
  typedef float OutputPixelType;
  const unsigned int Dimension = 2;
  typedef itk::Image< WorkingPixelType, Dimension > workingImageType;
  typedef itk::Image< OutputPixelType, Dimension > outputImageType;
  typedef itk::ImageRegionIterator< workingImageType> IteratorType;

	  workingImageType::Pointer workingImage = workingImageType::New();
	  workingImageType::IndexType workingImageIndex;
	  workingImageIndex[0] = 0;
	  workingImageIndex[1] = 0;
	  workingImageType::SizeType workingImageSize;
	  workingImageSize[0]=nx;
	  workingImageSize[1]=ny;
	  workingImageType::SpacingType spacing;
	   spacing[0] = 1; // spacing
		spacing[1] = 1; // spacing
		workingImage -> SetSpacing( spacing );
	  workingImageType::RegionType workingImageRegion;
	  workingImageRegion.SetIndex( workingImageIndex );
	  workingImageRegion.SetSize( workingImageSize );
	  workingImage->SetRegions( workingImageRegion );
	  workingImage->Allocate();

	 IteratorType workingImageIt(workingImage, workingImageRegion );
	 unsigned long int kk = 0;
     for (workingImageIt.GoToBegin(); !workingImageIt.IsAtEnd(); ++workingImageIt){
		if ( !(input[kk]!=input[kk])){
        	workingImageIt.Set(input[kk]);
       	}else{
       	        	workingImageIt.Set(0);
       	}
        kk++;
    }
	workingImage->Update();


	workingImageType::Pointer blurredImage = workingImageType::New();
	if (sigma <= EPSILON_ZERO){
		 blurredImage = workingImage;
	}else{
		typedef itk::RecursiveGaussianImageFilter <workingImageType, workingImageType >  RecursiveGaussianFilterType;
		RecursiveGaussianFilterType::Pointer filterX = RecursiveGaussianFilterType::New();
		RecursiveGaussianFilterType::Pointer filterY = RecursiveGaussianFilterType::New();
		RecursiveGaussianFilterType::Pointer filterZ = RecursiveGaussianFilterType::New();
        blurredImage = workingImage;
        if ( nx > 4 && ny > 4){
            filterX->SetDirection( direction );   // 0 --> X direction; 1 --> Y direction
            filterX->SetOrder( RecursiveGaussianFilterType::FirstOrder );
            filterX->SetNormalizeAcrossScale( true );
            filterX->SetSigma( sigma );
            filterX->SetInput( workingImage );
            filterX->Update();
            blurredImage = filterX->GetOutput();
        }
	}

	 IteratorType outputImageIt(blurredImage , workingImageRegion );
	 kk=0;
     for (outputImageIt.GoToBegin(); !outputImageIt.IsAtEnd(); ++outputImageIt){
		output[kk]=outputImageIt.Get();
		kk++;
    }
}







// **********************
// **********************
// OTSU
// **********************
// **********************
// *************************************
/** build the grayscale histogram from an image */
// buildHistogram
// **************************************
template<typename T>
void buildHistogram(const T* image, //!< input image
        double * histogram, //!< output histogram
        const unsigned long int bins, //!< number of bins for the histogram
        const unsigned long int nxyz,  //!< size (in pixel) of the image
        const T * mask = NULL,
        double minMaskValue = 0.0
){
  if (mask==NULL){ //NO MASK
	double min=image[0], max=image[0];
	for (int ii = 1; ii< nxyz; ii++){
		if (min>image[ii]) min=image[ii];
		if (max<image[ii]) max=image[ii];
	}
    for (int ii = 0; ii< bins; ii++){
		histogram[ii]=0.0f;
	}
	for (int ii = 0; ii< nxyz; ii++){
		double val = 0;
		if (max-min>0){
		 val = (image[ii] - min)*(bins-1)/(max-min);
		}
		int index = floor(val);
		histogram[index]++;
	}
   }else{ //MASK
	double min=image[0], max=image[0];
	for (int ii = 1; ii< nxyz; ii++){
	   if(mask[ii]>minMaskValue){
		if (min>image[ii]) min=image[ii];
		if (max<image[ii]) max=image[ii];
           }
	}
        for (int ii = 0; ii< bins; ii++){
		histogram[ii]=0.0f;
	}
	for (int ii = 0; ii< nxyz; ii++){
	   if(mask[ii]>minMaskValue){
		double val = 0;
		if (max-min>0){
		 val = (image[ii] - min)*(bins-1)/(max-min);
		}
		int index = floor(val);
		histogram[index]++;
           }
	}
   }
 }


// *************************************
//
/** This function recursively implements a multiple Otsu thresholding algorithm for
 * multiple thresolds. The code is based on the paper:
 * Liao, P-S. & Chung, P-C. (September 2001), "A fast algorithm for multilevel thresholding", Journal of Information Science and Engineering 17 (5): 713-727
 */
//
// **************************************
void computeThresholdsRecursively(const double * H, //!< lookup table for the histogram
    double * t, //!< vector with tcomputed hreshold values
    const unsigned long int bins, //!< histogram bin size
    const unsigned long int NT, //!< number of threshold
    bool firstCall = true, int start=1, int end=1,  unsigned long int n=1, double * maxSig = NULL, double * Sq=NULL, unsigned int * m=NULL){

		if(firstCall){
			start=1;
			end=bins-NT;
			n=1;
			maxSig =new double (0);
			Sq=new double[NT+1];
			m=new unsigned int[NT+2];
		}

		for (unsigned long int i = start; i < end; i++) {
			m[n]=i;
			if(n==1) {
				Sq[0] = H[1+bins*i];
			}
			if(n>1 && n <= NT){
				Sq[n-1] = H[start+bins*i];
			}
			if(n==NT){
				Sq[n] = H[i+1+bins*(bins-1)];
				double SqVal = 0;
				for (int kk=0; kk <= NT; kk++){
					SqVal += Sq[kk];
				}
				//std::cerr<< "Sq="<<SqVal <<"\n";
				if (*maxSig < SqVal)	{
					for (unsigned int cc=0; cc<NT; cc++){
							t[cc] = m[cc+1];
					}
					*maxSig = SqVal;
				}

			}else{
				computeThresholdsRecursively(H,t, bins, NT, false, i+1,end+1, n+1, maxSig, Sq, m);
			}
		}

		if(firstCall){
			delete maxSig;
			delete [] Sq;
			delete [] m;
		}

}



// *************************************
/** Compute Otsu threshold in an image */
// computeOtsuThresholds
// **************************************
template<typename T>
void computeOtsuThresholds ( T * image, //!< input image
    double * t, //!< output vector of threshold values
    const unsigned long int nt, //!< number of threshold values
    const unsigned long int bins, //!< number of bins in the histogram
    const unsigned long int nxyz, //!< size of the image
    bool verbose = false  //!< verbose output
){
		double epsilon = 0.01;

		// ----------------------------
		// COMPUTE HISTOGRAM
		// std::cerr<<" COMPUTE HISTOGRAM\n";
	  	double * histogram = new double [bins];
	  	buildHistogram(image, histogram,  bins, nxyz);
		T minImageValue=image[0], maxImageValue=image[0];
		for (int ii = 1; ii< nxyz; ii++){
			if (minImageValue>image[ii]) minImageValue=image[ii];
			if (maxImageValue<image[ii]) maxImageValue=image[ii];
		}
		for (int ii = 0; ii< bins; ii++){
			histogram[ii]=0.0f;
		}
		double ratioVal = 0;
		if (maxImageValue > minImageValue){
		 ratioVal=((double)((double)bins-1.0f))/((double)((double)maxImageValue-(double)minImageValue));
		}
		for (int ii = 0; ii< nxyz; ii++){
			double val = (image[ii] - minImageValue)*ratioVal;
                        int index = floor(val);
			histogram[index]++;
		}

		//std::cerr<<"Histogram Statistics: min value="<< minImageValue <<"   max value=" <<  maxImageValue <<"\n";


		// ----------------------------
		//BUILD LOOKUP TABLES
		//std::cerr<<" BUILD LOOKUP TABLES\n";
	  	double * P = new double[bins*bins];
	  	double * S = new double[bins*bins];
	  	double * H = new double[bins*bins];
			for (unsigned long int i = 0; i < bins*bins; i++){
				P[i]=S[i]=H[i]=0.0f;
			}
			//diagonal (row 0 is all zero)
			//std::cerr<<" diagonal\n";
			for (unsigned long int i=1; i < bins; i++){
				unsigned long int diagonalIndex = i+i*bins;
				 P[diagonalIndex] = histogram[i];
				 S[diagonalIndex] = ((double) i)*histogram[i];
			}
			// calculate first row (row 0 is all zero)
			//std::cerr<<" calculate first row\n";
			unsigned long int firstRowIndex = bins;
			for (unsigned long int i=1; i < bins-1; i++) {
			  P [1+(i+1)*bins] = P[1+i*bins] + histogram[i+1];
			  S [1+(i+1)*bins] = S[1+i*bins] + ((double) (i+1))*histogram[i+1];
			}
			// using row 1 to calculate others
			//std::cerr<<" using row 1 to calculate others\n";
			for (unsigned long int i=2; i < bins; i++){
			  for (unsigned long int j=i+1; j < bins; j++) {
				P[i+j*bins] = P[1+j*bins] - P[1+(i-1)*bins];
				S[i+j*bins] = S[1+j*bins] - S[1+(i-1)*bins];
			  }
			}
			// calculate H
			//std::cerr<<" calculate H\n";
			for (unsigned long int i=1; i < bins; i++){
			  for (unsigned long int j=i+1; j < bins; j++){
			  		unsigned long int ijIndex = i+j*bins;
					if (pow (P[ ijIndex ],2) > epsilon)
					  H[ ijIndex ] = (S[ ijIndex ]*S[ ijIndex ])/P[ ijIndex ];
					else
					  H[ ijIndex ] = 0.0f;
			  }
			}
	  	delete [] P;
	  	delete [] S;


		// ----------------------------
		// COMPUTE THRESHOLD VALUES
		if (verbose) std::cerr<<"[COMPUTE THRESHOLD VALUES...";
		computeThresholdsRecursively(H, t, bins, nt);
		//computeThresholdsIteratively(H, t, bins, nt);
		if (verbose) std::cerr<<" DONE]\n";
	  	delete [] H;
	  	delete [] histogram;

		//std::cerr<<"Binned Thresholds=[";
		//for (int i =0; i<nt; i++){
			//std::cerr<<t[i];
			//if (i<nt-1)
				//std::cerr<<",";
		//}
		//std::cerr<<"]\n";

		// ----------------------------
		// NORMALIZE HISTOGRAM VALUES
		// TO ORIGINAL IMAGE VALUES
		//std::cerr<<"[NORMALIZE HISTOGRAM VALUES...";
		double ratioValInv = 0;
		if (maxImageValue > minImageValue && bins>1){
		         ratioValInv=(double)((double)maxImageValue-(double)minImageValue)/(double)((double)bins-1.0f);
		}
		for(unsigned int tt = 0; tt < nt; tt++ ){
			        t [ tt ] = minImageValue + t [ tt ]*ratioValInv;
		}
		//std::cerr<<" DONE]\n";
}

template<typename T>
std::vector<double> computeOtsuThresholdsVector ( std::vector<T> values, //!< input image
    const unsigned long int nt, //!< number of threshold values
    const unsigned long int bins //!< number of bins in the histogram
){
  T * image = new T [values.size()];
  for (unsigned long int ii=0;ii<values.size();ii++){
   image[ii]=values[ii];
  }
  double * t=new double[nt];
  computeOtsuThresholds ( image, t, nt, bins, values.size() );
  std::vector<double> vectorOut;
  for (unsigned long int ii=0;ii<nt;ii++){
   vectorOut.push_back( t[ii] );
  }
  delete [] t;
  delete [] image;
  return vectorOut;
}

// *************************************
//
/** Main parcellation procedure using multiple Ostu threshold. The pixel value of each
 * parcellated region is the mean value of that region in the original image.
 */
// computeOtsuMeanImage
// **************************************
template<typename T>
void computeOtsuMeanImage ( T * image, //!< Original Input Image
  T * outImage, //!< Output image
  const unsigned long int numberOfThresholds, //!< number of threshold values
 const unsigned long int bins, //!< number of bins in the histogram
 double sigmaSmooth, //!< Sigma value for the Gaussian blurring the original image (regularization useful for the Otsu method)
 unsigned long int nx, //!< size of image (X)
 unsigned long int ny, //!< size of image (Y)
 unsigned long int nz , //!< size of image (Z)
T * preprocessedImage = NULL, //!< Image to use for computing the otsu (intead of using the original image)
bool verbose = false   //!< check if verbose
){
    if (nx==0) return;
    if (ny==0) return;
    if (nz==0) return;
    
	const unsigned long int nxyz = nx*ny*nz;
	const double epsilon = 0.01;
	T * blurredImage = new T [nxyz];

	if ( sigmaSmooth<epsilon ){
		sigmaSmooth = epsilon;
	}


	if ( preprocessedImage ){
		if (sigmaSmooth > 0.0f){
			GaussianBlurring(preprocessedImage, blurredImage, sigmaSmooth,  nx,  ny,  nz);
		}else{
			for (unsigned long int ii = 0; ii< nxyz; ii++){
				blurredImage[ii]=preprocessedImage[ii];
			}
		}
	}else{
		if (sigmaSmooth > 0.0f){
			GaussianBlurring(image, blurredImage, sigmaSmooth,  nx,  ny,  nz);
		}else{
			for (unsigned long int ii = 0; ii< nxyz; ii++){
				blurredImage[ii]=image[ii];
			}
		}
	}


	double * t = new double [numberOfThresholds];
	computeOtsuThresholds( blurredImage, t,  numberOfThresholds,  bins, nxyz);

	double * countsClasses = new double [numberOfThresholds + 1];
	double * meansClasses = new double [numberOfThresholds + 1];
	for (int jj=0; jj<numberOfThresholds+1; jj++){
		countsClasses[jj]=0.0f;
		meansClasses[jj]=0.0f;
	}
	
	if (verbose){
        std::cerr<<"Threshold Values=[";
        for (int jj=0; jj<numberOfThresholds; jj++){
                std::cerr<<t[jj];
            if (jj<numberOfThresholds-1)
                std::cerr<<",";
        }
        std::cerr<<"]\n";
    }
	//compute mean for each threshold
	for (unsigned long int ii=0; ii<nxyz; ii++){
		unsigned int index = 0;
		bool found = false;
		for (int jj=0; jj<numberOfThresholds && !found; jj++){
			if( blurredImage[ii] < t[jj] ){
				index = jj+1;
				found = true;
			}
		}
		countsClasses[index]++;
		meansClasses[index]+=image[ii];
		outImage[ii]=index;

	}

	for (int jj=0; jj<numberOfThresholds+1; jj++){
		if (countsClasses[jj]>epsilon)
			meansClasses[jj]/=countsClasses[jj];
	}


	for (unsigned long int ii=0; ii<nxyz; ii++){
		unsigned int index = (unsigned int) outImage[ii];
		outImage[ii]=meansClasses[index];
	}

    if (verbose){
        std::cerr<<"Mean Values=[";
        for (int jj=0; jj<=numberOfThresholds; jj++){
                std::cerr<<meansClasses[jj];
            if (jj<numberOfThresholds)
                std::cerr<<",";
        }
        std::cerr<<"]\n";
    }
	delete [] countsClasses;
	delete [] t;
	delete [] meansClasses;
	delete [] blurredImage;

}



/** 
 * Split the Otsu mean image
 */
// computeOtsuMeanImage
// **************************************
template<typename T>
void computeOtsuMeanImageSplit ( T * image, //!< Original Input Image
 T * outImage, //!< Output image
 const unsigned long int numberOfThresholds, //!< number of threshold values
 const unsigned long int bins, //!< number of bins in the histogram
 double sigmaSmooth, //!< Sigma value for the Gaussian blurring the original image (regularization useful for the Otsu method)
 int regionSizeInPixel, //!< Size of the region in Pixel
 const unsigned long int nx, //!< size of image (X)
 const unsigned long int ny, //!< size of image (Y)
 const unsigned long int nz , //!< size of image (Z)
 T * preprocessedImage = NULL //!< Image to use for computing the otsu (intead of using the original image)
){

    const unsigned long int nxy = nx * ny;
    const unsigned long int nxyz = nx * ny * nz;
    
    //std::cerr<<"\n\n  regionSizeInPixel="<<regionSizeInPixel<<"\n";
    //int regionSizeInPixel = 10;
    unsigned int nRx = ceil(nx / regionSizeInPixel);
    unsigned int nRy = ceil(ny / regionSizeInPixel);
    unsigned int nRz = ceil(nz / regionSizeInPixel);
    
    if ( nx-nRx*regionSizeInPixel <=0 ) nRx--;
    if ( ny-nRy*regionSizeInPixel <=0 ) nRy--;
    if ( nz-nRz*regionSizeInPixel <=0 ) nRz--;
    
    
    T * tmpImage = new T [regionSizeInPixel*regionSizeInPixel*regionSizeInPixel];
    T * tmpImageOut = new T [regionSizeInPixel*regionSizeInPixel*regionSizeInPixel];
    int counter = 0;
    for ( unsigned long int rZ = 0; rZ <= nRz; rZ++ ){
        unsigned long int Z1 = rZ*regionSizeInPixel;
        unsigned long int Z2 = (rZ+1)*regionSizeInPixel;
        unsigned long int Z3 = (rZ+2)*regionSizeInPixel;
        if (Z2 > nz) Z2 = nz;
        for ( unsigned long int rY = 0; rY <= nRy; rY++ ){
            unsigned long int Y1 = rY*regionSizeInPixel;
            unsigned long int Y2 = (rY+1)*regionSizeInPixel;
            unsigned long int Y3 = (rY+2)*regionSizeInPixel;
            if (Y2 > ny) Y2 = ny;
            for ( unsigned long int rX = 0; rX <= nRx; rX++ ){
                unsigned long int X1 = rX*regionSizeInPixel;
                unsigned long int X2 = (rX+1)*regionSizeInPixel;
                unsigned long int X3 = (rX+2)*regionSizeInPixel;
                if (X2 > nx) X2 = nx;
                
                //std::cerr<< counter << " X=(" << X1 << "," << X2 << ")  Y=(" <<Y1<< "," << Y2 << ")  Z=(" << Z1 << ", " <<Z2<< ")\n";

                for (unsigned long int kk=Z1, tmpKK=0; kk<Z2; kk++, tmpKK++){
                    for (unsigned long int jj=Y1, tmpJJ=0; jj<Y2; jj++, tmpJJ++){
                        for (unsigned long int ii=X1, tmpII=0; ii<X2; ii++, tmpII++){
                            tmpImage[tmpII+tmpJJ*(X2-X1)+tmpKK*(X2-X1)*(Y2-Y1)]=image[ii+jj*nx+kk*nxy];
                        }
                    }
                }
                
                computeOtsuMeanImage ( tmpImage, tmpImageOut, numberOfThresholds, bins, sigmaSmooth, (X2-X1), (Y2-Y1), (Z2-Z1));
                
                for (unsigned long int kk=Z1, tmpKK=0; kk<Z2; kk++, tmpKK++){
                    for (unsigned long int jj=Y1, tmpJJ=0; jj<Y2; jj++, tmpJJ++){
                        for (unsigned long int ii=X1, tmpII=0; ii<X2; ii++, tmpII++){
                            outImage[ii+jj*nx+kk*nxy]=tmpImageOut[tmpII+tmpJJ*(X2-X1)+tmpKK*(X2-X1)*(Y2-Y1)];
                        }
                    }
                }                

                counter ++;
                
            }
        }
    }

    
    delete [] tmpImage;
    delete [] tmpImageOut;
    
//    computeOtsuMeanImage ( image, outImage, numberOfThresholds, bins, sigmaSmooth, nx, ny, nz, preprocessedImage);
//    std::cerr<<"\n\n TEST \n\n";
    
}

// ////////////////////////////////
// MEDIAN FILTER
template<typename T>
void medianFilter (T * input, //!< input image
                     T * output, //!< output image
                     unsigned long int nx, //!< nx
                     unsigned long int ny, //!< ny
                     unsigned long int nz, //!< nz
                     int radiusX,
                     int radiusY,
                     int radiusZ
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



// *********************************************
// GENERIC FFT FUNCTION
// *********************************************

// *********************************************
// Replace amplitude

template<typename T>
void  replaceAmplitude (const T * mapIn, const T * replacementVect, const T * maskVect, T * mapOut, const unsigned long int nx,  const unsigned long int ny, const  unsigned long int nz){
    const double small_epsilon = 0.0001;
    const unsigned long int nxy = nx * ny;
    const unsigned long int nxyz = nx * ny * nz;
    double nD = 0.0f;
    if (nx > 1 ) nD++;
    if (ny > 1 ) nD++;
    if (nz > 1 ) nD++;

    //compute mean
    double meanIn = 0;
    double meanReplacement=0;
    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
      meanIn+=mapIn[ijk];
      meanReplacement+=replacementVect[ijk];
    }
    meanIn/=nxyz;
    meanReplacement/=nxyz;


	fftw_complex *I, *F, *O;
		I = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
		F = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
        O = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
    

    //REPLACEMENT VECTOR MAP
    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
					I[ijk][0] = (replacementVect[nxyz-ijk-1]-meanReplacement);
					I[ijk][1] = 0.0; //no imaginary part
    }
    fftw_plan planForwardMap;    
    if (nz>1){
		  planForwardMap = fftw_plan_dft_3d(nz, ny, nx, I, F, FFTW_FORWARD, FFTW_ESTIMATE);
    }else{
      planForwardMap = fftw_plan_dft_2d(ny, nx, I, F, FFTW_FORWARD, FFTW_ESTIMATE);
    }
    fftw_execute(planForwardMap);
    fftw_destroy_plan(planForwardMap);
    for (unsigned long int zyxIndex=0; zyxIndex < nxyz; zyxIndex++){
        O[ zyxIndex ][0] = F[ zyxIndex ][0]/nxyz;
        O[ zyxIndex ][1] = F[ zyxIndex ][1]/nxyz;
    }


    //MAP_IN
    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
					I[ijk][0] = (mapIn[nxyz-ijk-1]-meanIn);
					I[ijk][1] = 0.0; //no imaginary part
    }
    if (nz>1){
		  planForwardMap = fftw_plan_dft_3d(nz, ny, nx, I, F, FFTW_FORWARD, FFTW_ESTIMATE);
    }else{
      planForwardMap = fftw_plan_dft_2d(ny, nx, I, F, FFTW_FORWARD, FFTW_ESTIMATE);
    }
    fftw_execute(planForwardMap);
    fftw_destroy_plan(planForwardMap);
    for (unsigned long int zyxIndex=0; zyxIndex < nxyz; zyxIndex++){
        F[ zyxIndex ][0] = F[ zyxIndex ][0]/nxyz;
        F[ zyxIndex ][1] = F[ zyxIndex ][1]/nxyz;
    }


    for (unsigned long int zyxIndex=0; zyxIndex < nxyz; zyxIndex++){
        double outputAmplitude = pow( O[ zyxIndex ][0] * O[ zyxIndex ][0] + O[ zyxIndex ][1] * O[ zyxIndex ][1],0.5);
        double inputAmplitude =  pow( F[ zyxIndex ][0] * F[ zyxIndex ][0] + F[ zyxIndex ][1] * F[ zyxIndex ][1],0.5);
        if (inputAmplitude < small_epsilon){
            O[ zyxIndex ][0] = F[ zyxIndex ][0]*outputAmplitude/small_epsilon;
            O[ zyxIndex ][1] = F[ zyxIndex ][1]*outputAmplitude/small_epsilon;
        } else{
            O[ zyxIndex ][0] = F[ zyxIndex ][0]*outputAmplitude/inputAmplitude;
            O[ zyxIndex ][1] = F[ zyxIndex ][1]*outputAmplitude/inputAmplitude;
        }
    }    
    
    if (maskVect){
      for (unsigned long int zyxIndex=0; zyxIndex < nxyz; zyxIndex++){
        if (maskVect[ zyxIndex ] <= MASK_THRESHOLD_VALUE){
            O[ zyxIndex ][0] = 0;
            O[ zyxIndex ][1] = 0;
        }
      }
    }
    
    
  fftw_plan planForwardOriginalImage;
  if (nz>1){
    planForwardOriginalImage = fftw_plan_dft_3d(nz, ny, nx, O, F, FFTW_BACKWARD, FFTW_ESTIMATE);
  }else{
    planForwardOriginalImage = fftw_plan_dft_2d(ny, nx, O, F, FFTW_BACKWARD, FFTW_ESTIMATE);
  }

  fftw_execute(planForwardOriginalImage);
  fftw_destroy_plan(planForwardOriginalImage);

  for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
          mapOut[ijk] = meanReplacement+(F[nxyz-ijk-1][0]);
  }
		fftw_free(I);
		fftw_free(F);
		fftw_free(O);

        
     

}

// *********************************************
// Replace amplitude

template<typename T>
void  PSF_correction (const T * mapIn, const T * psfI, T * mapOut, const unsigned long int nx,  const unsigned long int ny, const  unsigned long int nz){
    const double small_epsilon = 0.001;
    const unsigned long int nxy = nx * ny;
    const unsigned long int nxyz = nx * ny * nz;
    double nD = 0.0f;
    if (nx > 1 ) nD++;
    if (ny > 1 ) nD++;
    if (nz > 1 ) nD++;

    //compute mean
    double meanIn = 0;
    double meanReplacement=0;
    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
      meanIn+=mapIn[ijk];
      meanReplacement+=psfI[ijk];
    }
    meanIn/=nxyz;
    meanReplacement/=nxyz;


	fftw_complex *I, *F, *O;
		I = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
		F = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
        O = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
    

    //REPLACEMENT VECTOR MAP
    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
					I[ijk][0] = (psfI[nxyz-ijk-1]-meanReplacement);
					I[ijk][1] = 0.0; //no imaginary part
    }
    fftw_plan planForwardMap;    
    if (nz>1){
		  planForwardMap = fftw_plan_dft_3d(nz, ny, nx, I, F, FFTW_FORWARD, FFTW_ESTIMATE);
    }else{
      planForwardMap = fftw_plan_dft_2d(ny, nx, I, F, FFTW_FORWARD, FFTW_ESTIMATE);
    }
    fftw_execute(planForwardMap);
    fftw_destroy_plan(planForwardMap);
    for (unsigned long int zyxIndex=0; zyxIndex < nxyz; zyxIndex++){
        O[ zyxIndex ][0] = F[ zyxIndex ][0]/nxyz; //PSF
        O[ zyxIndex ][1] = F[ zyxIndex ][1]/nxyz;
    }


    //MAP_IN
    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
					I[ijk][0] = (mapIn[nxyz-ijk-1]-meanIn);
					I[ijk][1] = 0.0; //no imaginary part
    }
    if (nz>1){
		  planForwardMap = fftw_plan_dft_3d(nz, ny, nx, I, F, FFTW_FORWARD, FFTW_ESTIMATE);
    }else{
      planForwardMap = fftw_plan_dft_2d(ny, nx, I, F, FFTW_FORWARD, FFTW_ESTIMATE);
    }
    fftw_execute(planForwardMap);
    fftw_destroy_plan(planForwardMap);
    for (unsigned long int zyxIndex=0; zyxIndex < nxyz; zyxIndex++){
        F[ zyxIndex ][0] = F[ zyxIndex ][0]/nxyz;
        F[ zyxIndex ][1] = F[ zyxIndex ][1]/nxyz;
    }


    for (unsigned long int zyxIndex=0; zyxIndex < nxyz; zyxIndex++){
        double psfAmplitude = pow( O[ zyxIndex ][0] * O[ zyxIndex ][0] + O[ zyxIndex ][1] * O[ zyxIndex ][1],0.5);
        //double inputAmplitude =  pow( F[ zyxIndex ][0] * F[ zyxIndex ][0] + F[ zyxIndex ][1] * F[ zyxIndex ][1],0.5);
	O[ zyxIndex ][0] = F[ zyxIndex ][0]/(psfAmplitude+small_epsilon);
	O[ zyxIndex ][1] = F[ zyxIndex ][1]/(psfAmplitude+small_epsilon);
    }    
    

    
    
  fftw_plan planForwardOriginalImage;
  if (nz>1){
    planForwardOriginalImage = fftw_plan_dft_3d(nz, ny, nx, O, F, FFTW_BACKWARD, FFTW_ESTIMATE);
  }else{
    planForwardOriginalImage = fftw_plan_dft_2d(ny, nx, O, F, FFTW_BACKWARD, FFTW_ESTIMATE);
  }

  fftw_execute(planForwardOriginalImage);
  fftw_destroy_plan(planForwardOriginalImage);

  for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
          mapOut[ijk] = meanReplacement+(F[nxyz-ijk-1][0]);
  }
		fftw_free(I);
		fftw_free(F);
		fftw_free(O);

}



// Fill missing spectrum
template<typename T>
void  fillMissingSpectrum (const T * mapIn, const T * replacementVect, const T * maskVect, T * mapOut, const unsigned long int nx,  const unsigned long int ny, const  unsigned long int nz){
        
    
    const unsigned long int nxy = nx * ny;
    const unsigned long int nxyz = nx * ny * nz;
    double nD = 0.0f;
    if (nx > 1 ) nD++;
    if (ny > 1 ) nD++;
    if (nz > 1 ) nD++;

    //compute mean
    double meanIn = 0;
    double meanReplacement=0;
    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
      meanIn+=mapIn[ijk];
      meanReplacement+=replacementVect[ijk];
    }
    meanIn/=nxyz;
    meanReplacement/=nxyz;


	fftw_complex *I, *F, *O;
		I = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
		F = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
        O = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
					I[ijk][0] = (mapIn[nxyz-ijk-1]-meanIn);
					I[ijk][1] = 0.0; //no imaginary part
    }
    fftw_plan planForwardMap;
    if (nz>1){
		  planForwardMap = fftw_plan_dft_3d(nz, ny, nx, I, F, FFTW_FORWARD, FFTW_ESTIMATE);
    }else{
      planForwardMap = fftw_plan_dft_2d(ny, nx, I, F, FFTW_FORWARD, FFTW_ESTIMATE);
    }

    fftw_execute(planForwardMap);
    fftw_destroy_plan(planForwardMap);

    for (unsigned long int zyxIndex=0; zyxIndex < nxyz; zyxIndex++){
      if (maskVect[ zyxIndex ] > MASK_THRESHOLD_VALUE){
        O[ zyxIndex ][0] = F[ zyxIndex ][0]/nxyz;
        O[ zyxIndex ][1] = F[ zyxIndex ][1]/nxyz;
      }
    }

    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
					I[ijk][0] = (replacementVect[nxyz-ijk-1]-meanReplacement);
					I[ijk][1] = 0.0; //no imaginary part
    }
    //fftw_plan planForwardMap;
    if (nz>1){
		  planForwardMap = fftw_plan_dft_3d(nz, ny, nx, I, F, FFTW_FORWARD, FFTW_ESTIMATE);
    }else{
      planForwardMap = fftw_plan_dft_2d(ny, nx, I, F, FFTW_FORWARD, FFTW_ESTIMATE);
    }

    fftw_execute(planForwardMap);
    fftw_destroy_plan(planForwardMap);

    for (unsigned long int zyxIndex=0; zyxIndex < nxyz; zyxIndex++){
      if (maskVect[ zyxIndex ] <= MASK_THRESHOLD_VALUE){
        O[ zyxIndex ][0] = F[ zyxIndex ][0]/nxyz;
        O[ zyxIndex ][1] = F[ zyxIndex ][1]/nxyz;
      }
    }    
    
    
  fftw_plan planForwardOriginalImage;
 if (nz>1){
    planForwardOriginalImage = fftw_plan_dft_3d(nz, ny, nx, O, F, FFTW_BACKWARD, FFTW_ESTIMATE);
  }else{
    planForwardOriginalImage = fftw_plan_dft_2d(ny, nx, O, F, FFTW_BACKWARD, FFTW_ESTIMATE);
  }

  fftw_execute(planForwardOriginalImage);
  fftw_destroy_plan(planForwardOriginalImage);

    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
          mapOut[ijk] = meanIn+(F[nxyz-ijk-1][0]);
    }
		fftw_free(I);
		fftw_free(F);
		fftw_free(O);

    
}


// *********************************************
// MSA SPECIFIC FUNCTION
// *********************************************

// **********************
// **********************
//
// blur1D vector
//
// **********************
// **********************
std::vector<double> vectorFftBlur1D (std::vector<double> inputMSA, double sigmaBlur){

    unsigned long int nx=inputMSA.size();
    double meanIn=0;
        for (unsigned long int i = 0; i < nx; i++){
    					meanIn += inputMSA[i];
        }
    meanIn/=nx;

    fftw_complex *I, *F;
    I = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nx);
    F = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nx);
    //inputMSA[0]=1;
    for (unsigned long int i = 0; i < nx; i++){
					I[i][0] = (inputMSA[i]-meanIn);
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
    for (unsigned long int i = 0; i < nx; i++){
          double value=meanIn+(I[i][0]);
          outputMSA.push_back(value);
    }

		fftw_free(I);
		fftw_free(F);

    return outputMSA;

}
std::vector<double> msaBlur (std::vector<double> inputMSA, double sigmaBlur){
  return vectorFftBlur1D (inputMSA, sigmaBlur);
}




// **********************
// **********************
//
// computeResolutions
//
// **********************
// **********************
std::vector<double> computeResolutions (const std::vector<double> inputMSA, const double apix, const unsigned long int nx, const unsigned long int ny, const unsigned long int nz){
    const double sizeMAV =  __sizeMAV__;
    unsigned long int numShells = inputMSA.size();
    std::vector<double> listResolutions;
    listResolutions.push_back(0);
    double shellSizeA=SHELL_WIDTH/(2*apix*sizeMAV);
    for (int ii = 1; ii < numShells; ii++){
      double A = shellSizeA*(ii+1); //OK
      listResolutions.push_back(A);
    }
    return listResolutions;
 }

// **********************
// **********************
//
// msaCutoff
//
// **********************
// **********************
void msaCutoff (std::vector<double> & inputMSA, const unsigned long int nx, const unsigned long int ny, const unsigned long int nz){
//void msaSoftCutoff (std::vector<double> & inputMSA, const unsigned long int nx, const unsigned long int ny, const unsigned long int nz){

    const double sizeMAV = __sizeMAV__;
    unsigned long int numShells = inputMSA.size();
    //std::vector<double> msaCutoffVector=inputMSA;
    for (unsigned long int ShellIndex = 0; ShellIndex < numShells; ShellIndex++){
          if (inputMSA[ShellIndex] <= EPSILON_ZERO || ShellIndex*sqrt(2) > numShells  /*ShellIndex*sqrt(2)>numShells+2*3+1*/ ){
              inputMSA[ShellIndex] = EPSILON_ZERO;
          }
    }
 }

//soft cutoff
//void msaCutoff (std::vector<double> & inputMSA, const unsigned long int nx, const unsigned long int ny, const unsigned long int nz){
void msaSoftCutoff (std::vector<double> & originalMsa, std::vector<double> & inputMSA, const unsigned long int nx, const unsigned long int ny, const unsigned long int nz){
    const double sizeMAV = __sizeMAV__;
    unsigned long int numShells = inputMSA.size();
    //std::vector<double> msaCutoffVector=inputMSA;
    double minDiff1 = 1;
    for (unsigned long int ShellIndex = 0; ShellIndex < numShells; ShellIndex++){
        if (inputMSA[ShellIndex] <= EPSILON_ZERO || ShellIndex > numShells / 3.0f  ){
            if (inputMSA[ShellIndex] <= EPSILON_ZERO){
                    inputMSA[ShellIndex] = EPSILON_ZERO;
            }else if (inputMSA[ShellIndex]<= EPSILON_ZERO || originalMsa[ShellIndex]<= EPSILON_ZERO ){
                inputMSA[ShellIndex] = EPSILON_ZERO;
            }else if ( ShellIndex>1 ){
                if (inputMSA[ShellIndex-1]<= EPSILON_ZERO || originalMsa[ShellIndex-1]<= EPSILON_ZERO ){
                    inputMSA[ShellIndex] = EPSILON_ZERO;
                }else{
                    double diff1=inputMSA[ShellIndex]/originalMsa[ShellIndex];
                    double diff0=inputMSA[ShellIndex-1]/originalMsa[ShellIndex-1];
                    if (diff1<diff0 && diff1<minDiff1){
                        minDiff1=diff1;
                        inputMSA[ShellIndex]=inputMSA[ShellIndex];
                    }else{
                        inputMSA[ShellIndex]=originalMsa[ShellIndex]*minDiff1;
                    }
                }
            }
        }else{
            minDiff1=inputMSA[ShellIndex]/originalMsa[ShellIndex];
        }
    }
 }

// **********************
// **********************
//
// compute resolution at a certain threshold
//
// **********************
// **********************
double computeResolution(std::vector<double> & vectorResolution, std::vector<double> & vectorFSC, double threshold){

    //ideal: check where the first two values below the threshold + interpolation
    //std::cerr<<"************************\n"<<threshold<<" -->   ";
    //for (int ii=0;ii<vectorResolution.size();ii++){
    //    std::cerr<<" "<<vectorResolution[ii] <<"("<<vectorFSC[ii]<<") ";
    //}
    //std::cerr<<"\n      [";
    
    double X1;
    double X2;
    double Y1;
    double Y2;
    double Y3;
    if ( vectorResolution.size() > 0){
        X1=vectorResolution[0];
        Y1=vectorFSC[0];
    }
    
    bool secondConsecutive = false;
    for (int ii=1;ii<vectorResolution.size();ii++){
        //std::cerr<<" "<<vectorResolution[ii] <<"("<<vectorFSC[ii]<<") ";
        if (vectorFSC[ii]<threshold){
            X2=vectorResolution[ii];
            Y2=vectorFSC[ii];
            double interpolatedX=X1+(threshold-Y1)*(X2-X1)/(Y2-Y1);
            //std::cerr<<"]\n  Y1="<<Y1<<"   Y2="<<Y2<<"   interpolatedX="<<1.0/interpolatedX<<"\n";
            //return 1.0/vectorResolution[ii];
            return 1.0/X2;
                    
            /*
            if (secondConsecutive){
                if(vectorResolution[ii]>0){
                    double interpolatedX=X1+(threshold-Y1)*(X2-X1)/(Y2-Y1);
                    std::cerr<<"]\n  Y1="<<Y1<<"   Y2="<<Y2<<"   interpolatedX="<<1.0/interpolatedX<<"\n";
                    //return 1.0/vectorResolution[ii];
                    return 1.0/X2;
                }
            }else{
                if (vectorResolution[ii]>0){
                    X2=vectorResolution[ii];
                    Y2=vectorFSC[ii];
                    secondConsecutive = true;
                }
            }
            */
        }else{
            secondConsecutive = false;
            X1=vectorResolution[ii];
            Y1=vectorFSC[ii];
        }

    }
    
//    std::cerr<<" //]\n";
    if (X1>0){
        return 1.0/X1;
    } else {
        return 0;
    }
}

// **********************************
// **********************************
//
/**         compute the Mean Shell Amplitude (MSA) of the input map.
 *
 * @return Returns a vector of the computed Mean Shell Amplitude (MSA) of the input map.
 */
//
// **********************************
// ***********************************
template<typename T>
std::vector<double>  computeMsa (const T * map, const T * maskVect,  long int nx,  long int ny, long int nz){

//writeVectorMrcImage("map1.mrc", map1,  nx,  ny,  1);
//writeVectorMrcImage("map2.mrc", map2,  nx,  ny,  1);

		const long int nxy = nx * ny;
		const long int nxyz = nx * ny * nz;
		const double X0 = nx/2.0;
		const double Y0 = ny/2.0;
		double Z0 = nz/2.0;
		if (nz<=1){nz=1;Z0=0;}
    //T * shellIndexImage = new T[nxyz];


    const int sizeMAV =  __sizeMAV__; //floor(  pow( floor( floor(pow(nz/2.0-1, 2))+ floor(pow(ny/2.0-1, 2))+ floor(pow(nx/2.0-1, 2)) ), 0.5 ) / SHELL_WIDTH ); //OK
    //double increment=1.0f/(double)sizeMAV;
    //std::cerr<<"increment="<<increment<<"\n";

    std::vector<double> vectorResolutionMSA;
    std::vector<double> vectorShellSize;
    std::vector<double> vectorMSA;

    for (long int ShellIndex = 0; ShellIndex < sizeMAV; ShellIndex++){
      vectorShellSize.push_back(0);
      vectorMSA.push_back( EPSILON_ZERO );
      vectorResolutionMSA.push_back( ShellIndex );
    }

    //
    double meanMap = 0.0;
    for  (long int ijk = 0; ijk < nxyz; ijk++){
      meanMap+=map[ijk];
    }
    meanMap/=(double)nxyz;

  //compute fourier transform of map
   fftw_complex *in, *F;
   in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
   F = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
   for (long int ijk = 0; ijk < nxyz; ijk++){
      in[ijk][0] = (map[nxyz-ijk-1] -(double)meanMap);
      in[ijk][1] = 0.0; //no imaginary part
   }
   fftw_plan planForwardMap;
   if (nz>1){
       planForwardMap = fftw_plan_dft_3d(nz, ny, nx, in, F, FFTW_FORWARD, FFTW_ESTIMATE);
   }else{
       planForwardMap = fftw_plan_dft_2d(ny, nx, in, F, FFTW_FORWARD, FFTW_ESTIMATE);
   }
   fftw_execute(planForwardMap);
   fftw_destroy_plan(planForwardMap);
   fftw_free(in);

    for ( long int zyxIndex=0; zyxIndex < nxyz; zyxIndex++){
      F[ zyxIndex ][0] /= (double)nxyz;
      F[ zyxIndex ][1] /= (double)nxyz;
    }

    for ( long int kk = 0, zyxIndex=nxyz-1; kk < nz; kk++){
        //double ZSquared=pow((kk-floor(Z0))/Z0, 2);
        double Z1Squared=0;
        if(Z0>0) Z1Squared=kk<Z0?pow(kk/(double)Z0,2):pow((nz-kk-1)/(double)Z0,2); //pow((kk-floor(Z0))/Z0, 2);
        for ( long int jj = 0; jj < ny; jj++){
             //double YZSquared=pow((jj-floor(Y0))/Y0, 2)+ZSquared;
             double YZ1Squared=Z1Squared;
             if(Y0>0) YZ1Squared=jj<Y0?pow(jj/(double)Y0,2)+Z1Squared:pow((ny-jj-1)/(double)Y0,2)+Z1Squared;
             for ( long int ii = 0; ii < nx; ii++, zyxIndex--){
                 bool keepMaskedValue = true;
                 if(maskVect!=NULL){
                   if (maskVect[zyxIndex] <= MASK_THRESHOLD_VALUE ) {
                       keepMaskedValue=false;
                   }
                 }                        
                 if ( keepMaskedValue ) {
                      //double XYZSquared=(pow((ii-floor(X0))/X0, 2)+YZSquared);
                      double XYZ1Squared=YZ1Squared;
                      if (X0>0) XYZ1Squared=ii<X0?pow(ii/(double)X0,2)+YZ1Squared:pow((nx-ii-1)/(double)X0,2)+YZ1Squared;
                      double XYZ=0;
                      if (SHELL_WIDTH>0) XYZ=pow(XYZ1Squared, 0.5)/(double)SHELL_WIDTH;
                      long int ShellIndex = floor( XYZ*sizeMAV );
                      //testImage[ zyxIndex ]=XYZSquared;
                      //testImage1[ zyxIndex ]=XYZ1Squared;
                      //shellIndexImage[zyxIndex]=ShellIndex;
                      double mapAmplitudeSquared = F[ zyxIndex ][0]*F[ zyxIndex ][0]+F[ zyxIndex ][1]*F[ zyxIndex ][1];
                      double mapAmplitude = pow(mapAmplitudeSquared, 0.5);
                      if(ShellIndex < sizeMAV /*&&XYZ<0.5774*/){// 1/sqrt(3)=0.5774
                           vectorShellSize[ShellIndex]++;
                           vectorMSA[ShellIndex]+=mapAmplitude;
                      }
                 }
             }
         }
    }
    // writeVectorMrcImage("shellIndexImage.mrc", shellIndexImage,  nx,  ny,  1);
    fftw_free(F);

    for (long int ShellIndex = 0; ShellIndex < sizeMAV; ShellIndex++){
      double sizeShell = (double)vectorShellSize[ShellIndex];
        if (sizeShell>0.0 ){
          double value=vectorMSA[ShellIndex] / (double)sizeShell;
          if (value==value) vectorMSA[ShellIndex] = value;
          else vectorMSA[ShellIndex] = 0;
        }else{
          vectorMSA[ShellIndex]=0;
        }
    }	

  // delete [] shellIndexImage;
  //plot2d (vectorMSA, vectorMSA, vectorResolutionMSA);
  return vectorMSA;
}




// ****************************
//computeMsa2d

template<typename T>
std::vector<double>  computeMsa2D (const T * map,  const unsigned long int nx,  const unsigned long int ny){

//writeVectorMrcImage("map1.mrc", map1,  nx,  ny,  1);
//writeVectorMrcImage("map2.mrc", map2,  nx,  ny,  1);

        const unsigned long int nxy = nx * ny;
        const unsigned long int nxyz = nx * ny;
        const float X0 = nx/2.0;
        const float Y0 = ny/2.0;

    //T * shellIndexImage = new T [nxyz];

    const unsigned int sizeMAV =  __sizeMAV__2D__;

    std::vector<double> vectorResolutionMSA;
    std::vector<double> vectorShellSize;
    std::vector<double> vectorMSA;

        for (unsigned int ShellIndex = 0; ShellIndex < sizeMAV; ShellIndex++){
      vectorShellSize.push_back(0);
      vectorMSA.push_back(EPSILON_ZERO);
      vectorResolutionMSA.push_back(ShellIndex );
        }

    //
    double meanMap = 0;
    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
      meanMap+=map[ijk];
    }
    meanMap/=nxyz;

  //compute fourier transform of map
    fftw_complex *inF, *F;
        inF = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
        F = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
                    inF[ijk][0] = (map[nxyz-ijk-1]-meanMap);
                    inF[ijk][1] = 0.0; //no imaginary part
    }
        fftw_plan planForwardMap;

    planForwardMap = fftw_plan_dft_2d(ny, nx, inF, F, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(planForwardMap);
        fftw_destroy_plan(planForwardMap);
        fftw_free(inF);

    for (unsigned long int zyxIndex=0; zyxIndex < nxyz; zyxIndex++){
      F[ zyxIndex ][0] /= nxyz;
      F[ zyxIndex ][1] /= nxyz;
    }
        for (unsigned long int kk = 0, zyxIndex=nxyz-1; kk < 1; kk++){
                //double ZSquared=pow((kk-floor(Z0))/Z0, 2);
                double Z1Squared=0; //pow((kk-floor(Z0))/Z0, 2);
                for (unsigned long int jj = 0; jj < ny; jj++){
                    //double YZSquared=pow((jj-floor(Y0))/Y0, 2)+ZSquared;
                    double YZ1Squared=jj<Y0?pow(jj/Y0,2)+Z1Squared:pow((ny-jj-1)/Y0,2)+Z1Squared;
                    for (unsigned long int ii = 0; ii < nx; ii++, zyxIndex--){
                                //double XYZSquared=(pow((ii-floor(X0))/X0, 2)+YZSquared);
                                double XYZ1Squared=ii<X0?pow(ii/X0,2)+YZ1Squared:pow((nx-ii-1)/X0,2)+YZ1Squared;
                                double XYZ=pow(XYZ1Squared, 0.5)/SHELL_WIDTH;
                                unsigned int ShellIndex = floor( XYZ*sizeMAV );
                                //testImage[ zyxIndex ]=XYZSquared;
                                //testImage1[ zyxIndex ]=XYZ1Squared;
                                //shellIndexImage[zyxIndex]=ShellIndex;
                                double mapAmplitudeSquared = F[ zyxIndex ][0]*F[ zyxIndex ][0]+F[ zyxIndex ][1]*F[ zyxIndex ][1];
                                double mapAmplitude = pow(mapAmplitudeSquared, 0.5);
                        if(ShellIndex < vectorShellSize.size() /*&&XYZ<0.5774*/){// 1/sqrt(3)=0.5774
                                        vectorShellSize[ShellIndex]++;
                                        vectorMSA[ShellIndex]+=mapAmplitude;
                                    }

                    }
                }
            }
        fftw_free(F);

    
    //writeVectorMrcImage("shellIndexImage.mrc", shellIndexImage,  nx,  ny,  1);
    for (unsigned int ShellIndex = 0; ShellIndex < vectorShellSize.size(); ShellIndex++){
      double sizeShell = (double)vectorShellSize[ShellIndex];
        if (sizeShell>0 ){
          vectorMSA[ShellIndex] /= sizeShell;
        }else{
          vectorMSA[ShellIndex]=0;
        }
        }

//  delete [] shellIndexImage;
  //plot2d (vectorMSA, vectorMSA, vectorResolutionMSA);
  return vectorMSA;
}

template<typename T>
void applylMsa2D ( const T * mapIn, T * mapOut, std::vector<double> replacingMSA, const unsigned long int nx,  const unsigned long int ny ){

        const unsigned long int nxy = nx * ny;
        const unsigned long int nxyz = nxy;
        const float X0 = nx/2.0;
        const float Y0 = ny/2.0;


    unsigned long int sizeMAV =  __sizeMAV__2D__;
    if (sizeMAV != replacingMSA.size()){
      std::cerr<<"ERROR in applyMsa: unexpected number of shells    Details:   ";
      std::cerr<<"sizeMAV="<<sizeMAV<<"   replacingMSA="<< replacingMSA.size()<<"\n";
      return;
    }
    //std::cerr<<"sizeMAV="<<sizeMAV<<"\n";
    if (sizeMAV < 1){
      std::cerr<<"ERROR: empty msa list.\n";
      return;
    }

    //compute mean
    double meanIn = 0;
    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
      meanIn+=mapIn[ijk];
    }
    meanIn/=nxyz;
    //std::cerr<<"mean="<<meanIn<<"\n";

    std::vector<double> vectorShellSize;
    std::vector<double> originalMSA;
    std::vector<double> normalizationMSA;

    for (unsigned int ShellIndex = 0; ShellIndex < replacingMSA.size(); ShellIndex++){
      vectorShellSize.push_back(0);
      originalMSA.push_back(0);
      normalizationMSA.push_back(replacingMSA[ShellIndex]);
     }
    fftw_complex *I, *F;
        I = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
        F = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
                    I[ijk][0] = (mapIn[nxyz-ijk-1]-meanIn);
                    I[ijk][1] = 0.0; //no imaginary part
    }
    fftw_plan planForwardMap;
    planForwardMap = fftw_plan_dft_2d(ny, nx, I, F, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(planForwardMap);
    fftw_destroy_plan(planForwardMap);


    for (unsigned long int zyxIndex=0; zyxIndex < nxyz; zyxIndex++){
      F[ zyxIndex ][0] /= nxyz;
      F[ zyxIndex ][1] /= nxyz;
    }



    for (unsigned long int kk = 0, zyxIndex=nxyz-1; kk < 1; kk++){
        double ZSquared=0;
            for (unsigned long int jj = 0; jj < ny; jj++){
          double YZSquared=jj<Y0?pow(jj/Y0,2)+ZSquared:pow((ny-jj-1)/Y0,2)+ZSquared;
          for (unsigned long int ii = 0; ii < nx; ii++, zyxIndex--){

                double XYZSquared=ii<X0?pow(ii/X0,2)+YZSquared:pow((nx-ii-1)/X0,2)+YZSquared;
                double XYZ = pow(XYZSquared, 0.5)/SHELL_WIDTH;
                unsigned int ShellIndex = floor( XYZ*sizeMAV );
                double originalAmplitude = pow(F[ zyxIndex ][0]*F[ zyxIndex ][0]+F[ zyxIndex ][1]*F[ zyxIndex ][1], 0.5);
              if(ShellIndex < originalMSA.size() ){
                    originalMSA[ShellIndex]+=originalAmplitude;
                    vectorShellSize[ShellIndex]++;
                }

          }
        }
    }



    for (unsigned int ShellIndex = 0; ShellIndex < vectorShellSize.size(); ShellIndex++){
    unsigned long int sizeShell = vectorShellSize[ShellIndex];
    if ( sizeShell > 0 && originalMSA[ShellIndex] > 0.0f ){
      originalMSA[ShellIndex] /= sizeShell;
      normalizationMSA[ShellIndex] = normalizationMSA[ShellIndex]/originalMSA[ShellIndex];
    }else{
      normalizationMSA[ShellIndex]=0.0f;
    }

  }
    unsigned long int counterOut=0;
    unsigned long int counterTotal=0;
    for (unsigned long int kk = 0, zyxIndex=nxyz-1; kk < 1; kk++){
            double ZSquared=0; //pow((kk-floor(Z0))/Z0, 2);
        for (unsigned long int jj = 0; jj < ny; jj++){
             double YZSquared=jj<Y0?pow(jj/Y0,2)+ZSquared:pow((ny-jj-1)/Y0,2)+ZSquared;
            for (unsigned long int ii = 0; ii < nx; ii++, zyxIndex--){
                    double XYZSquared=ii<X0?pow(ii/X0,2)+YZSquared:pow((nx-ii-1)/X0,2)+YZSquared;
                    double XYZ=pow(XYZSquared, 0.5)/SHELL_WIDTH;
                    unsigned int ShellIndex = floor( XYZ*sizeMAV );
//                    counterTotal++;
                if (ShellIndex>=vectorShellSize.size() ){
                            I[zyxIndex][0]=0;
                            I[zyxIndex][1]=0;
                    }else{
                            I[zyxIndex][0]=F[ zyxIndex ][0]*normalizationMSA[ShellIndex];
                            I[zyxIndex][1]=F[ zyxIndex ][1]*normalizationMSA[ShellIndex];
                    }
          }
        }
    }
    
  fftw_plan planForwardOriginalImage;
  planForwardOriginalImage = fftw_plan_dft_2d(ny, nx, I, F, FFTW_BACKWARD, FFTW_ESTIMATE);


  fftw_execute(planForwardOriginalImage);
  fftw_destroy_plan(planForwardOriginalImage);

    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
          mapOut[ijk] = meanIn+(F[nxyz-ijk-1][0]);
    }
    fftw_free(I);
    fftw_free(F);
    
}




// **********************************
// **********************************
//
//         countShellSize
//
// **********************************
// ***********************************
/*
template<typename T>
std::vector<double>  countShellSize (const T * map,  const unsigned long int nx,  const unsigned long int ny, const  unsigned long int nz){

//writeVectorMrcImage("map1.mrc", map1,  nx,  ny,  1);
//writeVectorMrcImage("map2.mrc", map2,  nx,  ny,  1);

		const unsigned long int nxy = nx * ny;
		const unsigned long int nxyz = nx * ny * nz;
		const float X0 = nx/2.0;
		const float Y0 = ny/2.0;
		const float Z0 = nz/2.0;
    //T * shellIndexImage = new T[nxyz];
		double nD = 0.0f;
		if (nx > 1 ) nD++;
		if (ny > 1 ) nD++;
		if (nz > 1 ) nD++;

    double minDimensionSize = nx;
    if (minDimensionSize > ny )
       minDimensionSize = ny;
    if (minDimensionSize > nz && nz > 4)
       minDimensionSize = nz;

    double maxDimensionSize = nx;
    if (maxDimensionSize < ny )
       maxDimensionSize = ny;
    if (maxDimensionSize < nz && nz > 4)
       maxDimensionSize = nz;


    const unsigned int sizeMAV =  __sizeMAV__; //floor(  pow( floor( floor(pow(nz/2.0-1, 2))+ floor(pow(ny/2.0-1, 2))+ floor(pow(nx/2.0-1, 2)) ), 0.5 ) / SHELL_WIDTH ); //OK
    //double increment=1.0f/(double)sizeMAV;
    //std::cerr<<"increment="<<increment<<"\n";

    std::vector<double> vectorResolutionMSA;
    std::vector<double> vectorShellSize;
    std::vector<double> vectorMSA;

		for (unsigned int ShellIndex = 0; ShellIndex < sizeMAV; ShellIndex++){
      vectorShellSize.push_back(0);
      vectorMSA.push_back(0);
      vectorResolutionMSA.push_back(ShellIndex);
		}

    //
    double meanMap = 0;
    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
      meanMap+=map[ijk];
    }
    meanMap/=nxyz;

  //compute fourier transform of map
	fftw_complex *in, *F;
		in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
		F = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
					in[ijk][0] = (map[nxyz-ijk-1]-meanMap);
					in[ijk][1] = 0.0; //no imaginary part
    }
		fftw_plan planForwardMap;
    if (nz>1){
		  planForwardMap = fftw_plan_dft_3d(nz, ny, nx, in, F, FFTW_FORWARD, FFTW_ESTIMATE);
    }else{
		  planForwardMap = fftw_plan_dft_2d(ny, nx, in, F, FFTW_FORWARD, FFTW_ESTIMATE);
    }
    fftw_execute(planForwardMap);
		fftw_destroy_plan(planForwardMap);
		fftw_free(in);

    for (unsigned long int zyxIndex=0; zyxIndex < nxyz; zyxIndex++){
      F[ zyxIndex ][0] /= nxyz;
      F[ zyxIndex ][1] /= nxyz;
    }

		for (unsigned long int kk = 0, zyxIndex=nxyz-1; kk < nz; kk++){
				double ZSquared=pow((kk-floor(Z0))/Z0, 2);
        double Z1Squared=kk<Z0?pow(kk/Z0,2):pow((nz-kk-1)/Z0,2); //pow((kk-floor(Z0))/Z0, 2);
				for (unsigned long int jj = 0; jj < ny; jj++){
					double YZSquared=pow((jj-floor(Y0))/Y0, 2)+ZSquared;
          double YZ1Squared=jj<Y0?pow(jj/Y0,2)+Z1Squared:pow((ny-jj-1)/Y0,2)+Z1Squared;
					for (unsigned long int ii = 0; ii < nx; ii++, zyxIndex--){
							double XYZSquared=(pow((ii-floor(X0))/X0, 2)+YZSquared);
              double XYZ1Squared=ii<X0?pow(ii/X0,2)+YZ1Squared:pow((nx-ii-1)/X0,2)+YZ1Squared;
              double XYZ=pow(XYZ1Squared, 0.5)/SHELL_WIDTH;
							unsigned int ShellIndex = floor( XYZ*sizeMAV );
              //testImage[ zyxIndex ]=XYZSquared;
              //testImage1[ zyxIndex ]=XYZ1Squared;
              //shellIndexImage[zyxIndex]=ShellIndex;
              double mapAmplitudeSquared = F[ zyxIndex ][0]*F[ zyxIndex ][0]+F[ zyxIndex ][1]*F[ zyxIndex ][1];
							double mapAmplitude = pow(mapAmplitudeSquared, 0.5);
							if(ShellIndex < sizeMAV ){// 1/sqrt(3)=0.5774
								vectorShellSize[ShellIndex]++;
								vectorMSA[ShellIndex]+=mapAmplitude;
							}
					}
				}
			}
    // writeVectorMrcImage("shellIndexImage.mrc", shellIndexImage,  nx,  ny,  1);
		fftw_free(F);

		for (unsigned int ShellIndex = 0; ShellIndex < sizeMAV; ShellIndex++){
      double sizeShell = (double)vectorShellSize[ShellIndex];
        if (sizeShell>0 ){
          vectorMSA[ShellIndex] /= sizeShell;
        }else{
          vectorMSA[ShellIndex]=EPSILON_ZERO;
        }
		}

  // delete [] shellIndexImage;
  //plot2d (vectorMSA, vectorMSA, vectorResolutionMSA);
  return vectorShellSize;
}

*/

// **********************************
// **********************************
//
//         ApplyMSA
//
// **********************************
// ***********************************
template<typename T>
void applylMsa ( const T * mapIn, T * mapOut, const T * maskVect, const std::vector<double> & replacingMSA, const unsigned long int nx,  const unsigned long int ny, unsigned long int nz){

		const unsigned long int nxy = nx * ny;
		const unsigned long int nxyz = nx * ny * nz;
		const double X0 = nx/2.0;
		const double Y0 = ny/2.0;
		double Z0 = nz/2.0;
		double nD = 0.0f;
		if (nx > 1 ) nD++;
		if (ny > 1 ) nD++;
		if (nz > 1 ) nD++;
		if (nz<=1) {Z0=0;nz=1;}

//    float * testImage = new float [nxyz];


    unsigned long int sizeMAV =  __sizeMAV__;
    if (sizeMAV > replacingMSA.size()){
      std::cerr<<"ERROR in applyMsa: unexpected number of shells    Details:   ";
      std::cerr<<"sizeMAV="<<sizeMAV<<"   replacingMSA="<< replacingMSA.size()<<"\n";
      return;
    }
    //std::cerr<<"sizeMAV="<<sizeMAV<<"\n";
    if (sizeMAV < 1){
      std::cerr<<"ERROR: empty msa list.\n";
      return;
    }

    //compute mean
    double meanIn = 0;
    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
      meanIn+=mapIn[ijk];
    }
    meanIn/=nxyz;
    //std::cerr<<"mean="<<meanIn<<"\n";

    std::vector<double> vectorShellSize;
    std::vector<double> originalMSA;
    std::vector<double> normalizationMSA;

		for (unsigned int ShellIndex = 0; ShellIndex < sizeMAV; ShellIndex++){
            vectorShellSize.push_back(0.00001);
            originalMSA.push_back(0.00001);
            normalizationMSA.push_back(replacingMSA[ShellIndex]);
		}


	fftw_complex *I, *F;
		I = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
		F = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
					I[ijk][0] = (mapIn[nxyz-ijk-1]-meanIn);
					I[ijk][1] = 0.0; //no imaginary part
    }
		fftw_plan planForwardMap;
    if (nz>1){
		  planForwardMap = fftw_plan_dft_3d(nz, ny, nx, I, F, FFTW_FORWARD, FFTW_ESTIMATE);
    }else{
      planForwardMap = fftw_plan_dft_2d(ny, nx, I, F, FFTW_FORWARD, FFTW_ESTIMATE);
    }

		fftw_execute(planForwardMap);
		fftw_destroy_plan(planForwardMap);


    for (unsigned long int zyxIndex=0; zyxIndex < nxyz; zyxIndex++){
      F[ zyxIndex ][0] /= nxyz;
      F[ zyxIndex ][1] /= nxyz;
    }



	for (unsigned long int kk = 0, zyxIndex=nxyz-1; kk < nz; kk++){
	        double ZSquared=0;
	        if (Z0>0)ZSquared=kk<Z0?pow(kk/Z0,2):pow((nz-kk-1)/Z0,2); //pow((kk-floor(Z0))/Z0, 2);
			for (unsigned long int jj = 0; jj < ny; jj++){
          double YZSquared=ZSquared;
          if (Y0>0) YZSquared=jj<Y0?pow(jj/Y0,2)+ZSquared:pow((ny-jj-1)/Y0,2)+ZSquared;
          for (unsigned long int ii = 0; ii < nx; ii++, zyxIndex--){
              bool keepMaskedValue = true;
              if(maskVect!=NULL){
                  if (maskVect[zyxIndex] <= MASK_THRESHOLD_VALUE ) {
                      keepMaskedValue=false;
                  }
              }
              if ( keepMaskedValue ) {
                double XYZSquared=YZSquared;
                if (X0>0) XYZSquared=ii<X0?pow(ii/X0,2)+YZSquared:pow((nx-ii-1)/X0,2)+YZSquared;
                double XYZ = pow(XYZSquared, 0.5)/SHELL_WIDTH;
                unsigned int ShellIndex = floor( XYZ*sizeMAV );
                double originalAmplitude = pow(F[ zyxIndex ][0]*F[ zyxIndex ][0]+F[ zyxIndex ][1]*F[ zyxIndex ][1], 0.5);
                if(ShellIndex < sizeMAV ){
                    originalMSA[ShellIndex]+=originalAmplitude;
                    vectorShellSize[ShellIndex]++;
                }
            }
          }
        }
    }



  for (unsigned int ShellIndex = 0; ShellIndex < sizeMAV; ShellIndex++){
    unsigned long int sizeShell = vectorShellSize[ShellIndex];
    if ( sizeShell > 0 && originalMSA[ShellIndex] > 0.00000001f ){
      originalMSA[ShellIndex] /= sizeShell;
      normalizationMSA[ShellIndex] = normalizationMSA[ShellIndex]/originalMSA[ShellIndex];
    }else{
      normalizationMSA[ShellIndex]=EPSILON_ZERO;
    }
    if (!(normalizationMSA[ShellIndex]==normalizationMSA[ShellIndex])) normalizationMSA[ShellIndex]=EPSILON_ZERO;
  }

	unsigned long int counterOut=0;
	unsigned long int counterTotal=0;
	for (unsigned long int kk = 0, zyxIndex=nxyz-1; kk < nz; kk++){
        	double ZSquared=0;
        	if (Z0>0) ZSquared=kk<Z0?pow(kk/Z0,2):pow((nz-kk-1)/Z0,2); //pow((kk-floor(Z0))/Z0, 2);
		for (unsigned long int jj = 0; jj < ny; jj++){
 			double YZSquared=ZSquared;
 			if (Y0>0) YZSquared=jj<Y0?pow(jj/Y0,2)+ZSquared:pow((ny-jj-1)/Y0,2)+ZSquared;
			for (unsigned long int ii = 0; ii < nx; ii++, zyxIndex--){
                bool skipMaskedValue = false;
                if(maskVect!=NULL){
                  if (maskVect[zyxIndex] <= MASK_THRESHOLD_VALUE ) {
                      skipMaskedValue=true;
                  }
                }                
                if ( skipMaskedValue ){
                    I[zyxIndex][0]=F[ zyxIndex ][0];
                    I[zyxIndex][1]=F[ zyxIndex ][1];
                }else{
                    double XYZSquared=YZSquared;
                    if(X0>0)XYZSquared=ii<X0?pow(ii/X0,2)+YZSquared:pow((nx-ii-1)/X0,2)+YZSquared;
                    unsigned int ShellIndex = floor( pow(XYZSquared, 0.5)*sizeMAV/SHELL_WIDTH );
//                    counterTotal++;
                    if (ShellIndex>=sizeMAV){
//                        counterOut++;
    //					std::cerr<< "(x,y,k)=(" << ii/2<<","<<jj/2<<","<<kk/2<<")   SHELL="<<ShellIndex<<"  of "<< sizeMAV << "\n";
                        ShellIndex = sizeMAV - 1;
                        
                            I[zyxIndex][0]=0;
                            I[zyxIndex][1]=0;
                            

                    }else{
                            I[zyxIndex][0]=F[ zyxIndex ][0]*normalizationMSA[ShellIndex];
                            I[zyxIndex][1]=F[ zyxIndex ][1]*normalizationMSA[ShellIndex];
                    }
                }
                
				//testImage[zyxIndex]=ShellIndex;


          }
        }
    }
    
					//std::cerr<<"counterTotal = " << counterTotal << "   counterOut= "<< counterOut << "   RATIO="<< counterTotal/counterOut <<"\n";
  fftw_plan planForwardOriginalImage;
 if (nz>1){
    planForwardOriginalImage = fftw_plan_dft_3d(nz, ny, nx, I, F, FFTW_BACKWARD, FFTW_ESTIMATE);
  }else{
    planForwardOriginalImage = fftw_plan_dft_2d(ny, nx, I, F, FFTW_BACKWARD, FFTW_ESTIMATE);
  }

  fftw_execute(planForwardOriginalImage);
  fftw_destroy_plan(planForwardOriginalImage);

    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
          mapOut[ijk] = meanIn+(F[nxyz-ijk-1][0]);
    }
		fftw_free(I);
		fftw_free(F);

//	writeVectorMrcImage("index.mrc", testImage,  nx,  ny,  nz);
//	delete[]testImage;

    
}





// **********************************
// **********************************
//
/** Compute optimal Mean Shell Amplitude (MSA) and Fourier Shell Correlation (FSC)
 * of two half maps, using the Figure of Merit computation
 */
//
// **********************************
// ***********************************
template<typename T>
std::vector<double>  computeHalfMapsMsaFSC (const T * map1, //!< first half map
const T * map2, //!< second half map
const T * maskVect, //!< mask (only positive values of the mask are taken into account)
std::vector<double> & vectorFSC, //!< output vector of FSC values
const unsigned long int nx, //!< size of image (X)
const unsigned long int ny, //!< size of image (Y)
const  unsigned long int nz, //!< size of image (Z)
double apix=1 //!< Pixel size in Angstrom
){

//writeVectorMrcImage("map1.mrc", map1,  nx,  ny,  1);
//writeVectorMrcImage("map2.mrc", map2,  nx,  ny,  1);

		const unsigned long int nxy = nx * ny;
		const unsigned long int nxyz = nx * ny * nz;
		const float X0 = nx/2.0;
		const float Y0 = ny/2.0;
		const float Z0 = nz/2.0;

		double nD = 0.0f;
		if (nx > 1 ) nD++;
		if (ny > 1 ) nD++;
		if (nz > 1 ) nD++;

    double minDimensionSize = nx;
    if (minDimensionSize > ny )
       minDimensionSize = ny;
    if (minDimensionSize > nz && nz > 4)
       minDimensionSize = nz;

    double maxDimensionSize = nx;
    if (maxDimensionSize < ny )
       maxDimensionSize = ny;
    if (maxDimensionSize < nz && nz > 4)
       maxDimensionSize = nz;



    //increment = 2/(minDimensionSize+1);
		//unsigned int sizeMAV = floor(1.0/increment - 2);
		//const double epsilon=0.00001;
    //float * testImage = new float [nxyz];

    const unsigned int sizeMAV =  __sizeMAV__;//floor(  pow( floor( floor(pow(nz/2.0-1, 2))+ floor(pow(ny/2.0-1, 2))+ floor(pow(nx/2.0-1, 2)) ), 0.5 ) / SHELL_WIDTH ); //OK
    //double increment=1.0f/(double)sizeMAV;

    //std::cerr<<"increment="<<increment<<"\n";

		std::vector<double> vectorMap1MSA;
		std::vector<double> vectorMap2MSA;
    		std::vector<double> vectorMap12MSA;
		std::vector<double> vectorMap1SumSquaredAmplitude;
		std::vector<double> vectorMap2SumSquaredAmplitude;
    //std::vector<double> vectorFSC;
    vectorFSC.clear();
    std::vector<double> vectorResolutionMSA;
    std::vector<double> vectorShellSize;
    std::vector<double> vectorShellExpectedSize;
    std::vector<double> resultMSA;


		for (unsigned int ShellIndex = 0; ShellIndex < sizeMAV; ShellIndex++){
      vectorMap1MSA.push_back(0);
      vectorMap2MSA.push_back(0);
      vectorMap12MSA.push_back(0);
      vectorMap1SumSquaredAmplitude.push_back(0);
      vectorMap2SumSquaredAmplitude.push_back(0);
      vectorFSC.push_back(0);
      vectorShellSize.push_back(0);
      resultMSA.push_back(0);
      vectorResolutionMSA.push_back(ShellIndex );
		}

    //
    double meanMap1 = 0,  meanMap2 = 0;
    for (unsigned long int ijk = 0; ijk < nxyz; ijk++){
      meanMap1+=map1[ijk];
      meanMap2+=map2[ijk];
    }
    meanMap1/=nxyz;
    meanMap2/=nxyz;


  //Plan:
  //compute fourier transform of map1 and map2
  //compute mean shell amplitude of map1 and map2
  //compute fourier shell correlation
  //compute sqrt(fsc1/(1+fsc))

  //compute fourier transform of map1 and map2
	fftw_complex *in1, *in2, *in12, *F1, *F2, *F12;
		in1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
    in2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
		F1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
    F2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nxyz);
		for (unsigned long int kk = 0, zyxIndex=0; kk < nz; kk++){
			for (unsigned long int jj = 0; jj < ny; jj++){
				for (unsigned long int ii = 0; ii < nx; ii++, zyxIndex++){
					in1[zyxIndex][0] = (map1[nxyz-zyxIndex]-meanMap1);
					in1[zyxIndex][1] = 0.0; //no imaginary part
					in2[zyxIndex][0] = (map2[nxyz-zyxIndex]-meanMap2);
					in2[zyxIndex][1] = 0.0; //no imaginary part
			 }
			}
		}
		fftw_plan planForwardMap1;
    fftw_plan planForwardMap2;
    if (nz>1){
		  planForwardMap1 = fftw_plan_dft_3d(nz, ny, nx, in1, F1, FFTW_FORWARD, FFTW_ESTIMATE);
      planForwardMap2 = fftw_plan_dft_3d(nz, ny, nx, in2, F2, FFTW_FORWARD, FFTW_ESTIMATE);
    }else{
		  planForwardMap1 = fftw_plan_dft_2d(ny, nx, in1, F1, FFTW_FORWARD, FFTW_ESTIMATE);
      planForwardMap2 = fftw_plan_dft_2d(ny, nx, in2, F2, FFTW_FORWARD, FFTW_ESTIMATE);
    }
    fftw_execute(planForwardMap1);
		fftw_destroy_plan(planForwardMap1);
 		fftw_execute(planForwardMap2);
		fftw_destroy_plan(planForwardMap2);

		fftw_free(in1);
		fftw_free(in2);

    for (unsigned long int zyxIndex=0; zyxIndex < nxyz; zyxIndex++){
      F1[ zyxIndex ][0] /= nxyz;
      F1[ zyxIndex ][1] /= nxyz;
      F2[ zyxIndex ][0] /= nxyz;
      F2[ zyxIndex ][1] /= nxyz;
    }
  //WorkingPixelType * testImage = new WorkingPixelType [nxyz];
  //WorkingPixelType * testImage1 = new WorkingPixelType [nxyz];
  //compute mean shell amplitude of map1 and map2
  //compute fourier shell correlation
  //compute sqrt(fsc1/(1+fsc))
		for (unsigned long int kk = 0, zyxIndex=nxyz-1; kk < nz; kk++){
            double ZSquared=pow((kk-floor(Z0))/Z0, 2);
            double Z1Squared=kk<Z0?pow(kk/Z0,2):pow((nz-kk-1)/Z0,2); //pow((kk-floor(Z0))/Z0, 2);
				for (unsigned long int jj = 0; jj < ny; jj++){
					double YZSquared=pow((jj-floor(Y0))/Y0, 2)+ZSquared;
                    double YZ1Squared=jj<Y0?pow(jj/Y0,2)+Z1Squared:pow((ny-jj-1)/Y0,2)+Z1Squared;
					for (unsigned long int ii = 0; ii < nx; ii++, zyxIndex--){
                        bool keepMaskedValue = true;
                        if(maskVect!=NULL){
                            if (maskVect[zyxIndex] <= MASK_THRESHOLD_VALUE ) {
                                keepMaskedValue=false;
                            }
                        }                        
                        if ( keepMaskedValue ) {
							double XYZSquared=(pow((ii-floor(X0))/X0, 2)+YZSquared);
                            double XYZ1Squared=ii<X0?pow(ii/X0,2)+YZ1Squared:pow((nx-ii-1)/X0,2)+YZ1Squared;
                            double XYZ=pow(XYZ1Squared, 0.5)/SHELL_WIDTH;
                            unsigned int ShellIndex = floor( XYZ*sizeMAV );
                            //testImage[ zyxIndex ]=ShellIndex;
                            //testImage1[ zyxIndex ]=XYZ1Squared;
                            //shellIndexImage[zyxIndex]=ShellIndex;
                            double map1AmplitudeSquared = F1[ zyxIndex ][0]*F1[ zyxIndex ][0]+F1[ zyxIndex ][1]*F1[ zyxIndex ][1];
                            double map2AmplitudeSquared = F2[ zyxIndex ][0]*F2[ zyxIndex ][0]+F2[ zyxIndex ][1]*F2[ zyxIndex ][1];
                            double map12AmplitudeSquared = pow((F1[ zyxIndex ][0]+F2[ zyxIndex ][0])/2,2)+pow((F1[ zyxIndex ][1]+F2[ zyxIndex ][1])/2,2);
                            double map1Amplitude = pow(map1AmplitudeSquared, 0.5);
                            double map2Amplitude = pow(map2AmplitudeSquared, 0.5);
                            double map12Amplitude = pow(map12AmplitudeSquared, 0.5);
                            //unsigned int ShellIndex2 = floor( XYZ*sizeMAV );
                            //testImage[ zyxIndex ]=XYZ1Squared;
                            //testImage1[ zyxIndex ]=XYZSquared;

                            if(ShellIndex < sizeMAV /*&&XYZ<0.5774*/){// 1/sqrt(3)=0.5774
								vectorShellSize[ShellIndex]++;
								vectorMap1MSA[ShellIndex]+=map1Amplitude;
                                vectorMap2MSA[ShellIndex]+=map2Amplitude;
                                vectorMap12MSA[ShellIndex]+=map12Amplitude;
								vectorMap1SumSquaredAmplitude[ShellIndex]+=map1AmplitudeSquared;
                                vectorMap2SumSquaredAmplitude[ShellIndex]+=map2AmplitudeSquared;
                                vectorFSC[ShellIndex]+=F1[ zyxIndex ][0]*F2[ zyxIndex ][0]+F1[ zyxIndex ][1]*F2[ zyxIndex ][1];
							}
                        }
					}
				}
			}
//    writeVectorMrcImage("index.mrc", testImage,  nx,  ny,  nz);
//    writeVectorMrcImage("testOld.mrc", testImage1,  nx,  ny,  nz);
//    writeVectorMrcImage("testImage.mrc", testImage,  nx,  ny,  nz);
		fftw_free(F1);
		fftw_free(F2);


    for (unsigned int ShellIndex = 0; ShellIndex < sizeMAV; ShellIndex++){
        vectorResolutionMSA[ShellIndex]=1000*1.0/(ShellIndex*apix); //TEST
        double sizeShell = (double)vectorShellSize[ShellIndex];
        
        if (sizeShell>0){

          double denomFSC = pow(vectorMap1SumSquaredAmplitude[ShellIndex]*vectorMap2SumSquaredAmplitude[ShellIndex], 0.5);
          //std::cerr << ShellIndex << "  ->  " << vectorFSC[ShellIndex] << " / " << denomFSC << " = " << (double)(vectorFSC[ShellIndex]/(double)denomFSC)  << "\n";
          if (denomFSC> 0.0f){
            vectorFSC[ShellIndex] /= denomFSC;
          }else{
              vectorFSC[ShellIndex]=0;
          }
          vectorMap1MSA[ShellIndex] /= sizeShell;
          vectorMap2MSA[ShellIndex] /= sizeShell;
          vectorMap12MSA[ShellIndex] /= sizeShell;
          //double normFoM = sqrt(2);
          double figureOfMerit;
          if (vectorFSC[ShellIndex] > 0.0f ){
            figureOfMerit = pow(vectorFSC[ShellIndex]/(1.0f+vectorFSC[ShellIndex]), 0.5);
          //}else if (vectorFSC[ShellIndex] > -1.0f ){
          //    figureOfMerit = pow(-vectorFSC[ShellIndex]/(1.0f-vectorFSC[ShellIndex]), 0.5);
          }else{
              figureOfMerit = 0.0f;
          }
          if (figureOfMerit != figureOfMerit) figureOfMerit = 0.0f; //NaN check

          //figureOfMerit=pow(figureOfMerit,0.5);
          //std::cerr<<ShellIndex<<"("<<nx<<") ";

          //double map12MSA =(vectorMap1MSA[ShellIndex]+vectorMap2MSA[ShellIndex])/2.0f;
          resultMSA[ShellIndex]=vectorMap12MSA[ShellIndex]*figureOfMerit;
          if (resultMSA[ShellIndex]<0) resultMSA[ShellIndex] = 0;
          //resultMSA[ShellIndex]=vectorFSC[ShellIndex];
          //std::cerr<<vectorResolutionMSA[ShellIndex]<<","<<vectorFSC[ShellIndex]<<"\n";
			}
		}
//plot2d (vectorMap1MSA, resultMSA,  vectorResolutionMSA);

 // delete [] testImage;
 // delete [] testImage1;

  //plot2d (vectorMap1MSA, vectorMap1MSA, vectorResolutionMSA);
  return resultMSA;
}


/** 
 * Split the Otsu mean image
 */
// computeOtsuMeanImage
// **************************************
template<typename T>
std::vector<double> computeHalfMapsMsaFSCSplit ( T * map1, //!< first half map
const T * map2, //!< second half map
const T * maskVect, //!< mask (only positive values of the mask are taken into account)
std::vector<double> & vectorFSC, //!< output vector of FSC values
T * outImage, //!< Output image
T * imageFSC_at_0143, //!< Output image
T * mageFSC_at_05, //!< Output image
 int regionSizeInPixel, //!< Size of the region in Pixel
 const unsigned long int nx, //!< size of image (X)
 const unsigned long int ny, //!< size of image (Y)
 const unsigned long int nz , //!< size of image (Z)
double apix=1 //!< Pixel size in Angstrom
){

    const unsigned long int nxy = nx * ny;
    const unsigned long int nxyz = nx * ny * nz;
    
    //std::cerr<<"\n\n  regionSizeInPixel="<<regionSizeInPixel<<"\n";
//int regionSizeInPixel = 10;
    
    unsigned int nRx = ceil(nx / regionSizeInPixel);
    unsigned int nRy = ceil(ny / regionSizeInPixel);
    unsigned int nRz = ceil(nz / regionSizeInPixel);
    
    if ( nx-nRx*regionSizeInPixel <=0 ) nRx--;
    if ( ny-nRy*regionSizeInPixel <=0 ) nRy--;
    if ( nz-nRz*regionSizeInPixel <=0 ) nRz--;    
    
    T * tmpImage1 = new T [regionSizeInPixel*regionSizeInPixel*regionSizeInPixel];
    T * tmpImage2 = new T [regionSizeInPixel*regionSizeInPixel*regionSizeInPixel];
    T * tmpImageAverage = new T [regionSizeInPixel*regionSizeInPixel*regionSizeInPixel];
    T * tmpImageOut = new T [regionSizeInPixel*regionSizeInPixel*regionSizeInPixel];
    int counter = 0;
    for ( unsigned long int rZ = 0; rZ <= nRz; rZ++ ){
        unsigned long int Z1 = rZ*regionSizeInPixel;
        unsigned long int Z2 = (rZ+1)*regionSizeInPixel;
        if (Z2 > nz) Z2 = nz;
        for ( unsigned long int rY = 0; rY <= nRy; rY++ ){
            unsigned long int Y1 = rY*regionSizeInPixel;
            unsigned long int Y2 = (rY+1)*regionSizeInPixel;
            if (Y2 > ny) Y2 = ny;
            for ( unsigned long int rX = 0; rX <= nRx; rX++ ){
                unsigned long int X1 = rX*regionSizeInPixel;
                unsigned long int X2 = (rX+1)*regionSizeInPixel;
                if (X2 > nx) X2 = nx;
                
                //std::cerr<< counter << " X=(" << X1 << "," << X2 << ")  Y=(" <<Y1<< "," << Y2 << ")  Z=(" << Z1 << ", " <<Z2<< ")\n";

                for (unsigned long int kk=Z1, tmpKK=0; kk<Z2; kk++, tmpKK++){
                    for (unsigned long int jj=Y1, tmpJJ=0; jj<Y2; jj++, tmpJJ++){
                        for (unsigned long int ii=X1, tmpII=0; ii<X2; ii++, tmpII++){
                            unsigned long int indexTmp = tmpII+tmpJJ*(X2-X1)+tmpKK*(X2-X1)*(Y2-Y1);
                            unsigned long int indexOriginal = ii+jj*nx+kk*nxy;
                            tmpImage1[indexTmp]=map1[indexOriginal];
                            tmpImage2[indexTmp]=map2[indexOriginal];
                            tmpImageAverage[indexTmp]=(tmpImage1[indexTmp]+tmpImage2[indexTmp])/2.0;
                        }
                    }
                }
                

                std::vector<double> vectorFSCTmp;
                std::vector<double> vectorMSA;
                vectorMSA=computeHalfMapsMsaFSC ( tmpImage1, tmpImage2, (const T *) NULL, vectorFSCTmp, (X2-X1), (Y2-Y1), (Z2-Z1), apix);
                msaCutoff (vectorMSA, (X2-X1), (Y2-Y1), (Z2-Z1));
                std::vector<double> listResolutions = computeResolutions (vectorMSA, apix, (X2-X1), (Y2-Y1), (Z2-Z1));
                applylMsa ( tmpImageAverage, tmpImageOut, (const T *) NULL, vectorMSA, (X2-X1), (Y2-Y1), (Z2-Z1));
                double resAt05;
                double resAt0143;
                if(imageFSC_at_0143){
                    resAt0143 = computeResolution(listResolutions, vectorFSCTmp, 0.143);
                }
                if(mageFSC_at_05){
                    resAt05 = computeResolution(listResolutions, vectorFSCTmp, 0.5);
                }
        
                for (unsigned long int kk=Z1, tmpKK=0; kk<Z2; kk++, tmpKK++){
                    for (unsigned long int jj=Y1, tmpJJ=0; jj<Y2; jj++, tmpJJ++){
                        for (unsigned long int ii=X1, tmpII=0; ii<X2; ii++, tmpII++){
                            outImage[ii+jj*nx+kk*nxy]=tmpImageOut[tmpII+tmpJJ*(X2-X1)+tmpKK*(X2-X1)*(Y2-Y1)];
                            if(imageFSC_at_0143){
                                imageFSC_at_0143[ii+jj*nx+kk*nxy]=resAt0143;
                            }
                            if(mageFSC_at_05){
                                mageFSC_at_05[ii+jj*nx+kk*nxy]=resAt05;
                            }
                        }
                    }
                }                

                counter ++;
                
            }
        }
    }
    std::vector<double> processedMSA = computeMsa (outImage, maskVect, nx, ny, nz);
    
    delete [] tmpImage1;
    delete [] tmpImage2;
    delete [] tmpImageAverage;
    delete [] tmpImageOut;
    return processedMSA;
    
}






// **********************************
// **********************************
//
// computeLocResFSC
//
// **********************************
// ***********************************
template<typename T>
void  computeLocResFSC (const T * map1, //!< first half map
const T * map2, //!< second half map
const T * maskVect, //!< mask (only positive values of the mask are taken into account)
T * outMap, //!< mask (only positive values of the mask are taken into account)
const unsigned long int locResBoxSize, //!< size of image (X)
const unsigned long int nx, //!< size of image (X)
const unsigned long int ny, //!< size of image (Y)
const  unsigned long int nz, //!< size of image (Z)
double apix=1 //!< Pixel size in Angstrom
){
    //std::vector<double> & vectorFSC, //!< output vector of FSC values
    const unsigned long int nxy = nx * ny;
    const unsigned long int nxyz = nx * ny * nz;    
    std::cerr<<"computing local resolution boxSize="<< locResBoxSize<<"\n";
    for (unsigned long int ii=0;ii<nxyz; ii++){
        outMap[ii]=0;        
    }
    unsigned int radius=floor(locResBoxSize/2);
    unsigned int sizeB = 1+2*radius;
    
    T * tmpMap1 = new T [sizeB*sizeB*sizeB];
    T * tmpMap2 = new T [sizeB*sizeB*sizeB];

    for (unsigned long int kk=0, zyxIndex=0; kk<nz; kk++){
        unsigned long int Z1=kk-radius; unsigned long int Z2=kk+radius+1;
        if (Z1<0) Z1=0; if (Z2>nz) Z2=nz;
        for (unsigned long int jj=0; jj<ny; jj++){
             unsigned long int Y1=jj-radius; unsigned long int Y2=jj+radius+1;
             if (Y1<0) Y1=0; if (Y2>ny) Y2=ny;
            for (unsigned long int ii=0; ii<nx; ii++, zyxIndex++){
                unsigned long int X1=ii-radius; unsigned long int X2=ii+radius+1;
                if (X1<0) X1=0; if (X2>nx) X2=nx;

                for (unsigned long int k=Z1, tmpXYZ=0; k<Z2; k++){
                    for (unsigned long int j=Y1; j<Y2; j++){
                        for (unsigned long int i=X1; i<X2; i++, tmpXYZ++){
                            tmpMap1[tmpXYZ]=map1[i+j*nx+k*nxy];
                            tmpMap2[tmpXYZ]=map2[i+j*nx+k*nxy];
                        }
                    }
                }
                std::vector<double> computedFSC;
                std::vector<double> msaHalfmaps = computeHalfMapsMsaFSC (tmpMap1, tmpMap2, maskVect, computedFSC, X2-X1, Y2-Y1, Z2-Z1, apix);
                std::vector<double> listResolutions = computeResolutions (msaHalfmaps, apix, X2-X1, Y2-Y1, Z2-Z1);
                outMap[zyxIndex] = computeResolution(listResolutions, computedFSC, 0.143);
                //outMap[zyxIndex]=map1[zyxIndex]+map2[zyxIndex];                
            }
        }
    }
    
    delete [] tmpMap1;
    delete [] tmpMap2;
}


// **********************
// **********************
//
// compute pairwise MSA
//
// **********************
// **********************
template<typename T>
T *  computePairwiseMSA(T * mapIn, T * mapOut, const unsigned long int nx,  const unsigned long int ny, const  unsigned long int nz){

    const unsigned int sizeMAV = floor(  pow( floor( floor(pow(nz/2.0-1, 2))+ floor(pow(ny/2.0-1, 2))+ floor(pow(nx/2.0-1, 2)) ), 0.5 ) / SHELL_WIDTH ); //OK
    std::vector<double>  outputMSA;


  for (unsigned long int ii=0; ii<nx*ny*nz; ii++){
    mapOut[ii]=mapIn[ii];
  }


  std::cerr<<"FSC pairwise: ";
  for (unsigned long int ii=0; ii<nz; ii++){
    const T * map1= &mapIn[nx*ny*ii];
    unsigned long int counter = 0;
    std::vector<double>  averageMSA;
		for (unsigned int ShellIndex = 0; ShellIndex < sizeMAV; ShellIndex++){
      averageMSA.push_back(0);

		}
    for(unsigned long jj=0; jj<nz; jj++){
      if (ii != jj){
        std::cerr<<ii<<"-"<<jj<<", ";
        const T * map2=&mapIn[nx*ny*jj];
        std::vector<double>  tmpMSA;
        std::vector<double>  tmpFSC;
        tmpMSA=computeHalfMapsMsaFSC (map1, map2, (T *)NULL, tmpFSC, nx, ny, 1);
    		  for (unsigned int ShellIndex = 0; ShellIndex < sizeMAV; ShellIndex++){
            averageMSA[ShellIndex] += tmpMSA[ShellIndex];
    		  }
          counter++;
        }
      }
      std::cerr<<"counter="<<counter<<"\n";
      for (unsigned int ShellIndex = 0; ShellIndex < sizeMAV; ShellIndex++){
              averageMSA[ShellIndex] /= counter;
      }
      applylMsa ( map1, &mapOut[nx*ny*ii] , (T *)NULL, averageMSA, nx,  ny, 1);
      std::cerr<<"\n";
  }

  return mapOut;

}




// **********************
// **********************
//
// gaussian MSA
//
// **********************
// **********************
template<typename T>
std::vector<double> computeGaussianMSA(T * mapIn, double sigmaAngstrom, const unsigned long int nx,  const unsigned long int ny, const  unsigned long int nz){

  T* maskVect = NULL;
  T * outVect = new T [nx*ny*nz];
  GaussianBlurring(mapIn, outVect, sigmaAngstrom, nx, ny, nz);
  std::vector<double> outputMSA = computeMsa (outVect, maskVect, nx, ny, nz);
  delete [] outVect;
  return outputMSA;

}

std::vector<double> computeGaussianFunctionMSA(std::vector<double> msaOriginal, std::vector<double> msaProcessed, double sigma, const unsigned long int nx,  const unsigned long int ny, const  unsigned long int nz){

    int dimension=0;
    if (nx > 1) dimension++;
    if (ny > 1) dimension++;
    if (nz > 1) dimension++;
    std::vector<double> listResolutions = computeResolutions (msaOriginal, 1, nx, ny, nz);
    std::vector<double> msaBlurred=listResolutions;
    for (unsigned long int ii = 0; ii < listResolutions.size(); ii++){
        //msaBlurred[ii] = msaOriginal[ii]*exp(-0.5*listResolutions[ii]*pow(2*sigma,2)); //ITK
        //msaBlurred[ii] = msaOriginal[ii]*exp(-0.5*dimension*pow(listResolutions[ii],2)*pow(sigma,2)); //correct  
        msaBlurred[ii] = msaOriginal[ii]*exp(-0.5*pow(2*3.14*listResolutions[ii],2)*pow(sigma,2)); //correct 
    }
    return msaBlurred;
    
}



// **********************
// **********************
//
// compute Cref from FSC
//
// **********************
// **********************
std::vector<double> computCrefFromFsc (std::vector<double> vectorFSC, double sigmaBlur){
  std::vector<double> outputMSA;
  for (unsigned long int i = 0; i < vectorFSC.size(); i++){
          double figureOfMerit;
          if (vectorFSC[i] >  EPSILON_ZERO){
            figureOfMerit = pow(vectorFSC[i]/(1.0f+vectorFSC[i]), 0.5);
          }else{
              figureOfMerit =  EPSILON_ZERO;
          }
          if (figureOfMerit != figureOfMerit) figureOfMerit =  EPSILON_ZERO; //NaN check
          outputMSA.push_back(figureOfMerit);
  }
  if (sigmaBlur>0){
    outputMSA= msaBlur (outputMSA, sigmaBlur);
  }
  return outputMSA;
}


// **********************
// **********************
//
// Apply Cref
//
// **********************
// **********************
std::vector<double> applyCref (std::vector<double> inputMSA, std::vector<double> Cref, double sigmaBlur = -1){
  std::vector<double> outputMSA;
  for (unsigned long int i = 0; i < inputMSA.size(); i++){
          double resultMSA = inputMSA[i]*Cref[i];
          if (resultMSA<0.0f) resultMSA = EPSILON_ZERO;
          if (resultMSA != resultMSA) resultMSA = EPSILON_ZERO; //NaN check
          outputMSA.push_back(resultMSA);
  }
  if (sigmaBlur>0){
    outputMSA= msaBlur (outputMSA, sigmaBlur);
  }
  return outputMSA;
}





// **********************
// **********************
//
// gaussian Blurring MSA
//
// **********************
// **********************
template<typename T, typename U>
std::vector<double> GaussianBlurringFFT(T * mapIn, U * mapOut, double sigmaAngstrom, const unsigned long int nx,  const unsigned long int ny, const  unsigned long int nz){

  T * maskVect = NULL;
  GaussianBlurring(mapIn, mapOut, sigmaAngstrom, nx, ny, nz);
  std::vector<double> outputMSA = computeMsa (mapOut, maskVect, nx, ny, nz);


/*
  std::vector<double> inputMSA = computeMsa (mapIn,maskVect, nx, ny, nz);
  std::vector<double> resolutions = computeResolutions (inputMSA, apix, nx, ny, nz);
  std::vector<double> outputMSA = inputMSA;

  unsigned long int numShells = inputMSA.size();

  int nD = 3;
  if (nz <= 1) nD = 2;

    for (int ii = 1; ii < inputMSA.size(); ii++){
      double K=sqrt(SHELL_WIDTH)/(2*numShells);
      double gauss = exp(-0.5f*nD*pow(sigmaAngstrom/(apix),2.0f)*pow(ii*K,2.0f) );
      outputMSA[ii]=inputMSA[ii]*gauss;
    }
*/

  return outputMSA;

}




// **********************
// **********************
//
// Hat Function MSA
//
// **********************
// **********************
std::vector<double> hatFilter(std::vector<double> resolutions, std::vector<double> & msaOriginal, double sigma1, double sigma2, const unsigned long int nx,  const unsigned long int ny, const  unsigned long int nz){

     std::vector<double> outputMSA;
     std::cerr<<"PERFORMING HAT FILTER\n\n";
     float dimensions=2;
     if (nz>3)dimensions=3;

//     std::vector<double> msaGaussian = GaussianBlurring(mapIn, mapOut, sigmaAngstrom/apix, nx, ny, nz);
//  GaussianBlurring(mapIn, mapOut, sigmaAngstrom/apix, nx, ny, nz);
//  std::vector<double> outputMSA = computeMsa (mapOut, NULL, nx, ny, nz);

      for (unsigned long int ii = 0; ii<msaOriginal.size(); ii++){
       double x=resolutions[ii]*SHELL_WIDTH;
       double gaussian1 = exp(-(x*x*sigma1*sigma1)*dimensions/0.5f)-exp(-(x*x*sigma2*sigma2)*dimensions/0.5f);
       msaOriginal[ii]*=gaussian1;
       outputMSA.push_back(gaussian1);
      }


/*    std::vector<double> msaGaussian1D = msaGaussian;

      outputMSA.push_back(msaOriginal[0]);
      for (unsigned long int ii = 1; ii<msaOriginal.size(); ii++){
         if ( msaOriginal[ii] > EPSILON_ZERO){
           double gaussian = exp(ii*ii);//msaGaussian[ii] / msaOriginal[ii];
           double result = (gaussian * msaOriginal[ii]) + (1.0f-gaussian) * msaProcessed[ii];
           outputMSA.push_back(result);
         }else{
            outputMSA.push_back(EPSILON_ZERO);
         }
      }
*/
  return outputMSA;

}





// **********************
// **********************
//
// gaussian Blurring MSA
//
// **********************
// **********************
std::vector<double> gaussianWeightMSA(std::vector<double> msaOriginal, std::vector<double> msaProcessed, std::vector<double> msaGaussian, const unsigned long int nx,  const unsigned long int ny, const  unsigned long int nz){

     std::vector<double> outputMSA;
//     std::vector<double> msaGaussian = GaussianBlurring(mapIn, mapOut, sigmaAngstrom/apix, nx, ny, nz);
//  GaussianBlurring(mapIn, mapOut, sigmaAngstrom/apix, nx, ny, nz);
//  std::vector<double> outputMSA = computeMsa (mapOut, NULL, nx, ny, nz);
    std::vector<double> msaGaussian1D = msaGaussian;

      outputMSA.push_back(msaOriginal[0]);
      for (unsigned long int ii = 1; ii<msaOriginal.size(); ii++){
         if ( msaOriginal[ii] > EPSILON_ZERO){
           double gaussian = exp(ii*ii);//msaGaussian[ii] / msaOriginal[ii];
           double result = (gaussian * msaOriginal[ii]) + (1.0f-gaussian) * msaProcessed[ii];
           outputMSA.push_back(result);
         }else{
            outputMSA.push_back(EPSILON_ZERO);
         }
      }
  return outputMSA;

}

// **********************
// **********************
//
// gaussian Weight MSA
//
// **********************
// **********************
std::vector<double> gaussianWeightMSA(std::vector<double> msaOriginal, std::vector<double> msaProcessed, std::vector<double> resolutions,  std::vector<double> & gaussianVector, double sigma){

     std::vector<double> outputMSA;
gaussianVector=outputMSA;
//     std::vector<double> msaGaussian = GaussianBlurring(mapIn, mapOut, sigmaAngstrom/apix, nx, ny, nz);
//  GaussianBlurring(mapIn, mapOut, sigmaAngstrom/apix, nx, ny, nz);
//  std::vector<double> outputMSA = computeMsa (mapOut, NULL, nx, ny, nz);
   // std::vector<double> msaGaussian1D = msaGaussian;

      outputMSA.push_back(msaOriginal[0]);
      gaussianVector.push_back(1);
      for (unsigned long int ii = 1; ii<msaOriginal.size(); ii++){
         if ( msaOriginal[ii] > EPSILON_ZERO){
           double gaussian = exp(-0.5f*(resolutions[ii]*resolutions[ii])*sigma*sigma*sigma*sigma);//msaGaussian[ii] / msaOriginal[ii];
            gaussianVector.push_back(gaussian);
           double result = (gaussian * msaOriginal[ii]) + (1.0f-gaussian) * msaProcessed[ii];
           outputMSA.push_back(result);
         }else{
            outputMSA.push_back(EPSILON_ZERO);
         }
      }
  return outputMSA;

}


// **********************
// **********************
//
// msaStretch
//
// **********************
// **********************
/*std::vector<double> msaStretch(std::vector<double> msaOriginal, std::vector<double> msaProcessed, const unsigned long int nx,  const unsigned long int ny, const  unsigned long int nz){

     std::vector<double> outputMSA;
     double stretchFactor;
     if (msaOriginal[0]>EPSILON_ZERO){
        stretchFactor = msaOriginal[0]/msaProcessed[0];
     }else if (msaOriginal[1]>EPSILON_ZERO){
        stretchFactor = msaOriginal[1]/msaProcessed[1];
     } else{
      stretchFactor = 1.0f;
    }
      for (unsigned long int ii = 0; ii<msaOriginal.size(); ii++){
           double result = msaProcessed[ii]*stretchFactor;
           outputMSA.push_back(result);
      }
  return outputMSA;

}
*/


// **********************
// **********************
//
// msaShift
//
// **********************
// **********************
/*std::vector<double> msaShift(std::vector<double> msaOriginal, std::vector<double> msaProcessed){

     std::vector<double> outputMSA;
     double shiftFactor =  msaProcessed[0] - msaOriginal[0];
      for (unsigned long int ii = 0; ii<msaOriginal.size(); ii++){
           double result = msaProcessed[ii] - shiftFactor;
           outputMSA.push_back(result);
      }
    return outputMSA;

}*/


// **********************
// **********************
//
// msaRatio
//
// **********************
// **********************
std::vector<double> msaRatio(std::vector<double> msaOriginal, std::vector<double> msaProcessed){

     std::vector<double> outputMSA;
      for (unsigned long int ii = 0; ii<msaOriginal.size(); ii++){
           double result = 0;
           if (msaProcessed[ii] != 0 ){
               result = msaOriginal[ii] / msaProcessed[ii];
           }
           outputMSA.push_back(result);
      }
    return outputMSA;

}


#endif
