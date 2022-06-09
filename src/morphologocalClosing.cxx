/*
    (C) Mauro Maiorca, 2020. Birkbeck College, London University.

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

#include <iostream>
#include <complex>
#include <fstream>
#include <vector>
#include <iterator>     // std::back_inserter
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <sstream>
#include <string>
#include <cstdio>
#include <cmath>
#include <iomanip>      // std::setprecision
#include <arpa/inet.h> //used for htonl (big/little endian check)
#include "mrcIO.h"


//ITK
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkCastImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"


#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryThinningImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkConnectedComponentImageFilter.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef TWO_PI
#define TWO_PI 6.28318530717958647
#endif

/*
#include "mrcIO.h"
#include "scoringFunctionsLib.h"
#include "CsvStarReadWriteAnalyseLibs.h"
#include "eigenAnalysisLibs.h"
#include "imageProcessingLib.h"
*/


typedef float WorkingPixelType;
typedef float OutputPixelType;
const unsigned int Dimension = 3;



/* ******************************************
 *  USAGE
 ***************************************** */
void usage(  char ** argv ){
    std::cerr<<"\n";
    std::cerr<<"Usage: " << argv[0] << "  MaskInFilename.mrc MaskOutFilename.mrc radiusSize --separateRegions [preSeparationClosingPixelRadius=1 preSeparationErosionPixelRadius=0 separatedRegionsOutMap.mrc]\n";
    std::cerr<<"\n";
    exit(1);
}


// *************************************
//
// retrieveInputParameters
// *************************************
typedef struct inputParametersType {
    char * MaskInFilename;
    char * MaskOutFilename;
    char * separatedRegionsOutMap;
    int radiusSize;
    int preSeparationClosingPixelRadius;
    int preSeparationErosionPixelRadius;
    bool separateRegions;
}inputParametersType;



void retrieveInputParameters(inputParametersType * parameters, int argc, char** argv){
    parameters->MaskInFilename= NULL;
    parameters->MaskOutFilename= NULL;
    parameters->separatedRegionsOutMap=NULL;
    parameters->radiusSize= 1;
    parameters->preSeparationClosingPixelRadius = 1;
    parameters->preSeparationErosionPixelRadius = 1;
    parameters->separateRegions = false;

    if (argc >= 2){
     parameters->MaskInFilename=argv[1];
    }

    if (argc >= 3){
     parameters->MaskOutFilename=argv[2];
    }

    if (argc >= 4){
     parameters->radiusSize=atoi(argv[3]);
    }

    for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--separateRegions") ){
          parameters->separateRegions = true;
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->preSeparationClosingPixelRadius=atoi(argv[jj]);
            if (idx==1) parameters->preSeparationErosionPixelRadius=atoi(argv[jj]);
            if (idx==2) parameters->separatedRegionsOutMap=argv[jj];
          }
          ii+=idx;
        }
    }

    
    if (!parameters->MaskInFilename || !parameters->MaskOutFilename){
      usage(argv);
    }

    
  }



// **********************************
//
//    INT MAIN
//
// **********************************
int main( int argc, char **argv ){
	  inputParametersType parameters;
	  retrieveInputParameters(&parameters, argc, argv);	
		
	  MRCHeader imageHeader;
	  readHeaderMrc(parameters.MaskInFilename, imageHeader);
	  unsigned long int nx=imageHeader.nx, ny=imageHeader.ny, nz=imageHeader.nz;
	  unsigned long int nxyz=nx*ny*nz;

    
    
	  WorkingPixelType * I = new WorkingPixelType [nxyz];
	  WorkingPixelType * I0 = new WorkingPixelType [nxyz];
	  readMrcImage(parameters.MaskInFilename, I, imageHeader);
	  typedef float OutputPixelType;
	  typedef itk::Image< WorkingPixelType, Dimension > workingImageType;
	  typedef itk::Image< OutputPixelType, Dimension > outputImageType;
	  typedef itk::ImageRegionIterator< workingImageType> IteratorType;
      typedef itk::BinaryBallStructuringElement<workingImageType::PixelType, workingImageType::ImageDimension> StructuringElementType;
      typedef int labelPixelType;
      typedef itk::Image< labelPixelType, Dimension > labelImageType;
    
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
	  workingImageIt.GoToBegin();
	  for (unsigned long int ii=0;ii<nxyz&&!workingImageIt.IsAtEnd();ii++,++workingImageIt){
	        if (I[ii]>=0.01){
        	  workingImageIt.Set(1);
        	  I0[ii]=1;
        	}else{
        	   workingImageIt.Set(0);
        	   I0[ii]=0;
        	}
	  }
  
      if (!parameters.separateRegions || (parameters.separateRegions && parameters.preSeparationClosingPixelRadius>0) ){
        

	          StructuringElementType structuringElement;
	          if (parameters.separateRegions && parameters.preSeparationClosingPixelRadius>0){
	              structuringElement.SetRadius(parameters.preSeparationClosingPixelRadius);
	          } else {
	              structuringElement.SetRadius(parameters.radiusSize);
	          }
	          structuringElement.CreateStructuringElement();
	          typedef itk::BinaryMorphologicalClosingImageFilter <workingImageType, workingImageType, StructuringElementType> BinaryMorphologicalClosingImageFilterType;
	          BinaryMorphologicalClosingImageFilterType::Pointer closingFilter = BinaryMorphologicalClosingImageFilterType::New();
	          closingFilter->SetInput(workingImage);
	          closingFilter->SetForegroundValue (1);
	          closingFilter->SetKernel(structuringElement);
	          closingFilter->Update();
	          IteratorType outImageIt(closingFilter->GetOutput(), workingImageRegion );
	          outImageIt.GoToBegin();
	          workingImageIt.GoToBegin();
	          for (unsigned long int ii=0;ii<nxyz&&!outImageIt.IsAtEnd();ii++,++outImageIt, ++workingImageIt){
	              if (I[ii]>=0.01f){
                      //I[ii]=1.0f;
                	  workingImageIt.Set(1);
                      }else if (outImageIt.Get() >= 0.01f ) {
	                      I[ii]=1;
	                      workingImageIt.Set(1);
                      } else{
                          I[ii]=0;
                          workingImageIt.Set(0);
                      }
              }
        }
    /*
	  if (parameters.separateRegions){
        	WorkingPixelType * labelI = new WorkingPixelType [nxyz];
	        WorkingPixelType * outI = new WorkingPixelType [nxyz];        	
        	workingImage->Update();
          //
          StructuringElementType structuringElement2;
          if (parameters.preSeparationErosionPixelRadius > 0){
              std::cerr<<"radius="<<parameters.preSeparationErosionPixelRadius<<"\n";
              structuringElement2.SetRadius(parameters.preSeparationErosionPixelRadius);
              typedef itk::BinaryErodeImageFilter <workingImageType, workingImageType, StructuringElementType> BinaryMorphologicalErosionImageFilterType;
              BinaryMorphologicalErosionImageFilterType::Pointer erosionFilter = BinaryMorphologicalErosionImageFilterType::New();
              erosionFilter->SetInput(workingImage);
              erosionFilter->SetForegroundValue (1);
              erosionFilter->SetKernel(structuringElement2);
              erosionFilter->Update();
              typedef itk::ImageRegionIterator< labelImageType> IteratorLabelType;
              IteratorType erosionImageIt(erosionFilter->GetOutput(), workingImageRegion );
              
               workingImageIt.GoToBegin();
              erosionImageIt.GoToBegin();
              
              for (unsigned long int ii=0;!erosionImageIt.IsAtEnd();++ii, ++erosionImageIt, ++workingImageIt){
                  if (erosionImageIt.Get()>0){
                    workingImageIt.Set(1);
                    outI[ii]=1;
                  }else{
                      workingImageIt.Set(0);
                      outI[ii]=0;
                  }
              }
              writeMrcImage("eroded.mrc", outI, imageHeader);
              workingImage->Update();
          }
          


        	typedef itk::ConnectedComponentImageFilter <workingImageType, labelImageType > ConnectedComponentImageFilterType;
                ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New ();
                connected->SetFullyConnected (true);
                connected->SetInput(workingImage);
                connected->Update();
                
                typedef itk::ImageRegionIterator< labelImageType> IteratorLabelType;
                IteratorLabelType labelsImageIt(connected->GetOutput(), (connected->GetOutput())->GetLargestPossibleRegion() );
                labelsImageIt.GoToBegin();
	        if (parameters.separatedRegionsOutMap){
                   for (unsigned long int ii=0;ii<nxyz&&!labelsImageIt.IsAtEnd();ii++,++labelsImageIt){
        	        if (I0[ii]>0.1){
        	          labelI[ii]=labelsImageIt.Get();
        	        }else{
         	          labelI[ii]=0;
        	        }
	           }
	           writeMrcImage(parameters.separatedRegionsOutMap, labelI, imageHeader);
	        }
	        int NumSeparatedRegions = connected->GetObjectCount();
	        std::cerr<<"Num Separated Regions="<<NumSeparatedRegions<<"\n";
	        for ( int ii=0; ii<nxyz; ii++){
        	    outI[ii]=0;
	        }

	        for ( int kk=0; kk<NumSeparatedRegions; kk++){
	             
	          StructuringElementType structuringElement;
	          structuringElement.SetRadius(parameters.radiusSize);
	          structuringElement.CreateStructuringElement();
	          typedef itk::BinaryMorphologicalClosingImageFilter <labelImageType, labelImageType, StructuringElementType> BinaryMorphologicalClosingImageFilterType;
	          BinaryMorphologicalClosingImageFilterType::Pointer closingFilter = BinaryMorphologicalClosingImageFilterType::New();
	          closingFilter->SetInput( connected->GetOutput() );
	          closingFilter->SetForegroundValue (kk+1);
	          closingFilter->SetKernel(structuringElement);
	          closingFilter->Update();

                
	          IteratorLabelType outImageIt(closingFilter->GetOutput(), workingImageRegion );
	          outImageIt.GoToBegin();
	          for ( int ii=0; ii<nxyz; ii++, ++outImageIt){
        	    if (outImageIt.Get()>0){
        	       outI[ii]=kk;
        	    }
	          }
	        }
            for ( int ii=0; ii<nxyz; ii++){
        	   I[ii]=outI[ii];
	        }
	        delete []outI;
	        delete []labelI;
	  }
*/
     writeMrcImage(parameters.MaskOutFilename, I, imageHeader);
          delete [] I;
          delete [] I0;
}






