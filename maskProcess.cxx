/*
    (C) Mauro Maiorca, 2022

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
#include "imageProcessingLib.h"

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
    //std::cerr<<"Usage: " << argv[0] << "  MaskInFilename.mrc MaskOutFilename.mrc radiusSize --separateRegions [preSeparationClosingPixelRadius=1 preSeparationErosionPixelRadius=0 separatedRegionsOutMap.mrc]\n";
    std::cerr<<"Usage: " << argv[0] << "\n";
    std::cerr<<"             --i MaskInFilename.mrc \n";
    std::cerr<<"             --o MaskOutFilename.mrc \n";
    std::cerr<<"             --closing closingRadiusSize \n";
    std::cerr<<"             --opening openingRadiusSize \n";
    std::cerr<<"             --dilate dilateRadiusSize \n";
    std::cerr<<"             --erode erodeRadiusSize \n";

    std::cerr<<"             --getLargestComponent \n";

    
    std::cerr<<"             --contrastDifference inputContrastMap.mrc outputContrastFile.txt ContrastRadiusSize[=1] ContrastRadiusStep [=1] \n";
    std::cerr<<"             --maskedMean inputGrayscaleMap.mrc outputContrastFile.txt ContrastRadiusSize[=1] ContrastRadiusStep [=1] \n";

//    std::cerr<<"             --blurEdge blurEdgeRadiusSize \n";    
    std::cerr<<"      --h (help)\n";
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

    int closingRadiusSize;
    int blurEdgeRadiusSize;
    int openingRadiusSize;
    int dilateRadiusSize;
    int erodeRadiusSize;


    char * inputContrastMap;
    char * outputContrastFile;
    int ContrastRadiusSize;
    int ContrastRadiusStep;

    char * inputGrayscaleMap;
    bool getLargestComponent;




}inputParametersType;



void retrieveInputParameters(inputParametersType * parameters, int argc, char** argv){
    parameters->MaskInFilename= NULL;
    parameters->MaskOutFilename= NULL;
    parameters->separatedRegionsOutMap=NULL;
    parameters->radiusSize= 1;
    parameters->preSeparationClosingPixelRadius = 1;
    parameters->preSeparationErosionPixelRadius = 1;
    parameters->separateRegions = false;
    parameters-> closingRadiusSize = -1;
    parameters-> openingRadiusSize = -1;
    parameters-> dilateRadiusSize = -1;
    parameters-> erodeRadiusSize = -1;
    parameters-> blurEdgeRadiusSize = -1;

    parameters-> inputContrastMap = NULL;
    parameters-> outputContrastFile = NULL;
    parameters-> ContrastRadiusSize=1;
    parameters-> ContrastRadiusStep=1;

    parameters-> inputGrayscaleMap = NULL;
    
    parameters->getLargestComponent = false;

//    if (argc >= 2){
//     parameters->MaskInFilename=argv[1];
//    }
  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--i") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->MaskInFilename=argv[jj];
          }
          ii+=idx;
        }
  }


  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--o") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->MaskOutFilename=argv[jj];
          }
          ii+=idx;
        }
  }


  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--closing") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->closingRadiusSize=atoi(argv[jj]);
          }
          ii+=idx;
        }
  }


  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--opening") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->openingRadiusSize=atoi(argv[jj]);
          }
          ii+=idx;
        }
  }

    for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--dilate") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->dilateRadiusSize=atoi(argv[jj]);
          }
          ii+=idx;
        }
  }

    for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--erode") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->erodeRadiusSize=atoi(argv[jj]);
          }
          ii+=idx;
        }
  }

  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--blurEdge") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->blurEdgeRadiusSize=atoi(argv[jj]);
          }
          ii+=idx;
        }
  }

  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--contrastDifference") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->inputContrastMap=argv[jj];
            if (idx==1) parameters->outputContrastFile=argv[jj];
            if (idx==2) parameters->ContrastRadiusSize=atoi(argv[jj]);
            if (idx==3) parameters->ContrastRadiusStep=atoi(argv[jj]);
          }
          ii+=idx;
        }
  }

    for (unsigned int ii=1;ii<argc;ii++){
          std::string optionStr(argv[ii]);
          if ( !optionStr.compare("--getLargestComponent") ){
            int idx=0;
            parameters->getLargestComponent=true;
          }
    }
    
    

  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--maskedMean") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->inputGrayscaleMap=argv[jj];
          }
          ii+=idx;
        }
  }


//    if (argc >= 4){
//     parameters->radiusSize=atoi(argv[3]);
//    }
//    closingRadiusSize



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

    for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if (!optionStr.compare("--h") ){
            usage(argv);
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
    WorkingPixelType * M = new WorkingPixelType [nxyz];
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

if (parameters.closingRadiusSize > 0){
	          StructuringElementType structuringElement;
	          //if (parameters.separateRegions && parameters.preSeparationClosingPixelRadius>0){
	          //    structuringElement.SetRadius(parameters.preSeparationClosingPixelRadius);
	          //} else {
	              structuringElement.SetRadius(parameters.closingRadiusSize);
	          //}
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
              writeMrcImage(parameters.MaskOutFilename, I, imageHeader);

} 

else if (parameters.openingRadiusSize > 0){
	          StructuringElementType structuringElement;
	          //if (parameters.separateRegions && parameters.preSeparationClosingPixelRadius>0){
	          //    structuringElement.SetRadius(parameters.preSeparationClosingPixelRadius);
	          //} else {
	              structuringElement.SetRadius(parameters.openingRadiusSize);
	          //}
	          structuringElement.CreateStructuringElement();
	          typedef itk::BinaryMorphologicalOpeningImageFilter <workingImageType, workingImageType, StructuringElementType> BinaryMorphologicalOpeningImageFilterType;
	          BinaryMorphologicalOpeningImageFilterType::Pointer openingFilter = BinaryMorphologicalOpeningImageFilterType::New();
	          openingFilter->SetInput(workingImage);
	          openingFilter->SetForegroundValue (1);
	          openingFilter->SetKernel(structuringElement);
	          openingFilter->Update();
	          IteratorType outImageIt(openingFilter->GetOutput(), workingImageRegion );
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
              writeMrcImage(parameters.MaskOutFilename, I, imageHeader);
} 

else if (parameters.dilateRadiusSize > 0){
	          StructuringElementType structuringElement;
	          //if (parameters.separateRegions && parameters.preSeparationClosingPixelRadius>0){
	          //    structuringElement.SetRadius(parameters.preSeparationClosingPixelRadius);
	          //} else {
	              structuringElement.SetRadius(parameters.dilateRadiusSize);
	          //}
	          structuringElement.CreateStructuringElement();
	          typedef itk::BinaryDilateImageFilter <workingImageType, workingImageType, StructuringElementType> DilateImageFilterType;
	          DilateImageFilterType::Pointer dilateFilter = DilateImageFilterType::New();
	          dilateFilter->SetInput(workingImage);
	          dilateFilter->SetForegroundValue (1);
	          dilateFilter->SetKernel(structuringElement);
	          dilateFilter->Update();
	          IteratorType outImageIt(dilateFilter->GetOutput(), workingImageRegion );
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
              writeMrcImage(parameters.MaskOutFilename, I, imageHeader);

} 

else if (parameters.erodeRadiusSize > 0){
	          StructuringElementType structuringElement;
	          //if (parameters.separateRegions && parameters.preSeparationClosingPixelRadius>0){
	          //    structuringElement.SetRadius(parameters.preSeparationClosingPixelRadius);
	          //} else {
	              structuringElement.SetRadius(parameters.erodeRadiusSize);
	          //}
	          structuringElement.CreateStructuringElement();
	          typedef itk::BinaryErodeImageFilter <workingImageType, workingImageType, StructuringElementType> ErodeImageFilterType;
	          ErodeImageFilterType::Pointer erodeFilter = ErodeImageFilterType::New();
	          erodeFilter->SetInput(workingImage);
	          erodeFilter->SetForegroundValue (1);
            erodeFilter->SetBackgroundValue (0);
	          erodeFilter->SetKernel(structuringElement);
	          erodeFilter->Update();
	          IteratorType outImageIt(erodeFilter->GetOutput(), workingImageRegion );
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
              writeMrcImage(parameters.MaskOutFilename, I, imageHeader);

} 

else if (parameters.getLargestComponent){
    
        typedef itk::ConnectedComponentImageFilter <workingImageType, labelImageType > ConnectedComponentImageFilterType;
        ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New ();
        connected->SetFullyConnected (true);
        connected->SetInput(workingImage);
        connected->Update();
        int NumSeparatedRegions = connected->GetObjectCount();
        //std::vector<long int> regionCounter(NumSeparatedRegions+1,0);
        long int * regionCounter = new long int [NumSeparatedRegions+1];
        for (unsigned long int ii=1; ii<NumSeparatedRegions+1; ii++){
            regionCounter[ii]=0;
        }
        
    
        typedef itk::ImageRegionIterator< labelImageType> IteratorLabelType;
        IteratorLabelType labelsImageIt(connected->GetOutput(), (connected->GetOutput())->GetLargestPossibleRegion() );
        labelsImageIt.GoToBegin();
        for (unsigned long int ii=0;ii<nxyz&&!labelsImageIt.IsAtEnd();ii++,++labelsImageIt){
            if (labelsImageIt.Get()>0){
              regionCounter[labelsImageIt.Get()]++;
              I[ii]=labelsImageIt.Get();
            }
        }
        //std::cerr<<"NumSeparatedRegions="<<NumSeparatedRegions <<"\n";

    
    
        //get the largest region
        int targetRegion = 0;
        for (unsigned long int ii=1; ii<NumSeparatedRegions; ii++){
            if ( regionCounter[ii] > regionCounter[targetRegion] ){
                targetRegion=ii;
            }
        }
    
    //std::cerr<<"targetRegion="<<targetRegion <<"\n";

        for (unsigned long int ii=0; ii<nxyz ; ii++ ){
            if ( I[ii] == targetRegion ){
                I[ii] = 1;
            }else{
                I[ii]=0;
            }
        }

    /*
        //labelsImageIt.GoToBegin();
        for (unsigned long int ii=0;ii<nxyz && !labelsImageIt.IsAtEnd();++ii,++labelsImageIt){
            if (labelsImageIt.Get() == targetRegion){
              //I[ii] = 1;
            }else{
              //I[ii]=0;
            }
        }
    */
    delete [] regionCounter;
        writeMrcImage(parameters.MaskOutFilename, I, imageHeader);
    //std::cerr<<"Num Separated Regions="<<NumSeparatedRegions<<"\n";
    //for ( int ii=0; ii<nxyz; ii++){
    //    outI[ii]=0;
   // }
    
}
    
    
else if (parameters.blurEdgeRadiusSize > 0){
//  for ()
//  void blurMaskEdge(I, M, I, double sigma, int kernelLength , const unsigned long int nx, const unsigned long int ny){
//  writeMrcImage(parameters.MaskOutFilename, I, imageHeader);
//    std::cerr<<"             --blurEdge blurEdgeRadiusSize \n";    
}else if (parameters.inputContrastMap && parameters.outputContrastFile){

/*
	    for (unsigned long int ii=0;ii<nxyz&&!workingImageIt.IsAtEnd();ii++,++workingImageIt){
	        if (I[ii]>=0.01){
        	  workingImageIt.Set(1);
        	  I0[ii]=1;
        	}else{
        	   workingImageIt.Set(0);
        	   I0[ii]=0;
        	}
	    }
*/

	          StructuringElementType structuringElement;
	          structuringElement.SetRadius(parameters.ContrastRadiusSize);
	          structuringElement.CreateStructuringElement();
	          typedef itk::BinaryErodeImageFilter <workingImageType, workingImageType, StructuringElementType> BinaryErodeFilterType;
	          BinaryErodeFilterType::Pointer erodeFilter = BinaryErodeFilterType::New();

	          typedef itk::BinaryDilateImageFilter <workingImageType, workingImageType, StructuringElementType> BinaryDilateImageFilter;
	          BinaryDilateImageFilter::Pointer dilateFilter = BinaryDilateImageFilter::New();

            erodeFilter->SetInput(workingImage);
	          erodeFilter->SetForegroundValue (1);
	          erodeFilter->SetKernel(structuringElement);
	          erodeFilter->Update();
	          IteratorType erodeImageIt(erodeFilter->GetOutput(), workingImageRegion );


            dilateFilter->SetInput(workingImage);
	          dilateFilter->SetForegroundValue (1);
	          dilateFilter->SetKernel(structuringElement);
	          dilateFilter->Update();
	          IteratorType dilateImageIt(dilateFilter->GetOutput(), workingImageRegion );

            WorkingPixelType * M1 = new WorkingPixelType [nxyz];
            WorkingPixelType * M2 = new WorkingPixelType [nxyz];

            WorkingPixelType * dilatedI = new WorkingPixelType [nxyz];
            WorkingPixelType * erodedI = new WorkingPixelType [nxyz];
            WorkingPixelType * boundary = new WorkingPixelType [nxyz];


            erodeImageIt.GoToBegin();
            dilateImageIt.GoToBegin();
            workingImageIt.GoToBegin();

            double meanOuter = 0;
            double meanInner = 0;
            double meanGlobal = 0;
            double counterInner = 0;
            double counterOuter = 0;
      	    readMrcImage(parameters.inputContrastMap, I, imageHeader);

            for (unsigned long int ii=0;ii<nxyz&&!dilateImageIt.IsAtEnd();ii++,++dilateImageIt,++erodeImageIt,++workingImageIt){
                boundary[ii]=0;
                if (erodeImageIt.Get() > 0.5){
                  erodedI[ii]=1;
                }else{
                  erodedI[ii]=0;
                }
                if (dilateImageIt.Get() > 0.5){
                  dilatedI[ii]=1;
                }else{
                  dilatedI[ii]=0;
                }


                if (workingImageIt.Get() < 0.5){
                    M1[ii]=0;
                } else if (erodeImageIt.Get() > 0.5){
                 M1[ii]=0;
                }else{
                  M1[ii]=1;
                  boundary[ii]=1;
                    meanInner+=I[ii];
                    counterInner++;
                }

                if (workingImageIt.Get() < 0.5 && dilateImageIt.Get() > 0.5){
                    M2[ii]=1;
                    boundary[ii]=2;
                    meanOuter+=I[ii];
                    counterOuter++;
                }else{
                  M2[ii]=0;
                }
            }
            writeMrcImage("erodedI.mrc", erodedI, imageHeader);
            writeMrcImage("dilatedI.mrc", dilatedI, imageHeader);
            writeMrcImage("boundaryI.mrc", boundary, imageHeader);
            if (counterInner > 0){
             meanInner=meanInner/ counterInner;
            }
            if (counterOuter > 0){
             meanOuter=meanOuter/ counterOuter;
            }
            meanGlobal+=meanInner+meanOuter;
            if (counterInner+counterOuter){
              meanGlobal=meanGlobal/ (counterInner+counterOuter);
            }
            //std::cerr<<"inner,outer,diff,percentage="<<
            std::cout<<meanInner<<","<<meanOuter<<","<< meanInner-meanOuter<<","<<(meanInner+meanOuter)/(meanInner-meanOuter) <<"\n";

            delete [] dilatedI;
            delete [] erodedI;
            delete [] boundary;
            delete [] M1;
            delete [] M2;



/*


        	  WorkingPixelType * I = new WorkingPixelType [nxyz];
	  WorkingPixelType * I0 = new WorkingPixelType [nxyz];
    WorkingPixelType * M = new WorkingPixelType [nxyz];
	  readMrcImage(parameters.MaskInFilename, I, imageHeader);
*/

/*
	          closingFilter->SetInput(workingImage);
	          closingFilter->SetForegroundValue (1);
	          closingFilter->SetKernel(structuringElement);
	          closingFilter->Update();
	          IteratorType outImageIt(closingFilter->GetOutput(), workingImageRegion );

*/


} else if (parameters.inputGrayscaleMap) {

    WorkingPixelType * gI = new WorkingPixelType [nxyz];
	  readMrcImage(parameters.inputGrayscaleMap, gI, imageHeader);

    unsigned long int maxLabel = 0;
    for (unsigned long int ii = 0; ii < nxyz; ii++){
      unsigned long int currentLabel = I[ii];
       if ( currentLabel > maxLabel ){
         maxLabel = currentLabel;
         //std::cerr<<"-> FOUND "<<maxLabel<<"\n";
       }
    }
    long double * labelsMean = new long double [maxLabel+1];
    unsigned long int * labelsCounter = new unsigned long int [maxLabel+1];
    for (unsigned long int ii = 0; ii < maxLabel+1; ii++){
      labelsMean[ii]=0;
      labelsCounter[ii]=0;
      //std::cerr<<"-> label "<<ii<<"\n";
    }

    for (unsigned long int ii = 0; ii < nxyz; ii++){
      unsigned long int currentLabel =int(I[ii]);
      M[ii]=currentLabel;
      labelsMean[currentLabel]+=gI[ii];
      labelsCounter[currentLabel]++;
    }

    for (unsigned long int ii = 0; ii < maxLabel+1; ii++){
      if ( labelsCounter[ii] > 0 ){
        labelsMean[ii]=labelsMean[ii]/((double)(labelsCounter[ii]));
      }
    }

    for (unsigned long int ii = 0; ii < nxyz; ii++){
      unsigned long int currentLabel =int(I[ii]);
      M[ii]=labelsMean[currentLabel];
    }

    //for (unsigned long int ii = 0; ii < maxLabel+1; ii++){
    //  std::cerr<<"label "<<ii<<"  mean="<<labelsMean[ii]<<"  counter=" <<labelsCounter[ii]<<"\n";
    //}


    writeMrcImage(parameters.MaskOutFilename, M, imageHeader);

    delete [] labelsMean;
    delete [] labelsCounter;
    delete [] gI;



} else{
  writeMrcImage(parameters.MaskOutFilename, I, imageHeader);
}



/*  
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
          delete [] I;
          delete [] I0;
          delete [] M;
}






