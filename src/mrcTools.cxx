#include <iostream>
#include <complex>
#include <fstream>
#include <list>
#include <iterator>     // std::back_inserter
#include <cstdlib>
#include <ctime>

//ITK
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileReader.h"
#include "itkMRCImageIO.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkMeanSquaresImageToImageMetric.h"
#include "itkTranslationTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkDiscreteGaussianDerivativeImageFilter.h"
#include "itkDerivativeImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include <itkThresholdMaximumConnectedComponentsImageFilter.h>
#include "itkConnectedComponentImageFilter.h"
#include "itkOtsuMultipleThresholdsImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryThinningImageFilter.h"
#include "itkMultiplyImageFilter.h"


//
//#include "fourierFunctions.h"
//#include "genericFilters.h"
#include "mrcIO.h"

  typedef float WorkingPixelType;
  typedef float OutputPixelType;
  const unsigned int Dimension = 3;
  typedef itk::Image< WorkingPixelType, Dimension > workingImageType;
  typedef itk::Image< OutputPixelType, Dimension > outputImageType;
  typedef itk::ImageRegionIterator< workingImageType> IteratorType;


// *************************************
//
// inputParametersType
// **************************************
typedef enum {COPY_HEADER, COPY_IMAGE} ActionType;
typedef struct inputParametersType {
	char * inputFilename;
	char * outputFilename;
  ActionType Action;
}inputParametersType;


/* ******************************************
 *  USAGE
 ***************************************** */
void usage(  char ** argv ){
    std::cerr<<"\n";
    std::cerr<<"Usage: " << argv[0] << "  InputFile.mrc outputFile.mrc  [Options]\n";
    std::cerr<<"     Options: \n";
    std::cerr<<"              --replaceHeader    (just replace the header)\n";
    std::cerr<<"              --copyImage  (Default)\n";
    std::cerr<<"\n";
  exit(1);
}



// *************************************
//
// retrieveInputParameters
// *************************************
void retrieveInputParameters(inputParametersType * parameters, int argc, char** argv){
	if ( argc < 2)
		usage(argv);
	parameters->inputFilename = argv[1];
	parameters->outputFilename = argv[2];
  bool replaceHeader = false;


  //check for replaceHeader
  for (unsigned int ii=3;ii<argc&&!replaceHeader;ii++){
    std::string optionStr(argv[ii]);
    if (!optionStr.compare("--replaceHeader")){
      std::cerr<<"--replaceHeader  \n";
        replaceHeader=true;
    }
  }

  if (replaceHeader){
    parameters->Action=COPY_HEADER;
  }else{
    parameters->Action=COPY_IMAGE;
  }

}




// **********************************
//
//    INT MAIN
//
// **********************************
int main( int argc, char **argv ){


//if (argc < 5){
//    usage( argv );
//  }
	inputParametersType parameters;
	retrieveInputParameters(&parameters, argc, argv);


  //check extension to see if input/output are mrc files
  bool isInMrcFile = true;
  bool isOutMrcFile = true;


  char extension;
  memcpy ( &extension, &((argv[1])[strlen(parameters.inputFilename)-1]), 1 );
  if (extension=='k' || extension=='f' || extension=='g' || extension=='m' || extension=='d'  || extension=='i') isInMrcFile = false;

  memcpy ( &extension, &((argv[2])[strlen(parameters.outputFilename)-1]), 1 );
  if (extension=='k' || extension=='f' || extension=='g' || extension=='m' || extension=='d'  || extension=='i') isOutMrcFile = false;

  if (parameters.Action==COPY_HEADER){
      std::cerr<<"COPY THE HEADER FILES \n";
//      copyHeaderMrcImage(parameters.inputFilename, parameters.outputFilename);

      MRCHeader HeaderReferenceMap;
      readHeaderMrc(parameters.inputFilename, HeaderReferenceMap);
      copyHeaderMrcImage(HeaderReferenceMap, parameters.outputFilename);

      //return;

  }else{
        std::cerr<<"COPY THE FILES \n";
        unsigned long int i = 0;
        unsigned long int kk=0;

        typedef itk::ImageFileReader< workingImageType >  ReaderType;
        ReaderType::Pointer readerOriginal = ReaderType::New();
        if (isInMrcFile) readerOriginal->SetImageIO( itk::MRCImageIO::New() );
        readerOriginal->SetFileName( parameters.inputFilename );
        readerOriginal->Update();

      	workingImageType::Pointer underProcessingImage = readerOriginal->GetOutput();
      	int nx=( underProcessingImage->GetLargestPossibleRegion() ).GetSize()[0];
      	int ny=( underProcessingImage->GetLargestPossibleRegion() ).GetSize()[1];
      	int nz=( underProcessingImage->GetLargestPossibleRegion() ).GetSize()[2];
      	if (nz<1) nz=1;




      	unsigned long int nxy = nx*ny;
      	unsigned long int nxyz = nx*ny*nz;
      	bool isImage2D = false;
      	if ( nz <= 1 ) isImage2D = true;

      //	std::cerr<<"Image Size: nx="<<nx<<"  ny="<<ny<<"  nz="<<nz<<"\n";

      	kk=0;
      	WorkingPixelType * inputImageVector = new WorkingPixelType [nxyz];
      	WorkingPixelType * underProcessingVector = new WorkingPixelType [nxyz];

      /*
      MRCHeader header;
      readMrcImage(parameters.inputFilename, inputImageVector, header);
      writeMrcImage(parameters.inputFilename, inputImageVector, header);
      */



      //std::cerr<<"header size="<<sizeof(header)<<"\n";
      /*std::cerr<<"MAURO Image Size: nx="<<header.nx<<"  ny="<<header.ny<<"  nz="<<header.nz<<"\n";
      std::cerr<<"\tNum of columns, rows,sections:     "<<header.nx<<",  "<<header.ny<<",  "<<header.nz<<"\n";
      std::cerr<<"\tMode:                              "<<header.mode<<"\n";
      std::cerr<<"\tNum of First column, row, section: "<<header.nxstart<<",  "<< header.nystart<<",  "<< header.nzstart<<"\n";
      std::cerr<<"\tNum of intervals along x, y, z:    "<<header.mx<<",  "<<header.my<<",  "<<header.mz<<"\n";
      std::cerr<<"\tCell dimensions in angstroms:      "<<header.cella[0]<<",  "<< header.cella[1]<<",  "<<header.cella[2]<<"\n";
      std::cerr<<"\tCell angles in degrees:            "<<header.cellb[0]<<",  "<< header.cellb[1]<<",  "<<header.cellb[2]<<"\n";
      std::cerr<<"\tAxis for cols, rows, sections:     "<<header.mapc<<",  "<< header.mapr<<",  "<<header.maps<<"\n";
      std::cerr<<"\tMin, max, mean density value:      "<<header.dmin<<",  "<< header.dmax<<",  "<< header.dmean<<"\n";
      std::cerr<<"\tSpace group number:                "<<header.ispg<<"\n";
      std::cerr<<"\tNum of bytes for symmetry data:    "<<header.nsymbt<<"\n";
      std::cerr<<"\tOrigin in X,Y,Z:                   "<<header.origin[0]<<",  "<< header.origin[1]<<",  "<< header.origin[2]<<"\n";
      std::cerr<<"\tFile type:                         "<<header.map[0]<<",  "<<header.map[1]<<",  "<<header.map[2]<<",  "<<header.map[3]<<"\n";
      std::cerr<<"\tMachine stamp:                     "<<header.machst<<"\n";
      std::cerr<<"\trms deviationfrom mean density:    "<<header.rms<<"\n";
      std::cerr<<"\tNum of labels being used:          "<<header.nlabels<<"\n";
      */
      	//std::cerr<<"Getting THE DATA\n";
        //get the data
      	IteratorType  underProcessingIt( underProcessingImage , underProcessingImage->GetLargestPossibleRegion() );
      	for (underProcessingIt.GoToBegin(); !underProcessingIt.IsAtEnd(); ++underProcessingIt){
      				inputImageVector[kk]=underProcessingIt.Get();
      				underProcessingVector[kk]=underProcessingIt.Get();
      				kk++;

      	}




      // *************************************
      // MAIN CODE
      // **************************************
      // writeVectorMrcImage(parameters.outputFilename, underProcessingVector,  nx,  ny,  nz)



      // *****************
      //  END OF ACTIONS
      // *****************
      //
      			kk=0;
      			for (underProcessingIt.GoToBegin(); !underProcessingIt.IsAtEnd(); ++underProcessingIt){
      				underProcessingIt.Set(underProcessingVector[kk]);
      				kk++;
      			}
      		underProcessingImage->Update();
      	  delete [] underProcessingVector;


      			typedef itk::CastImageFilter< workingImageType, outputImageType >  CastOutFilterType;
      			CastOutFilterType::Pointer       castOutFilter       = CastOutFilterType::New();
      			castOutFilter->SetInput( underProcessingImage );
      			castOutFilter->Update();

        if (!isOutMrcFile){
           typedef itk::ImageFileWriter< outputImageType > WriterType;
           WriterType::Pointer Writer = WriterType::New();
           Writer->SetFileName(parameters.outputFilename);
           Writer->SetInput( castOutFilter->GetOutput() );
           Writer->Update();
        }  else{
           typedef itk::ImageFileWriter< outputImageType > WriterType;
           WriterType::Pointer Writer = WriterType::New();
           Writer->SetFileName(parameters.outputFilename);
           Writer->SetImageIO( itk::MRCImageIO::New() );
           Writer->SetInput( castOutFilter->GetOutput() );
           Writer->Update();
        }

}

 return 0;

}
