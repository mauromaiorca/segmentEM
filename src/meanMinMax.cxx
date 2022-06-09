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


//ITK
/*#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkCastImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryThinningImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
*/

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
    std::cerr<<"Usage: " << argv[0] << " locresMapFileName.mrc maskFileName.mrc \n";
    std::cerr<<"\n";
    exit(1);
}


// *************************************
//
// retrieveInputParameters
// *************************************
typedef struct inputParametersType {
    char * locresMapFileName;
    char * maskFileName;
}inputParametersType;


void retrieveInputParameters(inputParametersType * parameters, int argc, char** argv){
    parameters->locresMapFileName= NULL;
    parameters->maskFileName= NULL;

    if (argc >= 2){
     parameters->locresMapFileName=argv[1];
    }

    if (argc >= 3){
     parameters->maskFileName=argv[2];
    }

    
    if (!parameters->locresMapFileName && !parameters->maskFileName){
      usage(argv);
    }

    
  }

// ******************************
// ******************************
//   MRC I/O
enum  DataType {UChar, SChar, Short, Int, UShort, UInt, Float,  RFloat, Double, Bool};


//https://www.sciencedirect.com/science/article/pii/S104784771500074X
struct MRCHeader
{             // file header for MRC data
    int nx;              //  0   0       image size
    int ny;              //  1   4
    int nz;              //  2   8
    int mode;            //  3           0=char,1=short,2=float
    int nxstart;         //  4           unit cell offset
    int nystart;         //  5
    int nzstart;         //  6
    int mx;              //  7           Number of intervals along X of the “unit cell” in voxels
    int my;              //  8		 Number of intervals along Y of the “unit cell” in voxels
    int mz;              //  9           Number of intervals along Z of the “unit cell” in voxels
    float a;             // 10   40      cell dimensions in A
    float b;             // 11
    float c;             // 12
    float alpha;         // 13           cell angles in degrees
    float beta;          // 14
    float gamma;         // 15
    int mapc;            // 16           column axis
    int mapr;            // 17           row axis
    int maps;            // 18           section axis
    float dmin;          // 19           minimum density value
    float dmax;          // 20   80      maximum density value
    float dmean;         // 21           average density value
    int ispg;            // 22           space group number
    int nsymbt;          // 23           bytes used for sym. ops. table
    float extra[25];     // 24           user-defined info
    float xOrigin;       // 49           phase origin in pixels
    float yOrigin;       // 50
    float zOrigin;       // 51
    char map[4];         // 52       identifier for map file ("MAP ")
    char machst[4];      // 53           machine stamp
    float drms;          // 54       RMS deviation
    int nlabl;           // 55           number of labels used
    char labels[800];    // 56-255       10 80-character labels
    MRCHeader(){
        nx=0;
        ny=0;
        nz=0;
        mode=2;
        nxstart=0;
        nystart=0;
        nzstart=0;
        mx=0;
        my=0;
        mz=0;
        a=0;
        b=0;
        c=0;
        alpha=90;
        beta=90;
        gamma=90;
        mapc=1;
        mapr=2;
        maps=3;
        dmin=0;
        dmax=0;
        dmean=0;
        ispg=0;
        nsymbt=0;
        xOrigin=0;
        yOrigin=0;
        zOrigin=0;
        drms=0;
        nlabl=1;
        }
    MRCHeader(int __nx, int __ny, int __nz):nx(__nx),ny(__ny),nz(__nz){
        mode=2;
        nxstart=0;
        nystart=0;
        nzstart=0;
        mx=__nx;
        my=__ny;
        mz=__nz;
        a=__nx;
        b=__ny;
        c=__nz;
        alpha=90;
        beta=90;
        gamma=90;
        mapc=1;
        mapr=2;
        maps=3;
        dmin=0;
        dmax=0;
        dmean=0;
        ispg=0;
        nsymbt=0;
        xOrigin=0;
        yOrigin=0;
        zOrigin=0;
        drms=0;
        nlabl=1;
        }

} ;
void readHeaderMrc(const char * filenameMRCinput, MRCHeader & header){

  //read input header
  FILE *input_fp = NULL;
	input_fp=fopen(filenameMRCinput,"r");
	if(input_fp==NULL) return;
	rewind(input_fp);
	if(fread(&header,1,1024,input_fp)<1024) return;
  fclose(input_fp);

}
template<typename T>
int readMrcImage(const char * filenameMRC, T * I, MRCHeader & header){
  //std::cerr<<"MRC read image\n";

  //constant
  const int headerSizeMRC = 1024;
  const unsigned int _MAX_UNSIGNED_SHORT_ = 65535;

  //header
  FILE *m_fp = NULL;
  m_fp=fopen(filenameMRC,"r+");
  if(m_fp==NULL) return 1;

  //read file header
  rewind(m_fp);
  if(fread(&header,1,headerSizeMRC,m_fp)<headerSizeMRC)
	return 1;


  //READ THE FILE
  int         i;
  if ( ( abs( header.mode ) > _MAX_UNSIGNED_SHORT_ ) || ( abs(header.nx) > _MAX_UNSIGNED_SHORT_ ) ){

        std::cerr<<"Warning: Swapping header byte order for 4-byte types\n";
        int     extent = headerSizeMRC - 800; // exclude labels from swapping
        char _tmp;
	for ( i=0; i<extent; i+=4 ){
	    //swaps bytes
	    char* v = (char *) (&(header)+i);
	    for ( int iTmp=0; iTmp<2; iTmp++ )
	    {
		_tmp = v[iTmp];
		v[iTmp] = v[3-iTmp];
		v[3-iTmp] = _tmp;
	    }
	}
    }



    const unsigned long int nx =(int) header.nx;
    const unsigned long int ny =(int) header.ny;
    const unsigned long int nz =(int) header.nz;
    const unsigned long int nxyz = nx * ny * nz;

    DataType datatype;
    size_t datatypesize;
    unsigned long int blockMemory=268435456; //256 Mb
    void * buf;


    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
	//std::cerr<<"File mode 0: UChar\n";
	datatype=UChar;
	datatypesize = 1;
        buf = new unsigned char [blockMemory/datatypesize];
    }else if (header.mode==1){ //image : 16-bit halfwords
	//std::cerr<<"File mode 1: Short\n";
	datatype=Short;
	datatypesize = 2;
	buf = new short [blockMemory/datatypesize];
    }else if (header.mode==2){ //image : 32-bit reals
	//std::cerr<<"File mode 2: Float\n";
	datatype=Float;
        datatypesize = 4;
	buf = new float [blockMemory/datatypesize];
    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
	//std::cerr<<"File mode 5: Char\n";
	datatype=SChar;
	datatypesize = 1;
	buf = new char [blockMemory/datatypesize];
    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
	//std::cerr<<"File mode 6: UShort\n";
	datatype=UShort;
	datatypesize = 2;
	buf = new unsigned short [blockMemory/datatypesize];
    }else{
	//std::cerr<<"Unsupported type\n";
	datatype=UChar;
	datatypesize = 1;
	buf = new unsigned char [blockMemory/datatypesize];
    }

    unsigned long int blockSize=(unsigned long int)blockMemory/(double)datatypesize; 

    //std::cerr<<"read header\n";
    int error_fseek = fseek( m_fp, headerSizeMRC, SEEK_SET );
    if (error_fseek != 0)
           return -1;
    
    //std::cerr<<"read data\n";
    for ( unsigned long int tt=0; tt<nxyz; tt+=blockSize ){
	unsigned long int endRead=tt+blockSize;
	if ( endRead >= nxyz ){
		endRead=nxyz;
	}
	unsigned long int toRead=endRead-tt;

	// IT MAY BE SLIGHTLY MORE EFFICIENT TO PUT IT IN DIRECTLY
	// BUT NEED TO MANAGE POLYMORPHISM in EXECUTION
	/*if (sizeof(T)==datatypesize ){
		std::cerr<<"same type t="<<tt<<"->"<<endRead<<"   (max="<<nxyz<<")\n";
		size_t result = fread( &I[tt], datatypesize, toRead, m_fp );
		if (result < 0 ){
			std::cerr<<"Error: wrong read\n";
        	        return -1;
		}
	}*/	
		size_t result = fread( buf, datatypesize, toRead, m_fp );
		if (result < 0 ){
			std::cerr<<"Error: wrong read\n";
        	        return -1;
		}
	
		//std::cerr<<" block ="<<tt<<"->"<<endRead<<"   (max="<<nxyz<<")\n";

		    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
			for (unsigned long int ii=0;ii<toRead;ii++)
				I[tt+ii]=(T)(((unsigned char*) buf)[ii]);
		    }else if (header.mode==1){ //image : 16-bit halfwords
			for (unsigned long int ii=0;ii<toRead;ii++)
				I[tt+ii]=(T)(((short*) buf)[ii]);
		    }else if (header.mode==2){ //image : 32-bit reals
			for (unsigned long int ii=0;ii<toRead;ii++)
				I[tt+ii]=(T)(((float*) buf)[ii]);
		    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
			for (unsigned long int ii=0;ii<toRead;ii++)
				I[tt+ii]=(T)(((char*) buf)[ii]);
		    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
			for (unsigned long int ii=0;ii<toRead;ii++)			
				I[tt+ii]=(T)(((unsigned short*) buf)[ii]);
		    }else{
			//std::cerr<<"Unsupported type\n";
			for (unsigned long int ii=0;ii<toRead;ii++)
				I[tt+ii]=(T)(((unsigned char*) buf)[ii]);
		    }
	}

    
    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
        delete [] (unsigned char *) buf;
    }else if (header.mode==1){ //image : 16-bit halfwords
	delete [] (short *) buf;
    }else if (header.mode==2){ //image : 32-bit reals
	delete [](float *) buf;
    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
	delete [](char *) buf;
    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
	delete [](unsigned short *) buf;
    }else{
	delete [] (unsigned char *) buf;
    }

  //close the file
  fclose(m_fp);

}
template<typename T>
T * writeMrcImage(const char * filenameMRC, T * I, MRCHeader & header){
  //check the header

    //detect if big endilan or little endian
    if ( htonl(47) == 47 ) { // Big endian
       header.machst[0] = header.machst[1] = 17;
    } else {// Little endian  
	header.machst[0] = 68;
	header.machst[1] = 65;
    }


    unsigned long int nx = header.nx;
    unsigned long int ny = header.ny;
    unsigned long int nz = header.nz;
    if (nx<1) nx=1;
    if (ny<1) ny=1;
    if (nz<1) nz=1;
    unsigned long int nxyz = nx*ny*nz;


    DataType datatype;
    size_t datatypesize;
    unsigned long int blockMemory=268435456; //256 Mb
    void * buf = NULL;

    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
	//std::cerr<<"File mode 0: UChar\n";
	datatype=UChar;
	datatypesize = 1;
        buf = new unsigned char [blockMemory/datatypesize];
    }else if (header.mode==1){ //image : 16-bit halfwords
	//std::cerr<<"File mode 1: Short\n";
	datatype=Short;
	datatypesize = 2;
	buf = new short [blockMemory/datatypesize];
    }else if (header.mode==2){ //image : 32-bit reals
	//std::cerr<<"File mode 2: Float\n";
	datatype=Float;
        datatypesize = 4;
	buf = new float [blockMemory/datatypesize];
    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
	//std::cerr<<"File mode 5: Char\n";
	datatype=SChar;
	datatypesize = 1;
	buf = new char [blockMemory/datatypesize];
    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
	//std::cerr<<"File mode 6: UShort\n";
	datatype=UShort;
	datatypesize = 2;
	buf = new unsigned short [blockMemory/datatypesize];
    }else{
	//std::cerr<<"Unsupported type\n";
	datatype=UChar;
	datatypesize = 1;
	buf = new unsigned char [blockMemory/datatypesize];
    }

   if ( !buf ){
	std::cerr<<"failed to allocate memory ("<< (int)blockMemory/1024 <<"Mb) for the buffer\n";
	return NULL;
   }

    //compute and write file header with amended statistical values
    double dmean = 0;
    double dRMS = 0;
    double dmin = (double)I[0];
    double dmax = (double)I[0];
    for (unsigned long int ii=0;ii<nxyz;ii++){
	dmean+=(double)I[ii];
    }
    dmean/=nxyz;
    for (unsigned long int ii=0;ii<nxyz;ii++){
	double tmpI=(double)I[ii];
	dRMS+=pow(tmpI-dmean,2.0);
	if (dmin>tmpI) dmin=tmpI;
	if (dmax<tmpI) dmax=tmpI;
    }
    dRMS/=(nxyz);
	
    header.dmean = (float)dmean;
    header.drms   = (float)dRMS;
    header.dmin  = (float)dmin;
    header.dmax  = (float)dmax;

     FILE *m_fp = NULL;
     m_fp=fopen(filenameMRC,"w");
     if(m_fp==NULL) return NULL;

     //write file header
     //rewind(m_fp);
     fseek(m_fp, 0, SEEK_SET);
     fwrite(&header,1,1024,m_fp);

     unsigned long int blockSize=(unsigned long int)blockMemory/(double)datatypesize;


    for ( unsigned long int tt=0; tt<nxyz; tt+=blockSize ){
	unsigned long int endWrite=tt+blockSize;
	if ( endWrite >= nxyz ){
		endWrite=nxyz;
	}
	unsigned long int toWrite=endWrite-tt;

		//std::cerr<<" block ="<<tt<<"->"<<endWrite<<"   (max="<<nxyz<<")\n";

		    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
			for (unsigned long int ii=0;ii<toWrite;ii++)
				((unsigned char*)buf)[ii]=(unsigned char)I[tt+ii];
		    }else if (header.mode==1){ //image : 16-bit halfwords
			for (unsigned long int ii=0;ii<toWrite;ii++)
				((short*)buf)[ii]=(unsigned long int)I[tt+ii];
		    }else if (header.mode==2){ //image : 32-bit reals
			for (unsigned long int ii=0;ii<toWrite;ii++)
				((float*)buf)[ii]=(float)I[tt+ii];
		    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
			for (unsigned long int ii=0;ii<toWrite;ii++)
				((char*)buf)[ii]=(char)I[tt+ii];
		    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
			for (unsigned long int ii=0;ii<toWrite;ii++)			
				((unsigned short*)buf)[ii]=(unsigned short)I[tt+ii];
		    }else{
			//std::cerr<<"Unsupported type\n";
			for (unsigned long int ii=0;ii<toWrite;ii++)
				((unsigned char*)buf)[ii]=(unsigned char)I[tt+ii];
		    }
		size_t result = fwrite( buf, datatypesize, toWrite, m_fp );
		if (result < 0 ){
			std::cerr<<"Error: wrong read\n";
        	        return NULL;
		}
	}



    if (header.mode==0){ //image : signed 8-bit bytes range -128 to 127
        delete [] (unsigned char *) buf;
    }else if (header.mode==1){ //image : 16-bit halfwords
	delete [] (short *) buf;
    }else if (header.mode==2){ //image : 32-bit reals
	delete [](float *) buf;
    }else if (header.mode==5){ //image : unsigned 8-bit range 0 to 255
	delete [](char *) buf;
    }else if (header.mode==6){ //image : unsigned 16-bit range 0 to 65535
	delete [](unsigned short *) buf;
    }else{
	delete [] (unsigned char *) buf;
    }
  //close the file
  fclose(m_fp);
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
           readHeaderMrc(parameters.locresMapFileName, imageHeader);
           unsigned long int nx=imageHeader.nx, ny=imageHeader.ny, nz=imageHeader.nz;
           unsigned long int nxyz=nx*ny*nz;
           WorkingPixelType * I = new WorkingPixelType [nxyz];
           WorkingPixelType * maskI = new WorkingPixelType [nxyz];
           readMrcImage(parameters.locresMapFileName, I, imageHeader);
           readMrcImage(parameters.maskFileName, maskI, imageHeader);
           long int firstIdx = 0;
           double mean=0;
           double min=0;
           double max=0;
           double counter=1;
           std::vector<double> values;
           for (unsigned long int ii=0, found = 0; ii<nxyz && found==0;ii++){
             if (maskI[ii]>0.2){
               found = 1;
               firstIdx=ii;
               mean=I[ii];
               min=I[ii];
               max=I[ii];
               values.push_back(I[ii]);
             }
           }
           for (unsigned long int ii=firstIdx+1; ii<nxyz; ii++){
             if (maskI[ii]>0.2){
               firstIdx=ii;
               mean+=I[ii];
               if (min>I[ii]) min=I[ii];
               if (max<I[ii]) max=I[ii];
               counter++;
               values.push_back(I[ii]);
             }
           }
           mean/=counter;
           
        double sd = 0;
        for (unsigned long int ii=0; ii<nxyz;ii++){
                    if (maskI[ii]>0.2){
                        sd+=pow(I[ii]-mean,2.0);
                    }
        }
        if ( counter - 1.0 > 0 ){
            sd = sd / (counter - 1.0 );
        }
        std::sort (values.begin(), values.end());


//        std::cout<<min<<","<< mean <<","<<max<<"\n";
        int firstQuartileIdx = values.size()/4.0;
        int lastQuartileIdx = values.size()*3.0/4.0;
        double minFinal = 0;
        double maxFinal = 0;
        if (values.size()>0){
            if (values.size()>4){
                minFinal=values[2];
                maxFinal=values[values.size()-3];
            }else{
                minFinal=values[0];
                maxFinal=values[values.size()-1];
            }

        }
        maxFinal = values[(unsigned long int)values.size()*(3.70/4.0)];
        std::cout<<minFinal<<","<<values[firstQuartileIdx]<<","<<mean<<","<<values[lastQuartileIdx]<<","<<maxFinal<<"\n";

        delete [] I;
        delete [] maskI;

}






