/*=========================================================================

  Program:	 
  Module:    NrrdToDWIReader.txx
  Language:  C++
  Date:      $Date: 2008/06/04 $
  Version:   $ $

=========================================================================*/
#ifndef _itkNrrdToDWIReader_txx
#define _itkNrrdToDWIReader_txx

#include "itkNrrdToDWIReader.h"
#include <stdlib.h>


namespace itk
{

template <class T >
void NrrdToDWIReader<T>::GenerateData( ){
	const char* fileName = this->GetFileName();
	if( fileName==(char*)NULL )
		itkExceptionMacro( << "File name has not been assigned" );
	this->Superclass::GenerateData();
	//==========================================================================
	// Determine the number of actual DWI gradients and the number of baselines.
	DirectionVectorType     grad;
	DiffusionDirectionsType directions;
	unsigned int gradients = 0;
	unsigned int baselines = 0;
	std::vector<unsigned int> pGradients;
	std::vector<unsigned int> pBaselines;
	FILE* fid = fopen( fileName, "r" );
	if( fid==NULL )
		itkExceptionMacro( << "Cannot open " << fileName << " for reading");
	// Parse the header file...
	char* aux;
	char* aux2;
	char* test;
	char number[200];
	bool control              = true;
	unsigned int counter      = 0;
	const unsigned int LENGTH = 2000;
	char buffer[LENGTH];
	
	while(   control   ){ // Read until EOF
		if(   fgets( (char*)buffer, LENGTH, fid ) == (char*)NULL   ){
			control = false; // Check EOF
			break;
		}
		if( !strncmp(buffer,"DWMRI_gradient_",15) ){ // Is this line a gradient line?
			buffer[LENGTH-1] = '\0';
			//======================================================================================================================
			// FIRST GRADIENT COMPONENT
			// Find the first blank space:
			aux = (char*)buffer + 15*sizeof(char);
			while( *aux!='=' ){
				if( *aux=='\0'){
					fclose(fid);
					itkExceptionMacro( << "Corrupt header; check for line: " << buffer );
				}
				aux += sizeof(char);
			}
			aux += sizeof(char);
			// Find the second blank space:
			aux2 = aux;                                               // Initiallise aux2
			while( *aux2==' ' || *aux2=='\t' ){aux2+=sizeof(char);}   // Avoid multiple spaces/tabs
			while( *aux2!=' ' && *aux2!='\t' ){
				if( *aux2=='\0'){
					fclose(fid);
					itkExceptionMacro( << "Corrupt header; check for line: " << buffer );
				}
				aux2 += sizeof(char);
			}
			// Copy the number into the number buffer:
			unsigned int l=0;
			while( aux<aux2 ){
				number[l++] = *aux;
				aux += sizeof(char);
			}
			number[l] = '\0';
			// Do not substitut strtod with atof!!! atof does not perform error checking. Fix in Windows
			grad[0] = (double)::strtod(   (const char*)number,   &test   );
			if( test <= (const char*)number ){
				fclose( fid );
				itkExceptionMacro( << "Corrupt header; check for line: " << buffer );
			}
			//======================================================================================================================
			// SECOND GRADIENT COMPONENT
			// Find the third blank space:
			while( *aux2==' ' || *aux2=='\t' ){aux2+=sizeof(char);}   // Avoid multiple spaces
			while( *aux2!=' ' && *aux2!='\t' ){
				if( *aux2=='\0'){
					fclose(fid);
					itkExceptionMacro( << "Corrupt header; check for line: " << buffer );
				}
				aux2 += sizeof(char);
			}
			// Copy the number into the number buffer:
			l=0;
			while( aux<aux2 ){
				number[l++] = *aux;
				aux += sizeof(char);
			}
			number[l] = '\0';
			// Do not substitut strtod with atof!!! atof does not perform error checking. Fix in Windows
			grad[1] = (double)::strtod(   (const char*)number,   &test   );
			if( test <= (const char*)number ){
				fclose( fid );
				itkExceptionMacro( << "Corrupt header; check for line: " << buffer );
			}
			//======================================================================================================================
			// THIRD GRADIENT COMPONENT
			while( *aux2==' ' || *aux2=='\t' ){aux2+=sizeof(char);}   // Avoid multiple spaces
			l=0;
			while( *aux2!=' ' && *aux2!='\t' && *aux2!='\n' && *aux2!='\r' && *aux2!='\0' ){
				number[l++] = *aux2;
				aux2 += sizeof(char);
			}
			number[l] = '\0';
			// Extract the last number:
			// Do not substitut strtod with atof!!! atof does not perform error checking. Fix in Windows
			grad[2] = (double)::strtod(   number,   &test   );
			if( test <= (const char*)number ){
				fclose( fid );
				itkExceptionMacro( << "Corrupt header; check for line: " << buffer );
			}
			//======================================================================================================================
			// CHECK FOR GRADIENT/BASELINE
			// Extract the norm:
			double norm = ::sqrt( grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2] );
			if( norm>1e-1 ){ //This is an actual gradient direction!
				++gradients;
				grad[0] /= norm;   grad[1] /= norm;  grad[2] /= norm;
				directions.push_back( grad );
				pGradients.push_back( counter++ );
			}
			else{ //This is a baseline
				++baselines;
				pBaselines.push_back( counter++ );
			}
			//======================================================================================================================
		}
		// Not a gradient line, ignore!
	}
	// The file has beend correctly parsed!
	fclose( fid );
	// Set the proper parameters:
	IndicatorType pG( gradients );
	IndicatorType pB( baselines );
	for( unsigned int p=0; p<gradients; ++p )
		pG[p] = pGradients[p];
	for( unsigned int p=0; p<baselines; ++p )
		pB[p] = pBaselines[p];
		
	if( baselines<1 )
		itkExceptionMacro( << "No T2 baselines were found in " << this->GetFileName() );
	if( gradients<6 )
		itkExceptionMacro( << "This DWI volume contains less than 6 gradient directions!!!" );
	if( baselines+gradients != this->GetOutput()->GetVectorLength() )
		itkExceptionMacro( << "The total number of channels in " << this->GetFileName() << " is different from the number of gradients and baselines described in its header" );
	this->GetOutput()->SetNumBaselines( baselines );
	this->GetOutput()->SetNumDWImages( gradients );
	this->GetOutput()->SetNumImages( baselines + gradients );
	this->GetOutput()->SetDiffusionDirections( directions );
	this->GetOutput()->SetIndexesDWImages( pG );
	this->GetOutput()->SetIndexesBaselines( pB );
	//==========================================================================
	
	//==========================================================================
	// Parse the input file to obtain the b-value:
	double bvalue = itk::NumericTraits<double>::Zero;
	bool   success = false;
	
	fid = (FILE*)NULL;
	fid = fopen( fileName, "r" );
	if( fid==NULL )
		itkExceptionMacro( << "Cannot open file: " << fileName << " for reading" );
	control = true;
	
	while(   control   ){ // Read until EOF
		if(   fgets( (char*)buffer, LENGTH, fid )   ==   (char*)NULL   ){
			control = false; // Check EOF
			break;
		}
		if( !strncmp(buffer,"DWMRI_b-value:=",15) ){
			// Do not substitut strtod with atof!!! atof does not perform error checking. Fix in Windows
			bvalue = (double)::strtod(   (char*)buffer + 15*sizeof(char),   &test   );
			if( test <= (char*)buffer + 15*sizeof(char) ){
				fclose( fid );
				itkExceptionMacro( << "Corrupt header; check for line: " << buffer );
			}
			success = true;
			control = false;
		}
	}
	fclose( fid );
	if( !success )
		itkExceptionMacro( << "Cannot determine the value of b in " << fileName );
	BValuesType bvalues( gradients + baselines );
	bvalues.Fill( bvalue );
	for( unsigned int k=0; k<baselines; ++k )
		bvalues[ pB[k] ] = 0;
	this->GetOutput()->SetBValues( bvalues );
	this->GetOutput()->SetB_Value( bvalue );
	this->GetOutput()->SetUniqueB( true ); // We only allow single values of B for Nrrd.
	//==========================================================================
	
	return;
}  


} // end namespace itk

#endif
