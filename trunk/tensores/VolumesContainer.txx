/*=========================================================================

  Program:   VolumesContainer.cxx
  Language:  C++
  Date:      7-05-2008
  Version:   1.0

  Copyright (c) 2008 Laboratoy of Image Processing, UVA. All rights reserved.
  See http://www.lpi.tel.uva.es/UsimagTool for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE. 

=========================================================================*/
#ifndef _VolumesContainer_txx
#define _VolumesContainer_txx

#include "VolumesContainer.h"


	
	template < class DataType >
	VolumesContainer<DataType>
	::VolumesContainer()
	{
		allChoices   = NULL;
		numChoices   = 0;
		baseChoices  = NULL;
		allBrowsers  = NULL;
		numBrowsers  = 0;
		baseBrowsers = NULL;
	}
	
	template < class DataType >
	VolumesContainer<DataType>
	::~VolumesContainer()
	{
		delete[] allChoices;
		delete[] baseChoices;
		delete[] allBrowsers;
		delete[] baseBrowsers;
	}
	
	template < class DataType >
	void VolumesContainer<DataType>
	::rename( unsigned int pos, const char* newname )
	{
		if( pos>=this->size() )
			return;
		(*this)[pos].nombre = std::string( newname );
		/** --------------------------------------------------------------------- */
		// We have to make sure that the name string is unique, so we need some extra logic:
		// If there is an older volume, we name it _old
		char name[400];
		int  allowed = pos;
		strcpy( name, newname );
		do{ allowed = this->uniqueString( allowed, name ); }while( allowed >= 0 );
		this->consolidateNames();
		/** --------------------------------------------------------------------- */
	}
	
	template < class DataType >
	void VolumesContainer<DataType>
	::consolidateNames()
	{
		for( unsigned int b=0; b<numBrowsers; ++b ){
			for( unsigned int k=0; k<this->size(); ++k )
				allBrowsers[b]->text( k+1+baseBrowsers[b], (*this)[k].nombre.c_str() );
			allBrowsers[b]->redraw();
		}
		for( unsigned int b=0; b<numChoices; ++b ){
			for( unsigned int k=0; k<this->size(); ++k )
				allChoices[b]->replace( k+baseChoices[b], (*this)[k].nombre.c_str() );
			allChoices[b]->redraw();
		}
	}
	
	template < class DataType >
	int VolumesContainer<DataType>
	::uniqueString( int allowed, char* unique )
	{
		char name[400];
		for(  int k=0; k<(int)(this->size()); ++k ){
			if(   k!=allowed   &&   !strcmp( (*this)[k].nombre.c_str(), unique )   ){
				sprintf( name, "%s_old", unique );
				(*this)[k].nombre = std::string( name );
				strcpy( unique, name );
				return k;
			}
		}
		return -1;
	}
	
	template < class DataType >
	void VolumesContainer<DataType>
	::push_back( const DataType& datum )
	{
		/** --------------------------------------------------------------------- */
		// We have to make sure that the name string is unique, so we need some extra logic:
		// If there is an older volume, we name it _old
		char name[400];
		int  allowed = -1;
		strcpy( name, datum.nombre.c_str() );
		do{ allowed = this->uniqueString( allowed, name ); }while( allowed >= 0 );
		this->consolidateNames();
		/** --------------------------------------------------------------------- */
		
		this->Superclass::push_back( datum );
		
		if( numBrowsers>0 ){
			for( unsigned int j=0; j<numBrowsers; ++j ){
				allBrowsers[j]->add( datum.nombre.c_str() );
				allBrowsers[j]->select( allBrowsers[j]->size() );
				allBrowsers[j]->redraw();
			}
		}
		
		if( numChoices>0 ){
			for( unsigned int j=0; j<numChoices; ++j ){
				allChoices[j]->add( datum.nombre.c_str(), 0, NULL, NULL, 0 );
				allChoices[j]->value( allChoices[j]->size()-2 );
				allChoices[j]->redraw();
			}
		}
	}
	
	template < class DataType >
	void VolumesContainer<DataType>
	::add( const DataType& datum )
	{
		this->push_back( datum );
	}
	
	template < class DataType >
	void VolumesContainer<DataType>
	::erase( unsigned int i )
	{
		if( i<this->size() ){
			(*this)[i].DestroyData();
			this->Superclass::erase( this->begin() + i );
			if( numChoices>0 ){
				for( unsigned int j=0; j<numChoices; ++j ){
					allChoices[j]->remove( i + baseChoices[j] );
					allChoices[j]->redraw();
				}
			}
			if( numBrowsers>0 ){
				for( unsigned int j=0; j<numBrowsers; ++j ){
					allBrowsers[j]->remove( i + baseBrowsers[j] + 1 );
					allBrowsers[j]->redraw();
				}
			}
		}
	}
	
	template < class DataType >
	unsigned int VolumesContainer<DataType>
	::set( unsigned int i, const DataType& datum )
	{
		if( i<this->size() ){
			*this[i] = datum;
			return i;
		}
		else{
			this->push_back( datum );
			return this->size();
		}
	}
	
	template < class DataType >
	DataType& VolumesContainer<DataType>
	::get( unsigned int i )
	{
		DataType data;
		if( i<this->size() ){
			return *this[i];
		}
		else if( this->size()>0 ){
			return *this[this->size()-1];
		}
		else{
			return data;
		}
	}
	
	template < class DataType >
	void VolumesContainer<DataType>
	::addChoice( Fl_Choice* viewer )
	{
		for( unsigned int i=0; i<this->size(); ++i )
			viewer->add( (*this)[i].nombre.c_str(), 0, NULL, NULL, 0 );
		viewer->redraw();
		
		++numChoices;
		
		unsigned int* bases = new unsigned int[numChoices];
		for( unsigned int j=0; j<numChoices-1; ++j )
			bases[j] = baseChoices[j];
		if( viewer->size() > 0 )
			bases[numChoices-1] = viewer->size()-1;
		else
			bases[numChoices-1] = 0;
		delete[] baseChoices;
		baseChoices = bases;
		
		
		Fl_Choice** viewers = new Fl_Choice*[numChoices];
		for( unsigned int j=0; j<numChoices-1; ++j )
			viewers[j] = allChoices[j];
		viewers[numChoices-1] = viewer;
		delete[] allChoices;
		allChoices = viewers;
	}
	
	template < class DataType >
	void VolumesContainer<DataType>
	::removeChoice( unsigned int j )
	{
		if( j<numChoices ){
			if( numChoices>1 ){
				--numChoices;
				Fl_Choice**   viewers = new Fl_Choice*[numChoices];
				unsigned int* bases   = new unsigned int[numChoices];
				for( unsigned int k=0; k<j; ++k ){
					viewers[k] = allChoices[k];
					bases[k]   = baseChoices[k];
				}
				for( unsigned int k=j; k<numChoices; ++k ){
					viewers[k] = allChoices[k+1];
					bases[k]   = baseChoices[k+1];
				}
				delete[] allChoices;
				delete[] baseChoices;
				allChoices  = viewers;
				baseChoices = bases;
			}
			else{
				delete[] allChoices;
				delete[] baseChoices;
				numChoices = 0;
			}
		}
	}
	
	template < class DataType >
	void VolumesContainer<DataType>
	::addBrowser( Fl_Browser* viewer )
	{
		for( unsigned int i=0; i<this->size(); ++i )
			viewer->add( (*this)[i].nombre.c_str() );
		viewer->redraw();
		
		++numBrowsers;
		
		unsigned int* bases = new unsigned int[numBrowsers];
		for( unsigned int j=0; j<numBrowsers-1; ++j )
			bases[j] = baseBrowsers[j];
		bases[numBrowsers-1] = viewer->size();
		delete[] baseBrowsers;
		baseBrowsers = bases;
		
		Fl_Browser** viewers = new Fl_Browser*[numBrowsers];
		for( unsigned int j=0; j<numBrowsers-1; ++j )
			viewers[j] = allBrowsers[j];
		viewers[numBrowsers-1] = viewer;
		delete[] allBrowsers;
		allBrowsers = viewers;
	}
	
	template < class DataType >
	void VolumesContainer<DataType>
	::removeBrowser( unsigned int j )
	{
		if( j<numBrowsers ){
			if( numBrowsers>1 ){
				--numBrowsers;
				Fl_Browser**  viewers = new Fl_Browser*[numBrowsers];
				unsigned int* bases   = new unsigned int[numBrowsers];
				for( unsigned int k=0; k<j; ++k ){
					viewers[k] = allBrowsers[k];
					bases[k]   = baseBrowsers[k];
				}
				for( unsigned int k=j; k<numBrowsers; ++k ){
					viewers[k] = allBrowsers[k+1];
					bases[k]   = baseBrowsers[k+1];
				}
				delete[] allBrowsers;
				delete[] baseBrowsers;
				allBrowsers  = viewers;
				baseBrowsers = bases;
			}
			else{
				delete[] allBrowsers;
				delete[] baseBrowsers;
				numBrowsers = 0;
			}
		}
	}
	
	template< class DataType >
	void VolumesContainer< DataType >
	::generateNewData( InputImageType::Pointer aux )
	{
		char nombre[200];
		sprintf( nombre, "New_%d", this->size() );
		
		DataType datum;
		if(   datum.generateNewData( aux, this->size(), nombre )   )
			this->push_back( datum );
	}
	
	template< class DataType >
	void VolumesContainer< DataType >
	::generateNewData( InputImageType::Pointer aux, const char* name )
	{
		DataType datum;
		if(   datum.generateNewData( aux, this->size(), name )   )
			this->push_back( datum );
	}
	
	template< class DataType >
	void VolumesContainer< DataType >
	::copyData( unsigned int j, InputImageType::Pointer aux )
	{
		if( j<0 || j>=this->size() )
			this->generateNewData( aux );
		else
			(*this)[j].copyData( aux );
	}
	
	
#endif



