// Evolutions_01.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "evoEngine.h"


#define kFilePath "dataFile"
#define kLineMaxLength 1024
#define kNumFields 9

/*
  <w lem="find" pos="vvb" xml:id="sha-1h410100201" eos="0" ord="12" reg="find" spe="Find" tok="Find">Find</w> 
  <w lem="we" pos="pns12" xml:id="sha-1h410100202" eos="0" ord="13" reg="we" spe="we" tok="we">we</w> 
  <w lem="a" pos="dt" xml:id="sha-1h410100203" eos="0" ord="14" reg="a" spe="a" tok="a">a</w> 
  <w lem="time" pos="n1" xml:id="sha-1h410100204" eos="0" ord="15" reg="time" spe="time" tok="time">time</w> 
  <w lem="for" pos="p-acp" xml:id="sha-1h410100205" eos="0" ord="16" reg="for" spe="for" tok="for">for</w> 
  <w lem="fright" pos="j-vvn" xml:id="sha-1h410100206" eos="0" ord="17" reg="frighted" spe="frighted" tok="frighted">frighted</w> 
  <w lem="peace" pos="n1" xml:id="sha-1h410100207" eos="0" ord="18" reg="peace" spe="peace" tok="peace">peace</w> 
  <w lem="to" pos="pc-acp" xml:id="sha-1h410100208" eos="0" ord="19" reg="to" spe="to" tok="to">to</w> 
  <w lem="pant" pos="vvi" xml:id="sha-1h410100209" eos="0" ord="20" reg="pant" spe="pant" tok="pant">pant</w> 
  <w lem="," pos="," xml:id="sha-1h410100210" eos="0" ord="21" reg="," spe="," tok=",">,</w> 
*/

// lem	text
// pos	text
// id	text
// ----
// eos	integer
// ord	integer
// ----
// reg	text
// spe	text
// tok	text
// word	text
//---------------------------------------------------

bool Configure() {
	DataPointDescriptor descriptor = {
		kNumFields,
		{ eShortText, eShortText, eShortText, eInteger, eInteger, eShortText, eShortText, eShortText, eShortText }
	};
	return InitEvoEngine( &descriptor );
}

//---------------------------------------------------

void ConvertAndSave( char * dataLine ) {
	//<w lem="," pos="," xml:id="sha-1h410100210" eos="0" ord="21" reg="," spe="," tok=",">,</w>
	char phase = 0;
	char * op = dataLine;
	char * end = null;	
	char * val = null;
	DataPoint * object = new DataPoint;

	object->fields = new FieldContent[ kNumFields ];
	object->nearest = null;
	object->next = null;
	object->prev = null;

	char ** tags = { "lem", "pos", "xml:id", "eos", "ord", "reg", "spe", "tok" };

	for (int i = 0; i < 8; i++ ) {
		op = strstr( dataLine, tags[i] );
		op = strchr( op, '"' );
		end = strchr( val = op+1, '"' );
		*end = 0;

		switch (i) {
			case 0 :
			case 1 :
			case 2 :
			case 5 :
			case 6 :
			case 7 : {
				// text
				object->fields[i].integer = IntFromText( val );
			} break;
			
			case 3 :
			case 4 : {
				// integer
				object->fields[i].integer = atol( val );
			} break;
		} // switch

		end = '"';
	} // for

	op = strchr( dataLine, '>' );
	if ( op ) {
		end = strrchr( val = op+1, '<' );
		if ( end ) {
			*end = 0;
			object->fields[8].integer = IntFromText( val );
		}
		
	}
}

//---------------------------------------------------

bool ReadData( const char * filePath ) {
	// Open the file 
	FILE * file = fopen( filePath, "r" );

	if ( file == null ) return false;

	char line[ kLineMaxLength ];
	char * str;
	
	do {
		// Read line
		str = fgets( line, kLineMaxLength, file );
		if ( str != null ) {
			// Check if this the line we need
			str = strchr( str, '<' );
			if ( str != null && str[1] == 'w' ) {
				// Save it to engine
				ConvertAndSave( str );
			}
		}		
	} while ( str != null );
}


int _tmain( int argc, _TCHAR* argv[] )
{
	if ( Configure() ) {
		// Read content of the file to application
		ReadData( kFilePath );
	}

	return 0;
}

