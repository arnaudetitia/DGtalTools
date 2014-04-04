#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <string>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/kernel/sets/SetPredicate.h"

#include "DGtal/topology/KhalimskySpaceND.h"

#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/kernel/sets/SetPredicate.h"

#include "DGtal/images/imagesSetsUtils/SimpleThresholdForegroundPredicate.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/ImageFromSet.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/ImageContainerBySTLVector.h"

#include "DGtal/io/colormaps/GrayscaleColorMap.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"

#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/geometry/volumes/distance/VoronoiMap.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/geometry/helpers/ContourHelper.h"

#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include "DGtal/io/boards/Board2D.h"

#include "DGtal/topology/helpers/Surfaces.h"



using namespace std;
using namespace DGtal;
typedef GrayscaleColorMap<unsigned char> Gray;
typedef ImageSelector<Z2i::Domain, unsigned int>::Type Image;
typedef IntervalThresholder<Image::Value> Binarizer; 
typedef SimpleThresholdForegroundPredicate<Image> PointPredicate;
typedef DistanceTransformation<Z2i::Space, PointPredicate, Z2i::L2Metric> DTL2;
typedef HueShadeColorMap<long int, 2> HueTwice;
void readData(char *filename,vector<int>&x ,vector<int>&y ,vector<int>&noiseLevel){
	ifstream data(filename);
	int entierPoubelle;
	char buffer[150];
	stringstream ss;
	data.getline(buffer,100);
	while( !(data.eof()) ){
		if (buffer[0] != '#'){
			ss << buffer;
			for (int i = 0 ; i<5 ; ++i){
				ss.getline(buffer,100,' ');
				if(i==1){noiseLevel.push_back(atoi(buffer));}
				if(i==3){x.push_back(atoi(buffer));}
				if(i==4){y.push_back(atoi(buffer));}
				else {entierPoubelle = atoi(buffer);}  
			}
			ss << "";
			ss.clear();
		}
		data.getline(buffer,100); 
	} 
	data.close();
} 

void getContours(vector< vector< Z2i::Point >  >&  contours,Image &image){
	
	Binarizer b(100,200);
	PointFunctorPredicate<Image,Binarizer> predicate(image, b);

	Z2i::KSpace ks;
	if (!ks.init( image.domain().lowerBound(), image.domain().upperBound(), true )) cout << "pb" << endl;
        
	SurfelAdjacency<2> sAdj( true );
	Surfaces<Z2i::KSpace>::extractAllPointContours4C( contours,ks, predicate, sAdj ); 
} 

void constructImage(vector< Z2i::Point >& contour,Image& image){
	for ( Image::Iterator it = image.begin(), itend = image.end();it != itend; ++it)
    		(*it)=128;

	for (int i=0;i<contour.size();i++){
		image.setValue(contour[i],0);
	}
	
}


int main(){
	char *filename ="../../../ib1.txt";
	vector<int>x ;
	vector<int>y ;
	vector<int>noiseLevel;
	readData(filename,x,y,noiseLevel);// récupération des données de contours 
	x.pop_back();
	y.pop_back();
	noiseLevel.pop_back();
	//boite englobante
	int valMaxX = x[0];
	int valMaxY = y[0];
	for (int i=1;i<x.size();i++){
		if (x[i] > valMaxX){
			valMaxX = x[i];
		} 
		if (y[i] > valMaxY){
			valMaxY = y[i];
		} 
	}

	//création du domaine 
	Z2i::Point lower(0,0);
	Z2i::Point upper(valMaxX+20,valMaxY+20);
	
	//Création de l'image
	
  	Image image ( Z2i::Domain(lower,upper));

	for ( Image::Iterator it = image.begin(), itend = image.end();it != itend; ++it)
    		(*it)=128;

	for (int i=0;i<x.size();i++){
		for (int a = max(0,x[i]-(noiseLevel[i]-1) ) ; a <= min(x[i]+(noiseLevel[i]-1),valMaxX); a++ ){
			for (int b = max(0,y[i]-(noiseLevel[i]-1) ) ; b <= min(y[i]+(noiseLevel[i]-1),valMaxY); b++ ){
				image.setValue( Z2i::Point(a,valMaxY-b+20),0);
			}
		}
	}
	// affichage de l'image dans un svg	
	

	Board2D board;
	Display2DFactory::drawImage<Gray>(board, image, (unsigned int)0, (unsigned int)129);
	board.saveSVG("../../../monResultat.svg");
	
	//récupération des contours intérieurs et extérieurs
	
	vector< vector< Z2i::Point >  >  contours;
	
	getContours(contours,image);
	Image image1 ( Z2i::Domain(lower,upper));
	Image image2 ( Z2i::Domain(lower,upper));
	constructImage(contours[0],image1);
	constructImage(contours[1],image2);
	
	Board2D board1,board2;
	Display2DFactory::drawImage<Gray>(board1, image1, (unsigned int)0, (unsigned int)129);
	Display2DFactory::drawImage<Gray>(board2, image2, (unsigned int)0, (unsigned int)129);
	board1.saveSVG("../../../contour1.svg");
	board2.saveSVG("../../../contour2.svg");

	//Calcul de la transformée en distance 
	

	PointPredicate predicate1(image1,0);
	PointPredicate predicate2(image2,0);
	
	DTL2 dt1(image1.domain(), predicate1, Z2i::l2Metric);
	DTL2 dt2(image2.domain(), predicate2, Z2i::l2Metric);

	DTL2::Value maxv1=0;
  	for ( DTL2::ConstRange::ConstIterator it = dt1.constRange().begin(), itend = dt1.constRange().end();it != itend; ++it)
    		if ( (*it) > maxv1)  maxv1 = (*it);

	DTL2::Value maxv2=0;
  	for ( DTL2::ConstRange::ConstIterator it = dt2.constRange().begin(), itend = dt2.constRange().end();it != itend; ++it)
    		if ( (*it) > maxv2)  maxv2 = (*it);
	
	Board2D boardDT1,boardDT2;
	boardDT1.clear();
	boardDT2.clear();
  	Display2DFactory::drawImage<HueTwice>(boardDT1, dt1, 0.0, maxv1 + 1);
	Display2DFactory::drawImage<HueTwice>(boardDT2, dt2, 0.0, maxv2 + 1);
	boardDT1.saveSVG ( "../../../DT_contour1.svg" );
	boardDT2.saveSVG ( "../../../DT_contour2.svg" );
 	
	return 0; 
}
