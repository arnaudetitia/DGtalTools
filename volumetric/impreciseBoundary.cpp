#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <string>
#include <math.h>

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
typedef DistanceTransformation<Z2i::Space, PointPredicate, Z2i::L2Metric> VML2;
typedef HueShadeColorMap<long int, 2> HueTwice;
typedef NotPointPredicate<Z2i::DigitalSet> NotPredicate;
typedef ExactPredicateLpSeparableMetric<Z2i::Space, 2> L2Metric;
typedef VoronoiMap<Z2i::Space, NotPredicate, L2Metric > Voronoi2D;

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

void constructImageBothContour(vector< Z2i::Point >& contour1,vector< Z2i::Point >& contour2,Image& image){
	for ( Image::Iterator it = image.begin(), itend = image.end();it != itend; ++it)
    		(*it)=128;

	for (int i=0;i<contour1.size();i++){
		image.setValue(contour1[i],0);
	}
	for (int i=0;i<contour2.size();i++){
		image.setValue(contour2[i],0);
	}
	
	
}

void computeDT(Image& image,char *outputfile){
	PointPredicate predicate(image,0);
	DTL2 dt(image.domain(), predicate, Z2i::l2Metric);

	DTL2::Value maxv=0;
  	for ( DTL2::ConstRange::ConstIterator it = dt.constRange().begin(), itend = dt.constRange().end();it != itend; ++it)
    		if ( (*it) > maxv)  maxv = (*it);

	Board2D board;
	board.clear();
  	Display2DFactory::drawImage<HueTwice>(board, dt, 0.0, maxv + 1); 

	board.saveSVG ( outputfile );
}

Image imageDT(Image& image,int maxX, int maxY){
	PointPredicate predicate(image,0);
	Z2i::Point lower(0,0);
	Z2i::Point upper(maxX,maxY);
	Image result(Z2i::Domain(lower,upper));
	DTL2 dt(image.domain(), predicate, Z2i::l2Metric);

	for (int i=0;i<maxX;i++){
		for (int j=0;j<maxY;j++){
			Z2i::Point p(i,j);
			result.setValue(p,dt(p));
		}
	}
	return result;
} 

void computeVoronoiMap(Image& image,vector< Z2i::Point >& contour,char *outputfile,char* outputfileCells){
	
	Z2i::DigitalSet set(image.domain());

	for (int i=0;i<contour.size();i++){
		set.insertNew(contour[i]);
	} 

	NotPredicate notSetPred(set);
	
	L2Metric l2;
	Voronoi2D voronoimap(image.domain(),notSetPred,l2);
	
	Board2D board;
	board.clear();
	board << image.domain() << set;
	for(Voronoi2D::Domain::ConstIterator it = voronoimap.domain().begin(),
      	itend = voronoimap.domain().end(); it != itend; ++it)
  	{
    		Voronoi2D::Value site = voronoimap( *it );   //closest site to (*it)
    		if (site != (*it))
      			Display2DFactory::draw( board,   site - (*it), (*it)); //Draw an arrow
  	}
	board.saveSVG(outputfile);

	board.clear();
  	for(Voronoi2D::Domain::ConstIterator it = voronoimap.domain().begin(),
      	itend = voronoimap.domain().end(); it != itend; ++it)
  	{
   		Voronoi2D::Value site = voronoimap( *it );   //closest site to (*it)
   		 unsigned char c = (site[1]*13 + site[0] * 7) % 256; //basic hashfunction
    		board << CustomStyle( (*it).className(), new CustomColors(Color(c,c,c),Color(c,c,c)))<< (*it);
  	}
	board.saveSVG(outputfileCells);
} 

void computeVoronoiMap(Image& image,vector< vector< Z2i::Point > >& contours,char *outputfile,char* outputfileCells){
	
	Z2i::DigitalSet set(image.domain());

	for (int i=0;i<contours.size();i++){
		for (int j=0;j<contours[i].size();j++){
			set.insertNew(contours[i][j]);
		}
	} 

	NotPredicate notSetPred(set);
	
	L2Metric l2;
	Voronoi2D voronoimap(image.domain(),notSetPred,l2);
	
	Board2D board;
	board.clear();
	board << image.domain() << set;
	for(Voronoi2D::Domain::ConstIterator it = voronoimap.domain().begin(),
      	itend = voronoimap.domain().end(); it != itend; ++it)
  	{
    		Voronoi2D::Value site = voronoimap( *it );   //closest site to (*it)
    		if (site != (*it))
      			Display2DFactory::draw( board,   site - (*it), (*it)); //Draw an arrow
  	}
	board.saveSVG(outputfile);

	board.clear();
  	for(Voronoi2D::Domain::ConstIterator it = voronoimap.domain().begin(),
      	itend = voronoimap.domain().end(); it != itend; ++it)
  	{
   		Voronoi2D::Value site = voronoimap( *it );   //closest site to (*it)
   		 unsigned char c = (site[1]*13 + site[0] * 7) % 256; //basic hashfunction
    		board << CustomStyle( (*it).className(), new CustomColors(Color(c,c,c),Color(c,c,c)))<< (*it);
  	}
	board.saveSVG(outputfileCells);
} 


int main(int argc,char **argv){
	if (argc != 2){
		cout << "Usage : ./impreciseBounadary <filename> " << endl;
		exit(1);
	}	 
	char *filename = argv[1];
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
	Image imageCB ( Z2i::Domain(lower,upper));
	
	constructImage(contours[0],image1);
	constructImage(contours[1],image2);

	constructImageBothContour(contours[0],contours[1],imageCB);

	Board2D boardCB;
	Display2DFactory::drawImage<Gray>(boardCB, imageCB, (unsigned int)0, (unsigned int)129);
	boardCB.saveSVG("../../../contourBoth.svg"); 
	
	Board2D board1,board2;
	Display2DFactory::drawImage<Gray>(board1, image1, (unsigned int)0, (unsigned int)129);
	Display2DFactory::drawImage<Gray>(board2, image2, (unsigned int)0, (unsigned int)129);
	board1.saveSVG("../../../contour1.svg");
	board2.saveSVG("../../../contour2.svg");

	//Calcul de la transformée en distance 
	
	computeDT(image1,"../../../DT_contour1.svg"); 
	computeDT(image2,"../../../DT_contour2.svg"); 
	computeDT(imageCB,"../../../DT_contourBoth.svg"); 

	// Calcul de la Voronoi Map 

	computeVoronoiMap(image1,contours[0],"../../../VM_contour1.svg","../../../VMCells_contour1.svg"); 
	computeVoronoiMap(image2,contours[1],"../../../VM_contour2.svg","../../../VMCells_contour2.svg"); 
	computeVoronoiMap(imageCB,contours,"../../../VM_contourBoth.svg","../../../VMCells_contourBoth.svg"); 

	//Calcul de la moyenne des 2DT
	Image DT1 = imageDT(image1,valMaxX,valMaxY);
	Image DT2 = imageDT(image2,valMaxX,valMaxY);
	Image DTAv(Z2i::Domain(lower,upper));
	Image DTDiff(Z2i::Domain(lower,upper));
	Board2D boardAv,boardDiff;

	for (int i = 0;i<valMaxX; i++){
		for (int j = 0;j<valMaxY; j++){
			Z2i::Point p(i,j);
			DTAv.setValue(p,(DT1(p)+DT2(p))/2);
			DTDiff.setValue(p,abs(DT1(p)-DT2(p)) );
		} 
	} 
	Display2DFactory::drawImage<Gray>(boardAv, DTAv, (unsigned int)0, (unsigned int)129);
	Display2DFactory::drawImage<Gray>(boardDiff, DTDiff, (unsigned int)0, (unsigned int)129);
	boardAv.saveSVG("../../../DT_contourAv.svg");
	boardDiff.saveSVG("../../../DT_contourDiff.svg");
	
	return 0; 
}
