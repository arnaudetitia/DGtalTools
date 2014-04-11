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
#include "DGtal/topology/MetricAdjacency.h"
#include "DGtal/topology/DomainMetricAdjacency.h"
#include "DGtal/topology/DomainAdjacency.h"
#include "DGtal/topology/DigitalTopology.h"
#include "DGtal/topology/Object.h"

#include "DGtal/topology/helpers/Surfaces.h"

#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/Shapes.h"

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

#include "DGtal/graph/Expander.h"


using namespace std;
using namespace DGtal;
typedef SpaceND<2> Z2;
typedef HyperRectDomain< Z2 > Domain; 
typedef DigitalSetSelector < Domain, BIG_DS + HIGH_BEL_DS >::Type DigitalSet;
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
typedef MetricAdjacency<Z2, 1> Adj4;
typedef MetricAdjacency<Z2, 2> Adj8;
typedef DigitalTopology< Adj8, Adj4 > DT8_4;
typedef DigitalTopology< Adj4, Adj8 > DT4_8;
typedef Object<DT4_8, DigitalSet> ObjectType;
typedef Object<DT4_8, DigitalSet> ObjectType48;
typedef typename DigitalSet::ConstIterator DigitalSetConstIterator;

void readData(char *filename,vector<int>&x ,vector<int>&y ,vector<int>&noiseLevel,vector<int>&freeman){
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
				if(i==2){freeman.push_back(atoi(buffer));}
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

Z2i::Point regard(int i){
	switch(i){
		case 0: return (Z2i::Point(1,0));
		case 1: return (Z2i::Point(0,1)); 
		case 2: return (Z2i::Point(-1,0)); 
		case 3: return (Z2i::Point(0,-1));
	}
}

void getContours(vector< vector< Z2i::Point >  >&  contours,Image &image,vector<int> & x,vector<int> & y,vector<int> &noise){
	
	DigitalSet image_set ( image.domain() );
	int maxX = image.domain().upperBound()[0];
	int maxY = image.domain().upperBound()[1];
	for (unsigned int i = 0;i<maxX; i++){
		for (unsigned int j = 0;j<maxY; j++){
			Z2i::Point p(i,j);
			if (image(p) == 0) image_set.insertNew(p); 
    		}
	}
	Adj4 adj4;
	Adj8 adj8;
	DT4_8 dt4_8 ( adj4, adj8, JORDAN_DT );
	ObjectType im ( dt4_8, image_set );
	ObjectType imageBorder = im.border();
	Image tmp(image.domain());
	Board2D board;
	for ( Image::Iterator it = tmp.begin(), itend = tmp.end();it != itend; ++it)
    		(*it)=128;
	DigitalSetConstIterator it_end = imageBorder.end();
	DigitalSetConstIterator it_begin = imageBorder.begin(); 
	int currentX,currentY;
	vector<Z2i::Point> cont1;
	vector<Z2i::Point> cont2;
	for ( DigitalSetConstIterator it = it_begin;it != it_end;++it ) {
		tmp.setValue((*it),0);	
	}
	Z2i::Point depart = Z2i::Point((*it_begin)[0],(*it_begin)[1]);
	Z2i::Point currentPT = depart;
	int direction = 0;
	do{
		for (int i=0;i<4;i++){ 
			if (tmp(currentPT+regard((direction+i)%4)) == 0) {
				currentPT += regard((direction+i)%4);
				cont1.push_back(currentPT); 
				direction = (direction +i+3)%4;
				break;
			}  
		} 
	}while (currentPT[0] != depart[0] || currentPT[1] != depart[1]);
	cont1.push_back(currentPT);
	int t;
	int nbVoisin;
	for ( DigitalSetConstIterator it = it_begin;it != it_end;++it ) {
		currentX = (*it)[0];
		currentY = (*it)[1];
		for (int i = 0; i<x.size();i++){
			if ((x[i] == currentX) && (y[i] == maxY-currentY)){
				if (noise[i] == 1){
					Z2i::Point pt = Z2i::Point((*it)); 
					if (tmp(pt) == 0){
						nbVoisin = 0;
						if (tmp(pt+regard(0)) == 0) nbVoisin++;
						if (tmp(pt+regard(1)) == 0) nbVoisin++;
						if (tmp(pt+regard(2)) == 0) nbVoisin++;
						if (tmp(pt+regard(3)) == 0) nbVoisin++;						
					} 
				} 
			}
			if  (nbVoisin == 2) {depart= (*it);break;}
		}   	
	} 
	currentPT = depart;
	direction = 0;
	do{
		for (int i=0;i<3;i++){ 
			if (tmp(currentPT+regard((direction+i)%4)) == 0) {
				t = i;
			}	 
		} 
		currentPT += regard((direction+t)%4);
		cont2.push_back(currentPT); 
		direction = (direction+t+3)%4;
	}while (currentPT[0] != depart[0] || currentPT[1] != depart[1]);
	cont2.push_back(currentPT); 
	/*int trouve1,trouve2; 
	for ( DigitalSetConstIterator it = it_begin;it != it_end;++it ) {
		currentX = (*it)[0];
		currentY = (*it)[1];
		trouve1 = 0;
		for (int i = 0;i<x.size();i++){
			if ((x[i] == currentX) && (y[i] == maxY-currentY)) {
				cont2.push_back((*it));
				trouve1 =1;
				break; 
			}
		}
		if (!trouve1){
			trouve2 = 0;
			for (int i = 0;i<cont1.size();i++){
				if ( (cont1[i][0] == currentX)  && (cont1[i][1] == maxYcurrentY) ) {
					trouve2 = 1;
					break;
				}
			}
			if (!trouve2) cont2.push_back((*it));	
		} 
	}*/
	contours.push_back(cont1);
	contours.push_back(cont2);
	Display2DFactory::drawImage<Gray>(board, tmp, (unsigned int)0, (unsigned int)129);
	
	board.saveSVG ( "../../../imageBorder.svg" );
	
} 

void constructImage(vector< Z2i::Point >& contour,Image& image){
	for ( Image::Iterator it = image.begin(), itend = image.end();it != itend; ++it)
    		(*it)=128;

	for (int i=0;i<contour.size();i++){
		image.setValue(contour[i],0);
	}
	
}

void contourBallReduced(vector<int>& x,vector<int>& y,vector<int>& noiseLevel,vector<int>& freeman,int maxX,int maxY){
	Z2i::Point lower(0,0);
	Z2i::Point upper(maxX,maxY);
	Image image ( Z2i::Domain(lower,upper));

	for ( Image::Iterator it = image.begin(), itend = image.end();it != itend; ++it)
    		(*it)=128;
	int cpt = 0;
	for (int i =1;i<x.size()-1;i++){
		if (
		   	( 
			  (noiseLevel[i] >= 2) &&  
		     	  (noiseLevel[i] > noiseLevel[i-1]) && 
		     	  (noiseLevel[i] > noiseLevel[i+1]) 
		   	)
		   || 
		   	(
		     		(
					(
						(noiseLevel[i] == noiseLevel[i-1]) && 
		       				(noiseLevel[i] - noiseLevel[i+1] > 2 )
					) 
		    			|| 
					(
						(noiseLevel[i] == noiseLevel[i+1]) && 
		       				(noiseLevel[i] - noiseLevel[i-1] > 2 )
					) 
		   		)
			)
		   ) cpt ++;
	}

	if (
		   ( 
			 (noiseLevel[0] >= 2) &&  
		     	 (noiseLevel[0] > noiseLevel[x.size()-1]) && 
		     	 (noiseLevel[0] > noiseLevel[1]) 
		   )
		   || 
		   (
		     	(
				(
					(noiseLevel[0] == noiseLevel[x.size()-1]) && 
		       			(noiseLevel[0] - noiseLevel[1] > 2 )
				) 
		    		|| 
				(
					(noiseLevel[0] == noiseLevel[1]) && 
		       			(noiseLevel[0] - noiseLevel[x.size()-1] > 2 )
				) 
		   	)
		)
	    ) cpt ++;
	
	if (
		   ( 
			 (noiseLevel[x.size()-1] >= 2) &&  
		     	 (noiseLevel[x.size()-1] > noiseLevel[x.size()-2]) && 
		     	 (noiseLevel[x.size()-1] > noiseLevel[0]) 
		   )
		   || 
		   (
		     	(
				(
					(noiseLevel[x.size()-1] == noiseLevel[x.size()-2]) && 
		       			(noiseLevel[x.size()-1] - noiseLevel[0] > 2 )
				) 
		    		|| 
				(
					(noiseLevel[x.size()-1] == noiseLevel[0]) && 
		       			(noiseLevel[x.size()-1] - noiseLevel[x.size()-2] > 2 )
				) 
		   	)
		)
	    ) cpt ++;

	cout << "Nombre de cas interessant : " << cpt << " sur " << x.size() << endl; 
	cout << "Soit un Ratio de "<< 100 * cpt / x.size() << "%" <<  endl;
		/*if (freeman[i-1] == freeman[i]){
			if (noiseLevel[i] == 1 || (noiseLevel[i-1] == 1 && noiseLevel[i+1] == 1) ){
				image.setValue( Z2i::Point(x[i],y[i]),0 );
			}
			else{
				if ((noiseLevel[i] <= noiseLevel[i-1]) && (noiseLevel[i] <= noiseLevel[i+1] ) ){
					for (int a = max(0,x[i]-(noiseLevel[i]-1) ) ; a <= min(x[i]+(noiseLevel[i]-1),maxX); a++ ){
					for (int b = max(0,y[i]-(noiseLevel[i]-1) ) ; b <= min(y[i]+(noiseLevel[i]-1),maxY); b++ ){
						image.setValue( Z2i::Point(a,maxY-b+20),0);
					}
					}
				}
				if ((noiseLevel[i] > noiseLevel[i-1]) && 
				    (noiseLevel[i] > noiseLevel[i+1] ) && 
				    (noiseLevel[i-1] ==  noiseLevel[i+1]) ){
				for (int a = max(0,x[i-1]-(noiseLevel[i-1]-1) ) ; a <= min(x[i+1]+(noiseLevel[i-1]-1),maxX); a++ ){
				for (int b = max(0,y[i-1]-(noiseLevel[i-1]-1) ) ; b <= min(y[i+1]+(noiseLevel[i-1]-1),maxY); b++ ){
					image.setValue( Z2i::Point(a,maxY-b+20),0);
				}
				}
				for (int a = max(0,x[i+1]-(noiseLevel[i+1]-1) ) ; a <= min(x[i+1]+(noiseLevel[i+1]-1),maxX); a++ ){
				for (int b = max(0,y[i+1]-(noiseLevel[i+1]-1) ) ; b <= min(y[i+1]+(noiseLevel[i+1]-1),maxY); b++ ){
					image.setValue( Z2i::Point(a,maxY-b+20),0);
				}
				}
				}
			}
		}
		else {
			if (noiseLevel[i] == 1 ){
				image.setValue( Z2i::Point(x[i],y[i]),0 );
			}
			if (noiseLevel[i-1] == 1 && noiseLevel[i+1] == 1) ){
				image.setValue( Z2i::Point(x[i],y[i]),0 );
				//image.setValue( Z2i::Point(x[i],y[i]),0 );
			}
			
		}
	} */
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

Image imageDT(Image& image,unsigned int maxX, unsigned int maxY){
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

void computeDTAverage(Image & image1,Image & image2,int maxX,int maxY,char* outputfile){
	Z2i::Point lower(0,0);
	Z2i::Point upper(maxX,maxY);
	Image DTAv(Z2i::Domain(lower,upper));
	Board2D boardAv;
	for (int i = 0;i<=maxX; i++){
		for (int j = 0;j<=maxY; j++){
			Z2i::Point p(i,j);
			DTAv.setValue(p,(image1(p)+image2(p))/2);
		} 
	} 

	unsigned int maxAv=0;
  	for ( Image::iterator it = DTAv.begin(), itend = DTAv.end();it != itend; ++it)
    		if ( (*it) > maxAv)  maxAv = (*it);

	Display2DFactory::drawImage<HueTwice>(boardAv, DTAv, (unsigned int)0, maxAv +1);
	boardAv.saveSVG(outputfile);
}

void computeDTDiff(Image & image1,Image & image2,int maxX,int maxY,char* outputfile){
	Z2i::Point lower(0,0);
	Z2i::Point upper(maxX,maxY);
	Image DTDiff(Z2i::Domain(lower,upper));
	Board2D boardDiff;
	for (unsigned int i = 0;i<maxX; i++){
		for (unsigned int j = 0;j<maxY; j++){
			Z2i::Point p(i,j);
			if (image1(p)>image2(p)) DTDiff.setValue(p,image1(p)-image2(p));
			else DTDiff.setValue(p,image2(p)-image1(p));
		} 
	} 

	unsigned int maxDiff=0;
  	for ( Image::iterator it = DTDiff.begin(), itend = DTDiff.end();it != itend; ++it)
    		if ( (*it) > maxDiff)  maxDiff = (*it);

	Display2DFactory::drawImage<Gray>(boardDiff, DTDiff, (unsigned int)0, maxDiff +1);
	boardDiff.saveSVG(outputfile);
}

int main(int argc,char **argv){
	if (argc != 2){
		cout << "Usage : ./impreciseBounadary <filename>" << endl;
		exit(1);
	}	 
	char *filename = argv[1];
	vector<int>x ;
	vector<int>y ;
	vector<int>noiseLevel;
	vector<int>freeman;
	readData(filename,x,y,noiseLevel,freeman);// récupération des données de contours 
	//x.pop_back();
	//y.pop_back();
	//noiseLevel.pop_back();
	//freeman.pop_back();
	

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
		for (int a = max(0,x[i]-(noiseLevel[i]-1) ) ; a <= min(x[i]+(noiseLevel[i]-1),valMaxX+20); a++ ){
			for (int b = max(0,y[i]-(noiseLevel[i]-1) ) ; b <= min(y[i]+(noiseLevel[i]-1),valMaxY+20); b++ ){
				image.setValue( Z2i::Point(a,valMaxY-b+20),0);
				//image.setValue( Z2i::Point(x[i],valMaxY-y[i]+20),0);
			}
		}
	}
	// affichage de l'image dans un svg	
	

	Board2D board;
	Display2DFactory::drawImage<Gray>(board, image, (unsigned int)0, (unsigned int)129);
	board.saveSVG("../../../monResultat.svg");
	
	//récupération des contours intérieurs et extérieurs
	
	vector< vector< Z2i::Point >  >  contours;
	
	getContours(contours,image,x,y,noiseLevel);
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
	
	/*computeDT(image1,"../../../DT_contour1.svg"); 
	computeDT(image2,"../../../DT_contour2.svg"); 
	computeDT(imageCB,"../../../DT_contourBoth.svg"); 

	// Calcul de la Voronoi Map 

	computeVoronoiMap(image1,contours[0],"../../../VM_contour1.svg","../../../VMCells_contour1.svg"); 
	computeVoronoiMap(image2,contours[1],"../../../VM_contour2.svg","../../../VMCells_contour2.svg"); 
	computeVoronoiMap(imageCB,contours,"../../../VM_contourBoth.svg","../../../VMCells_contourBoth.svg"); 

	//Calcul de la moyenne et la différence des 2DT
	Image DT1 = imageDT(image1,valMaxX+20,valMaxY+20);
	Image DT2 = imageDT(image2,valMaxX+20,valMaxY+20);
	computeDTAverage(DT1,DT2,valMaxX+20,valMaxY+20,"../../../DT_contourAv.svg");
	computeDTDiff(DT1,DT2,valMaxX+20,valMaxY+20,"../../../DT_contourDiff.svg");

	//Création d'un contour avec boules réduites
	contourBallReduced(x,y,noiseLevel,freeman,valMaxX+20,valMaxY+20); */
	
	return 0; 
}
