#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <string>
#include <math.h>
#include <queue>

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
#include "DGtal/images/CImage.h"

#include "DGtal/io/colormaps/GrayscaleColorMap.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"

#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpPowerSeparableMetric.h"
#include "DGtal/geometry/volumes/distance/VoronoiMap.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/geometry/helpers/ContourHelper.h"
#include "DGtal/geometry/volumes/distance/ReverseDistanceTransformation.h"
#include "DGtal/geometry/volumes/distance/ReducedMedialAxis.h"
#include "DGtal/geometry/volumes/distance/PowerMap.h"

#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include "DGtal/io/boards/Board2D.h"

#include "DGtal/graph/Expander.h"


using namespace std;
using namespace DGtal;
typedef SpaceND<2> Z2;
typedef HyperRectDomain< Z2 > Domain; 
typedef DigitalSetSelector < Domain, BIG_DS + HIGH_BEL_DS >::Type DigitalSet;
typedef GrayscaleColorMap<unsigned char> Gray;
typedef ImageSelector<Z2i::Domain, unsigned int>::Type Grille;// C'est une Image mais c'est pour éviter les conflits avec le fichier  ReducedMedialAxis.h
typedef ImageContainerBySTLMap<DigitalSetDomain<Z2i::DigitalSet> , DGtal::int64_t> MapImage;
typedef IntervalThresholder<Grille::Value> Binarizer; 
typedef SimpleThresholdForegroundPredicate<Grille> PointPredicate;
typedef DistanceTransformation<Z2i::Space, PointPredicate, Z2i::L2Metric> DTL2;
typedef DistanceTransformation<Z2i::Space, PointPredicate, Z2i::L2Metric> VML2;
typedef HueShadeColorMap<long int, 2> HueTwice;
typedef NotPointPredicate<Z2i::DigitalSet> NotPredicate;
typedef ExactPredicateLpSeparableMetric<Z2i::Space, 2> L2Metric;
typedef ExactPredicateLpPowerSeparableMetric<Z2i::Space,2 > LIPowerMetric;
typedef VoronoiMap<Z2i::Space, NotPredicate, L2Metric > Voronoi2D;
typedef MetricAdjacency<Z2, 1> Adj4;
typedef MetricAdjacency<Z2, 2> Adj8;
typedef DigitalTopology< Adj8, Adj4 > DT8_4;
typedef DigitalTopology< Adj4, Adj8 > DT4_8;
typedef Object<DT4_8, DigitalSet> ObjectType;
typedef Object<DT4_8, DigitalSet> ObjectType48;
typedef typename DigitalSet::ConstIterator DigitalSetConstIterator;
typedef ReverseDistanceTransformation< Grille , LIPowerMetric > RDT;


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

void getContours(vector< vector< Z2i::Point >  >&  contours,Grille &image,vector<int> & x,vector<int> & y,vector<int> &noise){
	
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
	Grille tmp(image.domain());
	Board2D board;
	for ( Grille::Iterator it = tmp.begin(), itend = tmp.end();it != itend; ++it)
    		(*it)=128;
	DigitalSetConstIterator it_end = imageBorder.end();
	DigitalSetConstIterator it_begin = imageBorder.begin(); 
	int currentX,currentY;
	vector<Z2i::Point> cont1;
	vector<Z2i::Point> cont2;
	for ( DigitalSetConstIterator it = it_begin;it != it_end;++it ) {
		tmp.setValue((*it),0);	
	}
	
	Z2i::Point depart = (*it_begin);
	Z2i::Point currentPT = depart;
	int flag;
	int nbVoisin;
	int direction = 0;
	// extraction du contour exterieur
	
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
	// extraction des point qui ne sont pas sur le contour exterieur
	vector<Z2i::Point> mand;
	for (int i = 0;i<=maxX;i++){
		for (int j = 0;j<=maxY;j++){
			if (tmp(Z2i::Point(i,j)) == 0){
				flag = 0;
				for (int k=0;k<cont1.size();k++){
					if (cont1[k][0] == i && cont1[k][1] == j){
						flag = 1;
						break; 
					}
				}
				if (!flag) mand.push_back(Z2i::Point(i,j)); 
			} 
		} 
	}
	// extraction des points sont sur le contour exterieur mais qui ne sont pas sur le contour original
	vector<Z2i::Point> forbid;
	for (int i = 0;i<cont1.size();i++){
		flag = 0;
		for (int j = 0;j<x.size();j++){
			if (cont1[i][0] == x[j] && cont1[i][1] == maxY-y[j]){
				flag = 1;
				break;
			} 
		}
		if (!flag) forbid.push_back(Z2i::Point(cont1[i]));
	}
	// recherche d'un point de départ pour le contour interieur
	int a=0;
	int t;
	flag = 0;
	Z2i::Point h;
	while (!flag){
		h[0] = x[a];
		h[1] = maxY -y[a];
		if (image(h) != 0){
			a++;
		} 
		else {
			nbVoisin = 0;
			if (image(h+regard(0)) == 0 ) nbVoisin++;
			if (image(h+regard(1)) == 0 ) nbVoisin++;
			if (image(h+regard(2)) == 0 ) nbVoisin++;
			if (image(h+regard(3)) == 0 ) nbVoisin++;
			if (nbVoisin == 2 ) {
				if (image(h+Z2i::Point(1,1)) == 0 || 
				    image(h+Z2i::Point(-1,1)) == 0 || 
				    image(h+Z2i::Point(-1,-1)) == 0 || 
				    image(h+Z2i::Point(1,-1)) == 0 ) a++;
				else flag = 1;
			}
			else a++;
		} 
	} 
	depart = Z2i::Point(x[a],maxY-y[a]);
	currentPT = depart;
	direction = 0;
	vector<Z2i::Point> pts;
	vector<int> dir;
	 
	do{
		pts.clear();
		dir.clear();
		for (int i=0;i<4;i++){ 
			if (tmp(currentPT+regard((direction+i)%4)) == 0) {
				pts.push_back(currentPT+regard((direction+i)%4));
				dir.push_back((direction+i)%4); 
			} 
		} 
		for (int i =0;i< pts.size();i++){
			flag = 0;
			for (int j =0;j< mand.size();j++){
				if (pts[i] == mand[j]){
					mand.erase(mand.begin()+j);
					cont2.push_back(pts[i]);
					currentPT = pts[i];
					direction = (dir[i]+3)%4;
					flag = 1;
					break;
				} 
			}
			if (flag) break;  
		}
		if (!flag){
			for (int i =0;i< pts.size();i++){
				flag = 0;
				for (int j =0;j< forbid.size();j++){
					if (pts[i] == forbid[j]){
						flag = 1;
						break;
					} 	
				}
				if (!flag) {
					cont2.push_back(pts[i]);
					currentPT = pts[i];
					direction = (dir[i]+3)%4;
					break;
				}
			}	
		}	
	}while (currentPT != depart && cont2.size() < image.domain().upperBound()[0] * image.domain().upperBound()[1] ); 
	cont2.push_back(currentPT); 
	contours.push_back(cont1);
	contours.push_back(cont2);
	Display2DFactory::drawImage<Gray>(board, tmp, (unsigned int)0, (unsigned int)129);
	
	board.saveSVG ( "../../../imageBorder.svg" );
	
}

bool isInside(Grille &image,Z2i::Point p){
	int limX = image.domain().upperBound()[0];
	int limY = image.domain().upperBound()[1];
	int L=0,R=0,U=0,D=0;
	int a,b;
	for (int i = p[0]+1 ; i<=limX ; i++){
		if (image(Z2i::Point(i,p[1])) == 128 ){
			if (( image(Z2i::Point(i,p[1]+1)) == 128 ) && ( image(Z2i::Point(i,p[1]-1)) == 128 )) R++;
			else {
				if ( image(Z2i::Point(i,p[1]+1)) == 128 ) a = p[1]+1;
				else a = p[1]-1;
				while ( image(Z2i::Point(i+1,p[1])) == 128 ) i++;
				if ( image(Z2i::Point(i,p[1]+1)) == 128 ) b = p[1]+1;
				else b = p[1]-1;
				if (a != b ) R++; 	 
			}
		}  
	}
	for (int i = p[0]-1 ; i>=0 ; i--){
		if (image(Z2i::Point(i,p[1])) == 128 ){
			if (( image(Z2i::Point(i,p[1]+1)) == 128 ) && ( image(Z2i::Point(i,p[1]-1)) == 128 )) L++;
			else {
				if ( image(Z2i::Point(i,p[1]+1)) == 128 ) a = p[1]+1;
				else a = p[1]-1;
				while ( image(Z2i::Point(i-1,p[1])) == 128 ) i--;
				if ( image(Z2i::Point(i,p[1]+1)) == 128 ) b = p[1]+1;
				else b = p[1]-1;
				if (a != b ) L++; 	 
			}
		}  
	}
	for (int i = p[1]+1 ; i<=limY ; i++){
		if (image(Z2i::Point(p[0],i)) == 128 ){
			if (( image(Z2i::Point(p[0]+1,i)) == 128 ) && ( image(Z2i::Point(p[0]-1,i)) == 128 )) D++;
			else {
				if ( image(Z2i::Point(p[0]+1,i)) == 128 ) a = p[0]+1;
				else a = p[0]-1;
				while ( image(Z2i::Point(p[0],i+1)) == 128 ) i++;
				if ( image(Z2i::Point(p[0]+1,i)) == 128 ) b = p[0]+1;
				else b = p[0]-1;
				if (a != b ) D++; 	 
			}
		}  
	}
	for (int i = p[1]-1 ; i>=0 ; i--){
		if (image(Z2i::Point(p[0],i)) == 128 ){
			if (( image(Z2i::Point(p[0]+1,i)) == 128 ) && ( image(Z2i::Point(p[0]-1,i)) == 128 )) U++;
			else {
				if ( image(Z2i::Point(p[0]+1,i)) == 128 ) a = p[0]+1;
				else a = p[0]-1;
				while ( image(Z2i::Point(p[0],i-1)) == 128 ) i--;
				if ( image(Z2i::Point(p[0]+1,i)) == 128 ) b = p[0]+1;
				else b = p[0]-1;
				if (a != b ) U++; 	 
			}
		}  
	}
	return ( (U%2==1) && (R%2==1) && (L%2==1) && (D%2==1) );	
}

void fill(Grille &image,Z2i::Point p){
	int limX = image.domain().upperBound()[0];
	int limY = image.domain().upperBound()[1];
	for (int i =0; i<4;i++){
		Z2i::Point pt = p+regard(i);
		if (image(pt) == 0 && pt[0]<=limX && pt[1]<=limY) {
			image.setValue(pt,128);
			fill(image,pt);
		} 
	} 
}  
  

void constructImage(vector< Z2i::Point >& contour,Grille& image){
	for ( Grille::Iterator it = image.begin(), itend = image.end();it != itend; ++it)
    		(*it)=0;
	 
	for (int i=0;i<contour.size();i++){
		image.setValue(contour[i],128);
	}
	srand (time(NULL));
	int limX = image.domain().upperBound()[0];
	int limY = image.domain().upperBound()[1];
	int x = rand()%limX,y=rand()%limY;
	while (!(isInside(image,Z2i::Point(x,y))) || (image(Z2i::Point(x,y)) == 128 ) ) {
		x = rand()%limX;
		y = rand()%limY;
	}
	image.setValue(Z2i::Point(x,y),128);
	fill(image,Z2i::Point(x,y)); 
	
	
}

void contourBallReduced(vector<int>& x,vector<int>& y,vector<int>& noiseLevel,vector<int>& freeman,int maxX,int maxY){
	Z2i::Point lower(0,0);
	Z2i::Point upper(maxX,maxY);
	Grille image ( Z2i::Domain(lower,upper));

	for ( Grille::Iterator it = image.begin(), itend = image.end();it != itend; ++it)
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

void constructImageBothContour(vector< Z2i::Point >& contour1,vector< Z2i::Point >& contour2,Grille& image){
	for ( Grille::Iterator it = image.begin(), itend = image.end();it != itend; ++it)
    		(*it)=128;

	for (int i=0;i<contour1.size();i++){
		image.setValue(contour1[i],0);
	}
	for (int i=0;i<contour2.size();i++){
		image.setValue(contour2[i],0);
	}
	
	
}

void computeDT(Grille& image,char *outputfile){
	PointPredicate predicate(image,0);
	DTL2 dt(image.domain(), predicate, Z2i::l2Metric);

	DTL2::Value maxv=0;
  	for ( DTL2::ConstRange::ConstIterator it = dt.constRange().begin(), itend = dt.constRange().end();it != itend; ++it)
    		if ( (*it) > maxv)  maxv = (*it);

	Board2D board;
	board.clear();
  	Display2DFactory::drawImage<HueTwice>(board, dt, 0, maxv + 1); 

	board.saveSVG ( outputfile );
}

Grille imageDT(Grille& image,unsigned int maxX, unsigned int maxY){
	PointPredicate predicate(image,0);
	Z2i::Point lower(0,0);
	Z2i::Point upper(maxX,maxY);
	Grille result(Z2i::Domain(lower,upper));
	DTL2 dt(image.domain(), predicate, Z2i::l2Metric);

	for (int i=0;i<maxX;i++){
		for (int j=0;j<maxY;j++){
			Z2i::Point p(i,j);
			result.setValue(p,dt(p));
		}
	}
	return result;
} 

void computeVoronoiMap(Grille& image,vector< Z2i::Point >& contour,char *outputfile,char* outputfileCells){
	
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

void computeVoronoiMap(Grille& image,vector< vector< Z2i::Point > >& contours,char *outputfile,char* outputfileCells){
	
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

void computeDTAverage(Grille & image1,Grille & image2,int maxX,int maxY,char* outputfile){
	Z2i::Point lower(0,0);
	Z2i::Point upper(maxX,maxY);
	Grille DTAv(Z2i::Domain(lower,upper));
	Board2D boardAv;
	for (int i = 0;i<=maxX; i++){
		for (int j = 0;j<=maxY; j++){
			Z2i::Point p(i,j);
			DTAv.setValue(p,(image1(p)+image2(p))/2);
		} 
	} 

	unsigned int maxAv=0;
  	for ( Grille::iterator it = DTAv.begin(), itend = DTAv.end();it != itend; ++it)
    		if ( (*it) > maxAv)  maxAv = (*it);

	Display2DFactory::drawImage<HueTwice>(boardAv, DTAv, (unsigned int)0, maxAv +1);
	boardAv.saveSVG(outputfile);
}

void computeDTDiff(Grille & image1,Grille & image2,int maxX,int maxY,char* outputfile){
	Z2i::Point lower(0,0);
	Z2i::Point upper(maxX,maxY);
	Grille DTDiff(Z2i::Domain(lower,upper));
	Board2D boardDiff;
	for (unsigned int i = 0;i<=maxX; i++){
		for (unsigned int j = 0;j<=maxY; j++){
			Z2i::Point p(i,j);
			if (image1(p)>image2(p)) DTDiff.setValue(p,image1(p)-image2(p));
			else DTDiff.setValue(p,image2(p)-image1(p));
		} 
	} 

	unsigned int maxDiff=0;
  	for ( Grille::iterator it = DTDiff.begin(), itend = DTDiff.end();it != itend; ++it)
    		if ( (*it) > maxDiff)  maxDiff = (*it);

	Display2DFactory::drawImage<Gray>(boardDiff, DTDiff, (unsigned int)0, maxDiff +1);
	boardDiff.saveSVG(outputfile);
}

void computeReverseDT(Grille& image,char* outputfile){
	
	Z2i::L2PowerMetric l2power;

	RDT reverseDT(&image.domain(),&image,&l2power);
	RDT::Value maxv=0;
	RDT::Value minv= 500000;
	for ( RDT::ConstRange::ConstIterator it = reverseDT.constRange().begin(), itend = reverseDT.constRange().end();it != itend; ++it)
    		{if ((*it) < minv)  minv = (*it);if ((*it) > maxv)  maxv = (*it);}
	
	Grille result (image.domain());
	//RDT::ConstRange::ConstIterator it = reverseDT.constRange().begin();
	for (int i = 0;i<=image.domain().upperBound()[1];i++){
		for (int j = 0;j<=image.domain().upperBound()[0];j++){
			Z2i::Point p(j,i); 
			if (image(p) != 0 ){
				for (int a = 0;a<=image.domain().upperBound()[1];a++){
					for (int b = 0;b<=image.domain().upperBound()[0];b++){
						if (abs(j-b)*abs(j-b) + abs(i-a)*abs(i-a)< image(p) ) result.setValue( Z2i::Point(b,a),1); 
					}
				}
			} 
		} 
	}	 
	Board2D board;
	board.clear();
	Display2DFactory::drawImage<Gray>(board, result ,0, 1);
	board.saveSVG(outputfile);
	
} 

void computeMA(Grille &image,char *outputfile){
	Z2i::Domain domain(Z2i::Point(0,0),image.domain().upperBound());
  	Z2i::Domain domainLarge(Z2i::Point(0,0),image.domain().upperBound());
	
	Grille im(domain);
	for (int i = 0;i<=image.domain().upperBound()[0];i++){
		for (int j = 0;j<=image.domain().upperBound()[1];j++){
			Z2i::Point p(i,j); 
			im.setValue(p,image(p) * image(p));  
		} 
	}
	 
	
	Z2i::L2PowerMetric l2power;
  	PowerMap<Grille, Z2i::L2PowerMetric> power(&domainLarge, &im, &l2power);
	
	ReducedMedialAxis<PowerMap<Grille, Z2i::L2PowerMetric> >::Type  rdma = ReducedMedialAxis< PowerMap<Grille, Z2i::L2PowerMetric> >::getReducedMedialAxisFromPowerMap(power);

	int maxv=0;
	for (int i = 0;i<=im.domain().upperBound()[0];i++){
		for (int j = 0;j<=im.domain().upperBound()[1];j++){
			if (rdma(Z2i::Point(i,j)) > maxv) maxv = rdma(Z2i::Point(i,j)); 
		} 
	}
	
	Board2D board;
	board.clear();
	Display2DFactory::drawImage<Gray>(board, rdma , (unsigned int)0, maxv +1);
	board.saveSVG(outputfile);
}

Grille imageAM(Grille& image){
	Z2i::Domain domain(Z2i::Point(0,0),image.domain().upperBound());
  	Z2i::Domain domainLarge(Z2i::Point(0,0),image.domain().upperBound());
	
	Grille im(domain);
	for (int i = 0;i<=image.domain().upperBound()[0];i++){
		for (int j = 0;j<=image.domain().upperBound()[1];j++){
			Z2i::Point p(i,j); 
			im.setValue(p,image(p) * image(p));  
		} 
	}
	 
	
	Z2i::L2PowerMetric l2power;
  	PowerMap<Grille, Z2i::L2PowerMetric> power(&domainLarge, &im, &l2power);
	
	ReducedMedialAxis<PowerMap<Grille, Z2i::L2PowerMetric> >::Type  rdma = ReducedMedialAxis< PowerMap<Grille, Z2i::L2PowerMetric> >::getReducedMedialAxisFromPowerMap(power);
	
	Grille result(Z2i::Domain(image.domain().lowerBound(),image.domain().upperBound()));
	for (int i = 0;i<=result.domain().upperBound()[0];i++){
		for (int j = 0;j<=result.domain().upperBound()[1];j++){
			Z2i::Point p(i,j);
			result.setValue(p,rdma(p));  
		} 
	}
	return result;
} 

void homotopicThinning(Grille &image,char *outputfile){
	PointPredicate predicate(image,0);
	Z2i::Point lower(0,0);
	Z2i::Point upper(image.domain().upperBound()[0],image.domain().upperBound()[1]);
	DTL2 dt(image.domain(), predicate, Z2i::l2Metric);
	Domain domain( lower, upper );
	DigitalSet set( domain );
	for (int i = 0; i<=image.domain().upperBound()[0] ; i++){
		for (int j = 0; j<=image.domain().upperBound()[1] ; j++){
			Z2i::Point p(i,j); 
			if (image(p) != 0) set.insertNew(p);   
		}
	}
	ObjectType48 shape(Z2i::dt4_8 , set);
	int nb_simple=0; 
  	int layer = 1;
  	std::queue<DigitalSet::Iterator> Q;
	do {
      		int nb=0;
      		DigitalSet & S = shape.pointSet();
 
      
      		for ( DigitalSet::Iterator it = S.begin(); it != S.end(); ++it ){
	  		if (dt( *it ) <= layer){
	      			if ( shape.isSimple( *it ) )
					Q.push( it );
	    		}
		}
      		nb_simple = 0;
      		while ( ! Q.empty() ){
          		DigitalSet::Iterator it = Q.front();
          		Q.pop();
          		if ( shape.isSimple( *it ) ){
              			S.erase( *it );
              			++nb_simple;
            		}
        	}
      		++layer;
     	}
 	while ( nb_simple != 0 );

  	Grille imageRes(domain);

	for ( Grille::Iterator it = imageRes.begin(), itend = imageRes.end();it != itend; ++it)
    		(*it)=0;

     	DigitalSet & S = shape.pointSet();
	for ( DigitalSet::Iterator it = S.begin(); it != S.end(); ++it )
		imageRes.setValue((*it),128);

	Board2D board;
	board.clear();
	Display2DFactory::drawImage<Gray>(board, imageRes , (unsigned int) 0 ,(unsigned int) 128);
	board.saveSVG(outputfile); 
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
	int ecartX = 10;
	int ecartY = 10;
	
	//création du domaine 
	Z2i::Point lower(0,0);
	Z2i::Point upper(valMaxX+ecartX,valMaxY+ecartY);
	
	//Création de l'image
	
  	Grille image ( Z2i::Domain(lower,upper));
	Grille test ( Z2i::Domain(lower,Z2i::Point(30,30)));
	Grille f ( Z2i::Domain(lower,upper));
	Grille imageN ( Z2i::Domain(lower,upper));

	for ( Grille::Iterator it = image.begin(), itend = image.end();it != itend; ++it)
    		(*it)=128;

	for ( Grille::Iterator it = test.begin(), itend = test.end();it != itend; ++it)
    		(*it)=0;

	test.setValue(Z2i::Point(15,15),25); 
	

	for (int i=0;i<x.size();i++){
		for (int a = max(0,x[i]-(noiseLevel[i]-1) ) ; a <= min(x[i]+(noiseLevel[i]-1),valMaxX+ecartX); a++ ){
			for (int b = max(0,y[i]-(noiseLevel[i]-1) ) ; b <= min(y[i]+(noiseLevel[i]-1),valMaxY+ecartY); b++ ){
				image.setValue( Z2i::Point(a,valMaxY-b+ecartY),0);
			}
		}
	}

	for ( Grille::Iterator it = f.begin(), itend = f.end();it != itend; ++it)
    		(*it)=0;
	
	for (int i=0;i<x.size();i++){
		f.setValue( Z2i::Point(x[i],valMaxY-y[i]+ecartY),noiseLevel[i]*noiseLevel[i]);
	}
	//f.setValue( Z2i::Point((valMaxX+ecartX)/2,(valMaxY+ecartY)/2),16);
	
	Z2i::Domain domain(Z2i::Point(0,0),image.domain().upperBound());
	Z2i::DigitalSet set(domain);

	for (int i=0;i<x.size();i++){
		set.insertNew(Z2i::Point(x[i],valMaxY-y[i]+20));
	}

	DigitalSetDomain<Z2i::DigitalSet> setDomain(set); 
	MapImage g(setDomain);
	for (int i=0;i<x.size();i++){
		g.setValue( Z2i::Point(x[i],valMaxY-y[i]+20),noiseLevel[i]*noiseLevel[i]);
	}

	// affichage de l'image dans un svg	
	
	
	Board2D board,boardN,boardG,bTest;
	Display2DFactory::drawImage<Gray>(board, image, (unsigned int)0, (unsigned int)129);
	board.saveSVG("../../../monResultat.svg");
	for (int i = 0; i < image.domain().upperBound()[0];i++){
		 for (int j = 0; j < image.domain().upperBound()[1];j++){
			imageN.setValue(Z2i::Point(i,j),128-image(Z2i::Point(i,j)) ); 
		}
	}
	Display2DFactory::drawImage<Gray>(boardN, imageN, (unsigned int)0, (unsigned int)129);
	Display2DFactory::drawImage<Gray>(bTest, test, (unsigned int)0, (unsigned int)25);
	boardN.saveSVG("../../../monResultatNegatif.svg");
	bTest.saveSVG("../../../TestOriginal.svg");
	homotopicThinning(imageN,"../../../skeleton.svg");
	//récupération des contours intérieurs et extérieurs
	
	vector< vector< Z2i::Point >  >  contours;
	
	getContours(contours,image,x,y,noiseLevel);
	Grille image1 ( Z2i::Domain(lower,upper));
	Grille image2 ( Z2i::Domain(lower,upper));
	Grille imageCB ( Z2i::Domain(lower,upper));
	
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
	computeDT(test,"../../../DT_test.svg");

	// Calcul de la Voronoi Map 

	computeVoronoiMap(image1,contours[0],"../../../VM_contour1.svg","../../../VMCells_contour1.svg"); 
	computeVoronoiMap(image2,contours[1],"../../../VM_contour2.svg","../../../VMCells_contour2.svg"); 
	computeVoronoiMap(imageCB,contours,"../../../VM_contourBoth.svg","../../../VMCells_contourBoth.svg"); 

	//Calcul de la moyenne et la différence des 2DT
	Grille DT1 = imageDT(image1,valMaxX+ecartX,valMaxY+ecartY);
	Grille DT2 = imageDT(image2,valMaxX+ecartX,valMaxY+ecartY);
	Grille DTB = imageDT(imageCB,valMaxX+ecartX,valMaxY+ecartY);
	Grille DTN = imageDT(imageN,valMaxX+ecartX,valMaxY+ecartY);
	Grille DTTest = imageDT(test,30,30);
	computeDTAverage(DT1,DT2,valMaxX+ecartX,valMaxY+ecartY,"../../../DT_contourAv.svg");
	computeDTDiff(DT1,DT2,valMaxX+ecartX,valMaxY+ecartY,"../../../DT_contourDiff.svg");

	// Calcul de la DT inverse
	cout << "reverse DT de F" << endl;
	computeReverseDT(f,"../../../RDT_contour1.svg"); 

	// Calcul de l'axe médian
	computeMA(DT1,"../../../MA_contour1.svg");
	//computeMA(DT2,"../../../MA_contour2.svg");
	Grille MA1 = imageAM(DT1);
	//Grille MA2 = imageAM(DT2);
	cout << "reverse DT de MA1" << endl;
	computeReverseDT(MA1,"../../../RDT_MA1.svg");
	computeReverseDT(test,"../../../RDT_MA_test.svg"); 
	//computeReverseDT(MA2,"../../../RDT_MA2.svg"); 
	

	//Création d'un contour avec boules réduites
	//contourBallReduced(x,y,noiseLevel,freeman,valMaxX+20,valMaxY+20); 
	
	return 0; 
}
