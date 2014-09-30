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
typedef ImageSelector<Z2i::Domain, float>::Type Grille;// C'est une Image mais c'est pour éviter les conflits avec le fichier  ReducedMedialAxis.h
typedef ImageContainerBySTLMap<DigitalSetDomain<Z2i::DigitalSet> , DGtal::int64_t> MapImage;
typedef IntervalThresholder<Grille::Value> Binarizer; 
typedef SimpleThresholdForegroundPredicate<Grille> PointPredicate;
typedef DistanceTransformation<Z2i::Space, PointPredicate, Z2i::L2Metric> DTL2;
typedef DistanceTransformation<Z2i::Space, PointPredicate, Z2i::L2Metric> VML2;
typedef HueShadeColorMap<long int, 2> HueTwice;
typedef NotPointPredicate<Z2i::DigitalSet> NotPredicate;
typedef ExactPredicateLpSeparableMetric<Z2i::Space, 2> L2Metric;
typedef ExactPredicateLpPowerSeparableMetric<Z2i::Space,2 > L2PowerMetric;
typedef VoronoiMap<Z2i::Space, NotPredicate, L2Metric > Voronoi2D;
typedef MetricAdjacency<Z2, 1> Adj4;
typedef MetricAdjacency<Z2, 2> Adj8;
typedef DigitalTopology< Adj8, Adj4 > DT8_4;
typedef DigitalTopology< Adj4, Adj8 > DT4_8;
typedef Object<DT4_8, DigitalSet> ObjectType;
typedef Object<DT4_8, DigitalSet> ObjectType48;
typedef typename DigitalSet::ConstIterator DigitalSetConstIterator;
typedef ReverseDistanceTransformation< Grille , L2PowerMetric > RDT;

/* A partir d'un fichier texte, récupère les informations. 
   Ce fichier doit contenir 5 entiers par ligne et les commentaires commencés par le symbole #
   Entreée : le nom d'un fichier
   Sortie : les liste des abscisse, ordonnées et niveau de bruit */
void readData(char *filename,vector<int>&x ,vector<int>&y ,vector<int>&noiseLevel,vector<int>&freeman){
	ifstream data(filename);
	if (data == NULL) exit(1);
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

/* Ces deux fonctions créent un fichier SVG à partir d'une grille. le facteur est optionnel et permet  de diviser toutes les valeurs d'une grille par le même nombre.
   NB : Remplissage ok , mais artéfact sur les contours
   Entrée : une grille
   Sortie : un fichier SVG
*/

void printImageNB(Grille g,char *outputfile,float facteur=1.0){
	int maxX = g.domain().upperBound()[0];
	int maxY = g.domain().upperBound()[1];
	float maxv = 0.0;
	for (int i =0;i<=maxX;i++){
		for (int j =0;j<=maxY;j++){
			Z2i::Point p(i,j);
			g.setValue(p,facteur*g(p));
			if (g(p) > maxv) maxv = g(p);
		}
	}
	Board2D board;
	board.clear();
	Display2DFactory::drawImage<Gray>(board, g , 0.0 ,maxv+0.1);
	board.saveSVG(outputfile);  
}

void printImageColor(Grille g,char *outputfile){
	int maxX = g.domain().upperBound()[0];
	int maxY = g.domain().upperBound()[1];
	float maxv = 0.0;
	for (int i =0;i<=maxX;i++){
		for (int j =0;j<=maxY;j++){
			Z2i::Point p(i,j);
			if (g(p) > maxv) maxv = g(p);
		}
	}
	Board2D board;
	board.clear();
	Display2DFactory::drawImage<HueTwice>(board, g , 0.0 ,maxv+0.1);
	board.saveSVG(outputfile);
}

Z2i::Point regard(int i){
	switch(i){
		case 0: return (Z2i::Point(1,0));
		case 1: return (Z2i::Point(0,1)); 
		case 2: return (Z2i::Point(-1,0)); 
		case 3: return (Z2i::Point(0,-1));
		case 4: return (Z2i::Point(1,1));
		case 5: return (Z2i::Point(-1,1));
		case 6: return (Z2i::Point(-1,-1));
		case 7: return (Z2i::Point(1,-1));
	}
}

/* 
A patir d'une grille, des coordonnées des points et de leur niveau de bruit, extrait les contours intérieur et extérieur
NB : artéfacts sur un des contours gommés lors du remplissage
Entrée : une grille du bord imprécis, la liste des coordonnées des points et leur niveau de bruit
Sortie : Deux listes correspondants aux bords intérieurs et extérieurs
*/
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
	// Recherche des bords des contours imprécis
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
	vector<Z2i::Point> cont1;
	vector<Z2i::Point> cont2;
	for ( DigitalSetConstIterator it = it_begin;it != it_end;++it ) {
		tmp.setValue((*it),0);	
	}
	printImageNB(tmp,"../../../Exemple/contours.svg");
	// recherche d'un point de départ pour le contour interieur
	int a=0;
	int flag = 0;
	Z2i::Point h;
	int nbVoisin;
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
	Z2i::Point depart = Z2i::Point(x[a],maxY-y[a]);
	Z2i::Point currentPT = depart;
	Z2i::Point currentPT2 = depart;
	int stop1 = 0, stop2 = 0;
	int direction = 0;
	int direction2 = 0;
	// extraction des contours
	do{
		//extraction du premier contour
		for (int i=0;i<3;i++){ 
			if (tmp(currentPT+regard((direction+i)%4)) == 0) {
				currentPT += regard((direction+i)%4);
				cont1.push_back(currentPT); 
				direction = (direction +i+3)%4;
				break;
			}  
		}
	}while(currentPT != depart);
	int nb;
	for(int u=0;u< x.size();u++){
		nb = 0;
		if (noise[u] == 1){
			Z2i::Point p(x[u],maxY-y[u]);
			for(int i =0;i<8;i++){
				nb += (image(p+regard(i)) == 0);
			}
			if (nb < 8) cont2.push_back(p);
		}
	}
	// Soustraction des points du premier contour
	for (int u = 0; u< cont1.size() ; u++){
		tmp.setValue(cont1[u],128);
	}
	for(int i = 0;i< maxX;i++){
		for(int j = 0;j< maxY;j++){
			Z2i::Point p(i,j);
			if (tmp(p) == 0) cont2.push_back(p);
		}
	}
	
	cont1.push_back(currentPT);
	 
	contours.push_back(cont1);
	contours.push_back(cont2);
	Display2DFactory::drawImage<Gray>(board, tmp, (unsigned int)0, (unsigned int)129);
	
	board.saveSVG ( "../../../imageBorder.svg" );
	cout << " Fini !!!!" << endl;
	
}

/* Détermine si un point est à l'intérieur d'une forme
   Entrée : une grille représentant le bord de la forme et un point
   Sortie : un booléen valant vrai si le point est à l'interieur de la forme
*/
bool isInside(Grille &image,Z2i::Point p){
	int limX = image.domain().upperBound()[0];
	int limY = image.domain().upperBound()[1];
	int L=0,R=0,U=0,D=0;
	int a,b; 
	//un point est à l'interieur si il coupe le contour un nombre impair de fois
	// Vers la droite
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
	// vers la gauche
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
	// vers le bas
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
	// vers le haut
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

/* Trouve le point de départ pour le remplissage
   Entrée : le contour de la  forme
   Sortie : un point à l'intérieur de a forme 
*/
Z2i::Point findGerme(Grille &image){
	srand (time(NULL));
	int limX = image.domain().upperBound()[0];
	int limY = image.domain().upperBound()[1];
	int x = rand()%limX,y=rand()%limY;
	while (!(isInside(image,Z2i::Point(x,y))) || (image(Z2i::Point(x,y)) == 128 ) ) {
		x = rand()%limX;
		y = rand()%limY;
	}
	return Z2i::Point(x,y);
} 

/* Remplit la forme définie par son contour 
   Entrée : le contour de la forme et un point de départ.
   Sortie : la forme remplie
*/

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

void fill(Grille &image){
	fill(image,findGerme(image));
}  
  

void constructImage(vector< Z2i::Point >& contour,Grille& image){
	for ( Grille::Iterator it = image.begin(), itend = image.end();it != itend; ++it)
    		(*it)=0;
	 
	for (int i=0;i<contour.size();i++){
		image.setValue(contour[i],128);
	}
	//fill(image);
	
	
}
/*Réduit l'épaisseur du contour imprécis en se basant sur le niveau de précision
  Entrée : les coordonnées des points , leur niveau de bruit et le code de freeman 
  Sortie : le contour imprécis réduit
*/
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

/* Calcule la transformée en distance d'une forme
   Entrée : la forme
   Sortie : sa transformée en distance approximative en float
*/

Grille imageDT(Grille& forme,char *outputfile){
	PointPredicate predicate(forme,0);
	int maxX = forme.domain().upperBound()[0];
	int maxY = forme.domain().upperBound()[1];
	Z2i::Point lower(0,0);
	Z2i::Point upper(maxX,maxY);
	Grille result(Z2i::Domain(lower,upper));
	DTL2 dt(forme.domain(), predicate, Z2i::l2Metric);
	
	int maxv = 0;
	for (int i=0;i<maxX;i++){
		for (int j=0;j<maxY;j++){
			Z2i::Point p(i,j);
			result.setValue(p,dt(p));
			if (dt(p) > maxv ) maxv = dt(p);
		}
	}
	Board2D board;
	board.clear();
  	Display2DFactory::drawImage<HueTwice>(board, dt, 0, maxv + 1); 
	board.saveSVG ( outputfile );
	return result;
}

/*Calcule le carré la transformée en distance d'une forme
   Entrée : le diagramme de Voronoi de la forme
   Sortie : sa transformée en distance au carré en entier
*/
Grille imageDT(Voronoi2D& voronoi,char *outputfile){
	Grille result(voronoi.domain());
	int maxv=0;
	for(Voronoi2D::Domain::ConstIterator it = voronoi.domain().begin(),
      	itend = voronoi.domain().end(); it != itend; ++it)
  	{
   		Voronoi2D::Value site = voronoi( *it );   //closest site to (*it)
		if (site != *it){
			Z2i::Point p(site[0]-(*it)[0],site[1]-(*it)[1]);
			result.setValue((*it),p[0]*p[0]+p[1]*p[1]);
			if (result(*it) > maxv) maxv=result(*it); 
		}
		
  	}
	Board2D board;
	board.clear();
  	Display2DFactory::drawImage<HueTwice>(board, result, 0, maxv + 1); 
	board.saveSVG ( outputfile );
	return result;
} 

/* Calcule le diagramme de Voronoi d'une forme
   Entrée : la forme
   Sortie : son diagramme de Voronoi
*/

Voronoi2D voronoiMap(Grille& forme,char *outputfile,char* outputfileCells){
	
	Z2i::DigitalSet set(forme.domain());

	for (int i=0;i<=forme.domain().upperBound()[0];i++){
		for (int j=0;j<=forme.domain().upperBound()[1];j++){
			Z2i::Point p(i,j);
			if (image(p) == 0) set.insertNew(p);
		}
	} 

	NotPredicate notSetPred(set);
	
	L2Metric l2;
	Voronoi2D voronoimap(forme.domain(),notSetPred,l2);
	
	Board2D board;
	board.clear();
	board << forme.domain() << set;
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
	return voronoimap;
} 

/* Calcule la moyenne des transformées en distances sur les formes intérieures et extérieures
   Entrée : 2 transformées en distance 
   Sortie : Leur moyenne
*/ 

void computeDTAverage(Grille & dt1,Grille & dt2,char* outputfile){
	Z2i::Point lower(0,0);
	int maxX = dt1.domain().upperBound()[0];
	int maxY = dt1.domain().upperBound()[1];
	Z2i::Point upper(maxX,maxY);
	Grille DTAv(Z2i::Domain(lower,upper));
	Board2D boardAv;
	for (int i = 0;i<=maxX; i++){
		for (int j = 0;j<=maxY; j++){
			Z2i::Point p(i,j);
			DTAv.setValue(p,(dt1(p)+dt2(p))/2);
		} 
	} 

	float maxAv=0.0;
  	for ( Grille::iterator it = DTAv.begin(), itend = DTAv.end();it != itend; ++it)
    		if ( (*it) > maxAv)  maxAv = (*it);

	Display2DFactory::drawImage<HueTwice>(boardAv, DTAv, 0.0, maxAv +0.1);
	boardAv.saveSVG(outputfile);
}

/* Calcule la différence entre les transformées en distances sur les formes intérieures et extérieures
   Entrée : 2 transformée en distance 
   Sortie : Leur différence
*/ 

void computeDTDiff(Grille & dt1,Grille & dt2,char* outputfile){
	Z2i::Point lower(0,0);
	int maxX = dt1.domain().upperBound()[0];
	int maxY = dt1.domain().upperBound()[1];
	Z2i::Point upper(maxX,maxY);
	Grille DTDiff(Z2i::Domain(lower,upper));
	Board2D boardDiff;
	for (unsigned int i = 0;i<=maxX; i++){
		for (unsigned int j = 0;j<=maxY; j++){
			Z2i::Point p(i,j);
			if (dt1(p)>dt2(p)) DTDiff.setValue(p,dt1(p)-dt2(p));
			else DTDiff.setValue(p,dt2(p)-dt1(p));
		} 
	} 

	float maxDiff=0.0;
  	for ( Grille::iterator it = DTDiff.begin(), itend = DTDiff.end();it != itend; ++it)
    		if ( (*it) > maxDiff)  maxDiff = (*it);

	Display2DFactory::drawImage<Gray>(boardDiff, DTDiff, 0.0, maxDiff +0.1);
	boardDiff.saveSVG(outputfile);
}

/* Calcule la reconstruction
   Entrée : une transformée en distance 
   Sortie : la reconstruction associée
*/

Grille imageRDT(Grille& dt,char* outputfile){
	
	Z2i::L2PowerMetric l2power;

	RDT reverseDT(&dt.domain(),&dt,&l2power);
	
	Grille result (dt.domain());
	for (int i = 0;i<=dt.domain().upperBound()[0];i++){
		for (int j = 0;j<=dt.domain().upperBound()[1];j++){
			Z2i::Point p(i,j);
			if (reverseDT(p) < 0 ) result.setValue(p,2);
		}
	}	 
	Board2D board;
	board.clear();
	Display2DFactory::drawImage<Gray>(board, result ,0, 2);
	board.saveSVG(outputfile);

	return result;
	
}

/* Calcule la reconstruction à partir d'un axe médian
   Entrée : une transformée en distance et un axe médian
   Sortie : la reconstruction associée
*/ 

Grille imageRDTWithAM(Grille& AM,Grille &DT,char* outputfile){
	
	Grille image(AM.domain());
	int maxX = AM.domain().upperBound()[0];
	int maxY = AM.domain().upperBound()[1];
	for (int i = 0;i<=maxX;i++){
		for (int j = 0;j<=maxY;j++){
			Z2i::Point p(i,j);
			image.setValue(p,AM(p)*DT(p));
		}
	}
	Z2i::L2PowerMetric l2power;

	RDT reverseDT(&image.domain(),&image,&l2power);
	
	Grille result (image.domain());
	for (int i = 0;i<=image.domain().upperBound()[0];i++){
		for (int j = 0;j<=image.domain().upperBound()[1];j++){
			Z2i::Point p(i,j);
			if (reverseDT(p) < 0 ) result.setValue(p,2);
		}
	}	 
	Board2D board;
	board.clear();
	Display2DFactory::drawImage<Gray>(board, result ,0, 2);
	board.saveSVG(outputfile);

	return result;
	
} 

/* Calcule l'axe médian à partir d'une forme
   Entrée : une transformée en distance 
   Sortie : la reconstruction associée
*/


Grille imageAM(Grille& dt,char *outputfile){
	Z2i::Domain domain(Z2i::Point(0,0),dt.domain().upperBound());
  	Z2i::Domain domainLarge(Z2i::Point(0,0),dt.domain().upperBound());
	
	/*Grille im(domain);
	for (int i = 0;i<=image.domain().upperBound()[0];i++){
		for (int j = 0;j<=image.domain().upperBound()[1];j++){
			Z2i::Point p(i,j); 
			im.setValue(p,image(p) * image(p));  
		} 
	}*/
	 
	
	Z2i::L2PowerMetric l2power;
  	PowerMap<Grille, Z2i::L2PowerMetric> power(&domainLarge, &dt, &l2power);
	
	ReducedMedialAxis<PowerMap<Grille, Z2i::L2PowerMetric> >::Type  rdma = ReducedMedialAxis< PowerMap<Grille, Z2i::L2PowerMetric> >::getReducedMedialAxisFromPowerMap(power);
	float maxv = 0;
	Grille result(Z2i::Domain(dt.domain().lowerBound(),dt.domain().upperBound()));
	for (int i = 0;i<=result.domain().upperBound()[0];i++){
		for (int j = 0;j<=result.domain().upperBound()[1];j++){
			Z2i::Point p(i,j);
			result.setValue(p,rdma(p));
			if (rdma(p) > maxv) maxv = rdma(p);  
		} 
	}

	Board2D board;
	board.clear();
	Display2DFactory::drawImage<Gray>(board, rdma , (unsigned int)0, maxv +1);
	board.saveSVG(outputfile);

	return result;
} 

/* Opère un amincissement homotopique
   Entrée : une forme
   Sortie : la forme amincie
*/

Grille homotopicThinning(Grille &forme,char *outputfile){
	PointPredicate predicate(forme,0);
	Z2i::Point lower(0,0);
	Z2i::Point upper(forme.domain().upperBound()[0],forme.domain().upperBound()[1]);
	DTL2 dt(forme.domain(), predicate, Z2i::l2Metric);
	Domain domain( lower, upper );
	DigitalSet set( domain );
	for (int i = 0; i<=forme.domain().upperBound()[0] ; i++){
		for (int j = 0; j<=forme.domain().upperBound()[1] ; j++){
			Z2i::Point p(i,j); 
			if (forme(p) != 0) set.insertNew(p);   
		}
	}
	ObjectType48 shape(Z2i::dt4_8 , set);
	int nb_simple=0; 
  	int layer = 1;
  	std::queue<DigitalSet::Iterator> Q;
	do {
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

	return imageRes;
}

#define PI 3.14159265
/* Determine l'angle dont l'origine est le point de la forme et les deux autres points sont les points les plus proches, un sur le bord intérieur, lautre sur le bord extérieur / utilisation du théorème d'Al-Kashi
   Entrée : Les deux formes accompagné de leur diagrammede Voronoi respectif
   Sortie : Carte des angles
*/
Grille imageAngle(Grille &i1,Grille &i2,Voronoi2D &v1,Voronoi2D &v2,char *outputfile){
	Grille result(i1.domain());
	int max = 0;
	for (int i = 0; i <= i1.domain().upperBound()[0] ; i++){
		for (int j = 0; j <= i1.domain().upperBound()[1] ; j++){
			Z2i::Point p(i,j);
			if (i1(p) != 0 && i2(p) != 0){
				Z2i::Point q1 = v1(p);
				Z2i::Point q2 = v2(p);
				int l = (abs(q1[0]- q2[0]) * abs(q1[0]- q2[0])) + (abs(q1[1]- q2[1]) * abs(q1[1]- q2[1]));
				int l1 = (abs(q1[0]- p[0]) * abs(q1[0]- p[0])) + (abs(q1[1]- p[1]) * abs(q1[1]- p[1]));
				int l2 = (abs(p[0]- q2[0]) * abs(p[0]- q2[0])) + (abs(p[1]- q2[1]) * abs(p[1]- q2[1]));
				result.setValue(p,acos( (l-l1-l2)/(-2*sqrt(l1)*sqrt(l2))) * 180/PI );
				if (result(p) > max ) max = result(p);
			} 
		} 
	}
	
	Board2D board;
	board.clear();
	Display2DFactory::drawImage<Gray>(board, result , (unsigned int) 0 ,max+1);
	board.saveSVG(outputfile);  
	return result;
}

/* Calcul de la distance entre 2 points
   Entrée : 2 points
   Sortie : la distance entre les deux points 
*/
double dist(Z2i::Point p1,Z2i::Point p2){
	return sqrt((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1]));
}

int distSquare(Z2i::Point p1,Z2i::Point p2){
	return ((p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1]));
}

Grille intersection(Grille &contour,Z2i::Point p,int r,char *outputfile= ""){
	Grille result(contour.domain());
	int maxX = contour.domain().upperBound()[0];
	int maxY = contour.domain().upperBound()[1];
	for (int i = 0;i<=maxX;i++){
		for (int j =0;j<= maxY;j++){
			Z2i::Point q(i,j);
			if (dist(p,q) <= sqrt(r) && contour(q) != 0){
				result.setValue(q,1);
			}  
		}
	}
	if (outputfile != ""){
		Board2D board;
		board.clear();
		Display2DFactory::drawImage<Gray>(board, result , 0 ,1);
		board.saveSVG(outputfile);
	}  
	return result;
}

/*Compte le nombre de composantes connexes
  Entrées : une grille avec des composantes connexes
  Sortie : le nombre de composantes connexes 
*/

int nbComposantesConnexes(Grille &image){
	DigitalSet image_set(image.domain());
	int maxX = image.domain().upperBound()[0];
	int maxY = image.domain().upperBound()[1];
	for (unsigned int i = 0;i<maxX; i++){
		for (unsigned int j = 0;j<maxY; j++){
			Z2i::Point p(i,j);
			if (image(p) != 0) image_set.insertNew(p); 
    		}
	}

	Adj4 adj4;
	Adj8 adj8;
	DT4_8 dt4_8 ( adj4, adj8, JORDAN_DT );

	ObjectType im ( dt4_8, image_set );
	vector<ObjectType> obs;
	obs.push_back(im);
  	back_insert_iterator< vector< ObjectType > > insert( obs );
  	unsigned int nbc0 = im.border().writeComponents( insert );
	return nbc0;
}



void parabolloide(Grille &g,Z2i::Point &point,float dt,float e = 0.0){
	int maxX = g.domain().upperBound()[0];
	int maxY = g.domain().upperBound()[1];
	for (int i = 0 ; i <= maxX ; i++ ){
		for (int j = 0 ; j <= maxY ; j++ ){
			Z2i::Point p(i,j);
			if (e == 0.0) g.setValue(p,distSquare(p,point)-dt);
			else g.setValue(p,distSquare(p,point)+e);
		}
	}
}

/* Calcul le négatif d'une image : un pixel blanc sera noir et inversement
   Entrée : une image 
   Sortie : son négatif
*/

Grille negatif(Grille &g){
	int maxX = g.domain().upperBound()[0];
	int maxY = g.domain().upperBound()[1];
	Grille result(g.domain());
	for (int i = 0 ; i <= maxX ; i++ ){
		for (int j = 0 ; j <= maxY ; j++ ){
			Z2i::Point p(i,j);
			if (g(p) == 0) result.setValue(p,1);
			else  result.setValue(p,0);
		}
	}
	return result;
}

/////////////////////////////////////////////Reconstruction : calcul des alphas ////////////////////////////////////

/* Determine les alphas min pour une partie de la zone hybride ( union des boules de plus grand rayon)
   Entrée : les 2 transformées en distance 
   Sortie : les alphas pour chaque point de la zone hybride ( voir page 14 )
*/

Grille imageAlphaMin(Grille &dt1,Grille &dt2,Grille &dtinv,char *outputfile){
	Z2i::L2PowerMetric l2power;
	Z2i::Domain dom = dt1.domain();
	Grille n1 = negatif(dt1);
	Grille n2 = negatif(dt2);
    Grille ninv = negatif(dtinv);
    PowerMap<Grille, Z2i::L2PowerMetric> powerinv(&dom, &ninv, &l2power);
  	PowerMap<Grille, Z2i::L2PowerMetric> power1(&dom, &n1, &l2power);
	PowerMap<Grille, Z2i::L2PowerMetric> power2(&dom, &n2, &l2power);
    Grille result(dt1.domain());
	int maxX = dt1.domain().upperBound()[0];
	int maxY = dt1.domain().upperBound()[1];
	Grille fm(dt1.domain()),fM(dt1.domain());
    float rm,rM,d1,d2,dinv;
	float x,e;
	//initialisation
	for (int i = 0 ; i <= maxX ; i++ ){
		for (int j = 0 ; j <= maxY ; j++ ){
			Z2i::Point p(i,j);
            result.setValue(p,1.4);
		}
	}
	for (int i = 0 ; i <= maxX ; i++ ){
		for (int j = 0 ; j <= maxY ; j++ ){
			Z2i::Point p(i,j);
            if (dt1(p) != 0 && dt2(p) != 0 ) {
                Z2i::Point powerp1 = powerinv(p);
                Z2i::Point powerp2 = power2(p);
                d1 = dist(p,powerp1);
                d2 = dist(p,powerp2);
                rm = min(d1,d2);
                rM = max(d1,d2);
                parabolloide(fm,p,dt1(p));
                parabolloide(fM,p,dt2(p));
                for (int a = 0 ; a <= maxX ; a++ ){
                    for (int b = 0 ; b <= maxY ; b++ ){
                        Z2i::Point q(a,b);
                        if (fm(q) <= 0 ) result.setValue(q,0.0);
                        else{
                            if (fM(q) < 0 ){
				//calcul de la puissance
                                x = (float)((float)(dist(p,q)-rm)/(float)(rM-rm));
                                if (x>=0.0 && x<=1.0) result.setValue(q,(float)(min(x,result(q))));
                            }
                        }
                    }
                }
			}
		}
	}
	for (int i = 0 ; i <= maxX ; i++ ){
		for (int j = 0 ; j <= maxY ; j++ ){
			Z2i::Point p(i,j);
            result.setValue(p,100 * result(p));
		}
	}
	Board2D board;
	board.clear();
    Display2DFactory::drawImage<Gray>(board, result ,0,140);
	board.saveSVG(outputfile);  
	return result;
}


/* Determine les alpha min pour l'intégralité des points de la zone hybride
   Entrée : les 2 transformées en distance 
   Sortie : les alphas pour chaque point de la zone hybride
*/


Grille imageAlphaMinDynamique(Grille &dtInt,Grille &dtExt,char *outputfile){
	Grille result(dtInt.domain());
	int maxX = dtInt.domain().upperBound()[0];
	int maxY = dtInt.domain().upperBound()[1];
	for (int i = 0 ; i <= maxX ; i++ ){
		for (int j = 0 ; j <= maxY ; j++ ){
			Z2i::Point p(i,j);
            result.setValue(p,5.0);
		}
	}
	int oldN=0,newN=1;
	float rm,rM,d,x;
	while (oldN != newN){ // tant qu'on augmente le nombre de points, on continue
		oldN = newN; 
		for (int i = 0 ; i <= maxX ; i++ ){
			for (int j = 0 ; j <= maxY ; j++ ){
				Z2i::Point p(i,j);
		   		if ( ( (dtInt(p) != 0) && (dtExt(p) != 0) )  || ( result(p) < 5.0 ) ) {
					rm = dtInt(p);
					rM = dtExt(p);
					for (int a = 0 ; a <= maxX ; a++ ){
		            			for (int b = 0 ; b <= maxY ; b++ ){
		                			Z2i::Point q(a,b);
							d = distSquare(p,q);
							x = float(float(d-rm)/float(rM - rm));
							if (d < rm) result.setValue(p,0.0);
							else {
								if (d < rM){
									if (dtInt(p) == 0.0 ) result.setValue(q,float(min(result(q),result(p)+x)));
									else result.setValue(q,float(min(result(q),x))); 
								}
							}
						}
					}	
				}
			}
		
		}
	newN = nbPoints(result,5.0);
	cout << oldN << " " << newN << endl;
	}
	float maxv = 0;
	for (int i = 0 ; i <= maxX ; i++ ){
		for (int j = 0 ; j <= maxY ; j++ ){
			Z2i::Point p(i,j);
            		if ((maxv < result(p)) && (result(p) < 5.0) ) maxv = result(p);
		}
	}
	cout << maxv << endl;
	maxv += 0.1; 
	for (int i = 0 ; i <= maxX ; i++ ){
		for (int j = 0 ; j <= maxY ; j++ ){
			Z2i::Point p(i,j);
			if (result(p) == 5.0) result.setValue(p,100);
            		else result.setValue(p,100 * float( float( result(p) ) / float( maxv )));
		}
	}
	Board2D board;
	board.clear();
    	Display2DFactory::drawImage<Gray>(board, result ,0,100);
	board.saveSVG(outputfile);  
	return result;
}

/* Détermine les alphas avec vitesse de grossissement de boules constante
   Entrées : 2 transformées en diastance
   Sortie : La carte des Alphas */

Grille imageAlphaMinEqual(Grille &dtInt,Grille &dtExt,char *outputfile){
	Z2i::L2PowerMetric l2power;
	Z2i::Domain dom = dtInt.domain();
	PowerMap<Grille, Z2i::L2PowerMetric> powerInt(&dom, &dtInt, &l2power);
  	PowerMap<Grille, Z2i::L2PowerMetric> powerExt(&dom, &dtExt, &l2power);
	Grille result(dtInt.domain());
	int maxX = dtInt.domain().upperBound()[0];
	int maxY = dtInt.domain().upperBound()[1];
	int delta=0;
	float r;
	for (int i = 0; i <= maxX;i++){
		for(int j = 0; j<= maxY; j++){
			Z2i::Point p(i,j);
			if (dtExt(p)-dtInt(p) > delta ) delta = dtExt(p)-dtInt(p);
			result.setValue(p,130);
		}
	}
	for (int i = 0; i <= maxX;i++){
		for(int j = 0; j<= maxY; j++){
			Z2i::Point p(i,j);
			if (dtExt(p) > 0 ){
				for (int a = 0; a <= maxX;a++){
					for(int b = 0; b<= maxY; b++){
						Z2i::Point q(a,b);
						if (dtInt(q) > 0){
							r = distSquare(p,q) - dtInt(q);
							if (r >= 0) {
								result.setValue(p,min(100*(float)(float(r)/float(delta)),(float)result(p)));
							}
							else result.setValue(p,0.0);
						}
					}
				}	
			}
		}
	}
	Board2D board;
	board.clear();
    	Display2DFactory::drawImage<Gray>(board, result ,0,130);
	board.saveSVG(outputfile);  
	return result;
}



/* Determine les alphas min à partir d'un axe médian
   Entrée : deux transformées en distance et l'axe médian
   Sortie : les Alphas min ( voir page 14 ) 
*/ 

Grille imageAlphaMinWithAM(Grille AM,Grille &dt1,Grille &dt2,Grille &dtinv,char *outputfile){
	Z2i::L2PowerMetric l2power;
	Z2i::Domain dom = dt1.domain();
	Grille n1 = negatif(dt1);
	Grille n2 = negatif(dt2);
    	Grille ninv = negatif(dtinv);
    	PowerMap<Grille, Z2i::L2PowerMetric> powerinv(&dom, &ninv, &l2power);
  	PowerMap<Grille, Z2i::L2PowerMetric> power1(&dom, &n1, &l2power);
	PowerMap<Grille, Z2i::L2PowerMetric> power2(&dom, &n2, &l2power);
    	Grille result(dt1.domain());
	int maxX = dt1.domain().upperBound()[0];
	int maxY = dt1.domain().upperBound()[1];
	Grille fm(dt1.domain()),fM(dt1.domain());
    	float rm,rM,d1,d2,dinv;
	float x,e;
	for (int i = 0 ; i <= maxX ; i++ ){
		for (int j = 0 ; j <= maxY ; j++ ){
			Z2i::Point p(i,j);
            result.setValue(p,1.4);
		}
	}
	for (int i = 0 ; i <= maxX ; i++ ){
		for (int j = 0 ; j <= maxY ; j++ ){
			Z2i::Point p(i,j);
            if (AM(p)*dt1(p) != 0 || AM(p)*dt2(p) != 0 ) {
                Z2i::Point powerp1 = power1(p);
                Z2i::Point powerp2 = power2(p);
                d1 = dist(p,powerp1);
                d2 = dist(p,powerp2);
                rm = min(d1,d2);
                rM = max(d1,d2);
                parabolloide(fm,p,dt1(p));
                parabolloide(fM,p,dt2(p));
                for (int a = 0 ; a <= maxX ; a++ ){
                    for (int b = 0 ; b <= maxY ; b++ ){
                        Z2i::Point q(a,b);
                        if (fm(q) <= 0 ) result.setValue(q,0.0);
                        else{
                            if (fM(q) < 0 ){
                                //x = (float)(fm(q)/(fm(q)-fM(q)));
				x = (float)((dist(p,q)-rm)/(rM-rm));
                                if (x>=0.0 && x<=1.0) result.setValue(q,(float)(min(x,result(q))));
                            }
                        }
                    }
                }
			}
		}
	}
	for (int i = 0 ; i <= maxX ; i++ ){
		for (int j = 0 ; j <= maxY ; j++ ){
			Z2i::Point p(i,j);
            result.setValue(p,100 * result(p));
		}
	}
	Board2D board;
	board.clear();
    Display2DFactory::drawImage<Gray>(board, result ,0,140);
	board.saveSVG(outputfile);  
	return result;
}


/* Compte le nombre de points d'une image
   Entrée : une image 
   Sortie : le nombre de point qui la compose
*/
int nbPoints(Grille &image,float val = 0.0){
	int cpt  =0;
	int maxX = image.domain().upperBound()[0];
	int maxY = image.domain().upperBound()[1];
	for (int i =0;i<= maxX;i++){
		for (int j =0;j<= maxY;j++){
			Z2i::Point p(i,j);
			cpt += (image(p) != val);
		}
    }
	return cpt;
}
/* Compte le nombre de zone qui contient au moins un point
   Entrée : une image et un rayon
   Sortie : le nombre de zones 
*/
int nbZones(Grille &image,int r,float val = 0.0){
	int cpt  =0;
	int X,Y,R=2*r+1,flag;
	int maxX = image.domain().upperBound()[0];
	int maxY = image.domain().upperBound()[1];
	for (int i =0;i< maxX/R;i++){
		for (int j =0;j< maxY/R;j++){
			X=r+i*R;
			Y=r+j*R;
			flag =0;
			for (int a = 0;a<R;a++){
				for (int b = 0;b<R;b++){
					Z2i::Point p(X+a-r,Y+b-r);
					if (image(p) != val) {cpt++;flag=1;break;}
				}
				if (flag) break;
			}
		}
    }
	return cpt;
}

/* Determine les alphas ( rapport entre la distance à la forme intérieure et la somme de celles à la forme intérieure et extérieure)
   Entrées : 2 transformées en distance + transformée en distance du complémentaire
   Sortie : Les alphas 

*/

Grille alphaFromIntGlobal(Grille &dtExt,Grille &dtIntInv,char *outputfile){
	int maxX = dtExt.domain().upperBound()[0];
	int maxY = dtExt.domain().upperBound()[1];
	Grille result(dtExt.domain());
	int maxv = 0;
	for (int i =0;i<= maxX;i++){
		for (int j =0;j<= maxY;j++){
			Z2i::Point p(i,j);
			if (dtExt(p) != 0 && dtIntInv(p) != 0){
				if (dtIntInv(p) > maxv) maxv = dtIntInv(p); 
			} 
		}
	}
	for (int i =0;i<= maxX;i++){
		for (int j =0;j<= maxY;j++){
			Z2i::Point p(i,j);
			if (dtExt(p) != 0 && dtIntInv(p) != 0){
				result.setValue(p,100*dtIntInv(p)/maxv); 
			} 
		}
	}
	
	Board2D b;
	b.clear();
	Display2DFactory::drawImage<Gray>(b, result ,0,100);
	b.saveSVG(outputfile);
	return result;
}

/* Calcule l'intégralité des alpha-formes 
   Entrée : 2 transformées en distance, un pas d'échantillonage
   Sortie : Les alphas-formes
*/

void allFAlpha(Grille &dt1,Grille &dt2,float step,char *outputfileAM,char *outputfileRDT,int r=0){
	int n = 1.0/step;
	Grille AM(dt1.domain());
	Grille RDT(dt1.domain());
	Grille resultAM(dt1.domain());
	Grille resultRDT(dt1.domain());
	int maxX = dt1.domain().upperBound()[0];
	int maxY = dt1.domain().upperBound()[1];
	int X,Y,R=2*r+1,flag;
	for (int s =0;s<=n;s++){
		AM = imageAlphaAM(dt1,dt2,s*step,"");
		RDT = imageRDT(AM,"");
		//if (s == 0) printImageNB(RDT,"../../../exemple.svg");
		for (int i = 0;i< maxX/R;i++){
			for (int j = 0;j< maxY/R;j++){
				X=r+i*R;
				Y=r+j*R;
				flag = 0;
				for (int a = 0;a<R;a++){
					for (int b = 0;b<R;b++){
						Z2i::Point p(X+a-r,Y+b-r);
						if (AM(p) != 0) {
						flag =1;
							for (int c = 0;c<R;c++){
								for (int d = 0;d<R;d++){
									Z2i::Point pp(X+c-r,Y+d-r);
									resultAM.setValue(pp,resultAM(pp)+1);
								}
							}
						break;
						}
					}
				if (flag) break;
				}
			}
		}
	}
	Board2D boardAM,boardRDT;
	boardAM.clear();
	boardRDT.clear();
	Display2DFactory::drawImage<Gray>(boardAM, resultAM , 0 ,n+1);
	Display2DFactory::drawImage<Gray>(boardRDT, resultRDT , 0 ,n+1);
	boardAM.saveSVG(outputfileAM); 
	boardRDT.saveSVG(outputfileRDT);
}



void preuve(Grille &contour,Grille &dt1,Grille &dt2,Grille &ma){
	int maxX = dt1.domain().upperBound()[0];
	int maxY = dt1.domain().upperBound()[1];
	int rm;
	int rM;
	Grille lesCC(dt1.domain());
	Grille exemple(dt1.domain());
	int maxv = 0;
	for (int a =0;a<= maxX;a++){
		for (int b =0;b<= maxY;b++){
			Z2i::Point x(a,b);
			if(contour(x) != 0) exemple.setValue(x,2);
		}
	}
	int flag = 0;
	for (int i =0;i<= maxX;i++){
		for (int j =0;j<= maxY;j++){
			Z2i::Point p(i,j);
			if (ma(p) != 0 ){
				cout << p << ":" << endl;
				rm = min(dt1(p),dt2(p));
				rM = max(dt1(p),dt2(p));
				for (int r = rm;r<=rM;r++){
					lesCC = intersection(contour,p,r*r);
					cout << nbPoints(lesCC) << " " << nbComposantesConnexes(lesCC)<< endl;
					if (nbComposantesConnexes(lesCC) == 1 && !flag){
						//flag = 1;
						Grille boule(dt1.domain());
						for (int a =0;a<= maxX;a++){
							for (int b =0;b<= maxY;b++){
								Z2i::Point x(a,b);
								if(distSquare(p,x) <= ma(p) ) boule.setValue(x,1);
							}
						}
						for (int a =0;a<= maxX;a++){
							for (int b =0;b<= maxY;b++){
								Z2i::Point x(a,b);
								int y=0;
								if(boule(x) == 1){
									for (int d = 0; d<8;d++){
										if (boule(x+regard(d)) == 1) y++;
									} 
									if (y < 8) {
										if (exemple(x) == 2 || exemple(x) == 3 ) exemple.setValue(x,3);
										else exemple.setValue(x,1);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	Board2D b;
	b.clear();
	Display2DFactory::drawImage<Gray>(b, exemple , (unsigned int) 0 ,3);
	b.saveSVG("../../../exemple.svg");
}

/* Distingue la DT de la forme intérieure de celle de la forme extérieure
   Entrées : 2 transformées en distance 
   Sortie : Les 2 transformées en distance distinguées
*/

void distingueIntExt(Grille &dt1,Grille &dt2,Grille &dtInt,Grille &dtExt){
	int maxX = dt1.domain().upperBound()[0];
	int maxY = dt1.domain().upperBound()[1];
	int flag = 0;
	for (int i =0;i<= maxX;i++){
		for (int j =0;j<= maxY;j++){
			Z2i::Point p(i,j);
			if ((dt1(p) == 0) && (dt2(p) != 0)){
				dtInt = dt1;
				dtExt = dt2;
				flag = 1;
				break;
			}
			if ((dt1(p) != 0) && (dt2(p) == 0)){
				dtInt = dt2;
				dtExt = dt1;
				flag = 1;
				break;
			}
		}
		if (flag) break;
	}	
}

/* Calcule la reconstruction dynamique 
   Entrées : les 2 tranformées en distance 
   Sortie : La reconstruction dynamique ( voir page 16 )
*/

Grille RDTDynamique(Grille dtInt,Grille dtExt){
	int maxX = dtInt.domain().upperBound()[0];
	int maxY = dtInt.domain().upperBound()[1];
	Grille result(dtInt.domain());
	Grille dttmp = dtInt;
	Grille oldrdt(dtInt.domain());
	Grille rdttmp = imageRDT(dttmp,"");
	int cpt=1;
	do{
		for (int i =0;i<= maxX;i++){
			for (int j =0;j<= maxY;j++){
				Z2i::Point p(i,j);
				if (rdttmp(p) != 0){
					dttmp.setValue(p,dtExt(p));
					if (result(p) == 0)result.setValue(p,cpt);
				}
			}
		}
		cpt++;
		oldrdt = rdttmp;
		rdttmp = imageRDT(dttmp,"");
	}while(nbPoints(oldrdt) != nbPoints(rdttmp)); 

	Board2D b;
	b.clear();
	Display2DFactory::drawImage<Gray>(b, result , (unsigned int) 0 ,cpt);
	b.saveSVG("../../../RDT_Dynamique.svg");
	
	return result;
}

/* Calcule la reconstruction dynamique à partir d'un axe médian 
   Entrées : les 2 tranformées en distance et un axe médian
   Sortie : La reconstruction dynamique ( voir page 16 )
*/

Grille RDTDynamiqueWithAM(Grille AM,Grille dtInt,Grille dtExt,char *outputfile){
	int maxX = dtInt.domain().upperBound()[0];
	int maxY = dtInt.domain().upperBound()[1];
	Grille result(dtInt.domain());
	Grille dttmp(dtInt.domain());
	for (int i =0;i<= maxX;i++){
			for (int j =0;j<= maxY;j++){
				Z2i::Point p(i,j);
				dttmp.setValue(p,dtInt(p)*AM(p));
			}
	}
	Grille oldrdt(dtInt.domain());
	Grille rdttmp = imageRDT(dttmp,"");
	int cpt=1;
	do{
		for (int i =0;i<= maxX;i++){
			for (int j =0;j<= maxY;j++){
				Z2i::Point p(i,j);
				if (rdttmp(p) != 0){
					dttmp.setValue(p,AM(p)*dtExt(p));
					if (result(p) == 0)result.setValue(p,cpt);
				}
			}
		}
		cpt++;
		oldrdt = rdttmp;
		rdttmp = imageRDT(dttmp,"");
	}while(nbPoints(oldrdt) != nbPoints(rdttmp)); 

	Board2D b;
	b.clear();
	Display2DFactory::drawImage<Gray>(b, result , (unsigned int) 0 ,cpt);
	b.saveSVG(outputfile);
	
	return result;
}

/*Calcule les axes médians à partir d'une reconstruction dynamique
  Entrée : une reconstruction dynamique
  Sortie : Les différents axes médians
*/ 

Grille MAfromRDTDym(Grille &rdtDym,char *outputfile){
	Grille forme(rdtDym.domain());
	int maxX = rdtDym.domain().upperBound()[0];
	int maxY = rdtDym.domain().upperBound()[1];
	int nb=0;
	for (int i=0;i<=maxX;i++){
		for (int j=0;j<=maxY;j++){
			Z2i::Point p(i,j);
			if (rdtDym(p) > nb) nb = rdtDym(p);  
		}
	}
	
	Grille DT(rdtDym.domain());
	Grille AM(rdtDym.domain());
	Grille result(rdtDym.domain());

	for (int a=1;a<=nb;a++) {
		for (int i=0;i<=maxX;i++){
			for (int j=0;j<=maxY;j++){
				Z2i::Point p(i,j);
				if (rdtDym(p) <= a && rdtDym(p) != 0  ) forme.setValue(p,1);	
			}
		}
		Voronoi2D VM = voronoiMap(forme,"","");
		DT = imageDT(VM,"");
		AM = imageAM(DT,"");
		for (int i=0;i<=maxX;i++){
			for (int j=0;j<=maxY;j++){
				Z2i::Point p(i,j);
				if (AM(p) > 0 ) result.setValue(p,result(p) + 1);	
			}
		} 
	}
	
	Board2D b;
	b.clear();
	Display2DFactory::drawImage<Gray>(b, result , (unsigned int) 0 ,nb);
	b.saveSVG(outputfile);
	return result;			
}



Z2i::Point findFirstCC(Grille &CC){
	srand (time(NULL));
	int limX = CC.domain().upperBound()[0];
	int limY = CC.domain().upperBound()[1];
	int x = rand()%limX,y=rand()%limY;
	while (CC(Z2i::Point(x,y)) != 100 ) {
		x = rand()%limX;
		y = rand()%limY;
	}
	return Z2i::Point(x,y);
}

void fillCC(Z2i::Point p,Grille &CC, int x){
	int limX = CC.domain().upperBound()[0];
	int limY = CC.domain().upperBound()[1];
	for (int i =0; i<4;i++){
		Z2i::Point pt = p+regard(i);
		if (CC(pt) == 100 && pt[0]<=limX && pt[1]<=limY) {
			CC.setValue(pt,x);
			fillCC(pt,CC,x);
		} 
	}
}

void fillCC(Grille &CC,int x){
	fillCC(findFirstCC(CC),CC,x);
}

/* Determine les alphas ( méthode des composantes connexes )
   Entrées : 2 transformées en distance + transformée en distance du complémentaire
   Sortie : Les alphas

*/

Grille alphaFromIntLocal(Grille &dtExt,Grille &dtIntInv,char *outputfile){
	int maxX = dtExt.domain().upperBound()[0];
	int maxY = dtExt.domain().upperBound()[1];
	Grille result(dtExt.domain());
	Grille CC(dtExt.domain());
	for (int i =0;i<= maxX;i++){
		for (int j =0;j<= maxY;j++){
			Z2i::Point p(i,j);
			if (dtExt(p) != 0 && dtIntInv(p) != 0){
				CC.setValue(p,100); 
			} 
		}
	}

	int nbCC = nbComposantesConnexes(CC) ;
	
	for (int i = 1;i<=nbCC;i++){
		fillCC(CC,i);
	}
	
	float maxv[nbCC+1];

	for (int i =0;i<= maxX;i++){
		for (int j =0;j<= maxY;j++){
			Z2i::Point p(i,j);
			if (dtExt(p) != 0){
				if (dtIntInv(p) > maxv[(int)CC(p)]) maxv[(int)CC(p)] = dtIntInv(p);
			}
		}
	}

	for (int i =0;i<= maxX;i++){
		for (int j =0;j<= maxY;j++){
			Z2i::Point p(i,j);
			if (dtExt(p) != 0){
				result.setValue(p,100*dtIntInv(p)/maxv[(int)CC(p)]);
			}
		}
	}
	
	Board2D b;
	b.clear();
	Display2DFactory::drawImage<Gray>(b,result,0,100);
	b.saveSVG(outputfile);
	return result;
}


/* Calcul le alpha max pour chaque étape d'une reconstruction dynamique
   Entrée : une reconstruction dynamique
   Sortie : le vecteur des alphas max
*/ 
void alphas(vector<float> &vect,Grille &rdtDym,Grille &alpha){
	int maxX = rdtDym.domain().upperBound()[0];
	int maxY = rdtDym.domain().upperBound()[1];
	int nb=0;
	for (int i=0;i<=maxX;i++){
		for (int j=0;j<=maxY;j++){
			Z2i::Point p(i,j);
			if (rdtDym(p) > nb) nb = rdtDym(p);  
		}
	}
	for (int i=1;i<=nb;i++) vect.push_back(0.0);

	for (int i=0;i<=maxX;i++){
		for (int j=0;j<=maxY;j++){
			Z2i::Point p(i,j);
			if (rdtDym(p) > 0){
				if (alpha(p) > vect[rdtDym(p)-1]) vect[rdtDym(p)-1] = alpha(p);
			}  
		}
	} 
}


///////////////////////////////// Calcul de l'axe médian ////////////////////////////////

/* Calcul de l'axe médian pour une union d'une boule de chaque hyperboule
   Entrée : 2 transformées en distance et un réel alpha compris entre 0 et 1
   Sortie : L'axe médian associée
*/

Grille imageAlphaAM (Grille &dt1,Grille &dt2,float alpha,char *outputfile){
	assert(alpha<=1.0);
	assert(alpha>=0.0);
	Grille input(dt1.domain());
	int maxX = dt1.domain().upperBound()[0];
	int maxY = dt1.domain().upperBound()[1];
	int rm,rM; 
	// Formule pour le calcul de l'axe médian
	for (int i = 0 ; i <= maxX ; i++ ){
		for (int j = 0 ; j <= maxY ; j++ ){
			Z2i::Point p(i,j);
			if (dt1(p) > 0 && dt2(p) >0 ){
				rm = min(dt1(p),dt2(p));
				rM = max(dt1(p),dt2(p));
				input.setValue(p,alpha*rM + (1.0-alpha) * rm);
			}
		}
	}
	Grille result(dt1.domain());
	int flag;
	int maxv=0;
	for (int ia = 0;ia<= maxX;ia++){
		for (int ja = 0;ja<= maxY;ja++){
			Z2i::Point p(ia,ja);
			if (input(p) != 0){
				flag = 0;
				for (int ib = 0;ib<= maxX;ib++){
					for (int jb = 0;jb<= maxY;jb++){
						Z2i::Point q(ib,jb);
						if (input(q) != 0 && p != q){
							//int distanceSquare = distSquare(p,q);
							int distance = dist(p,q);
							if ( distance >= input(q) ){
							}
							else {
								if ( distance+ input(p) <= input(q) ) {
									// la boule p est contenu dans la boule q
									flag = 1;
									break;
								}
							}
						}
					}
					if (flag) break;
				}
				if (!flag){
					result.setValue(p,input(p)*input(p));
					if (result(p) >= maxv) maxv = result(p); 
				}
			}
		}
	}

	Board2D board;
	board.clear();
	Display2DFactory::drawImage<Gray>(board, result , 0.0 ,maxv+0.1);
	board.saveSVG(outputfile);  
	return result;
}


/* Calcul de l'axe médian pour une union d'une boule de chaque hyperboule
   Entrée : 2 transformées en distance et deux réels alpha et béta compris entre 0 et 1
   Sortie : L'axe médian associé 
*/

Grille imageAlphaBetaAM (Grille &dt1,Grille &dt2,float alpha,float beta,char *outputfile){
	assert(alpha<=1.0);
	assert(alpha>=0.0);
	assert(beta<=1.0);
	assert(beta>=0.0);
	Grille inputAlpha(dt1.domain());
	Grille inputBeta(dt1.domain());
	int maxX = dt1.domain().upperBound()[0];
	int maxY = dt1.domain().upperBound()[1];
	int rm,rM; 
	for (int i = 0 ; i <= maxX ; i++ ){
		for (int j = 0 ; j <= maxY ; j++ ){
			Z2i::Point p(i,j);
			if (dt1(p) > 0 && dt2(p) >0 ){
				rm = min(dt1(p),dt2(p));
				rM = max(dt1(p),dt2(p));
				inputAlpha.setValue(p,alpha*rM + (1.0-alpha) * rm);
				inputBeta.setValue(p,beta*rM + (1.0-beta) * rm);
			}
		}
	}
	Grille result(dt1.domain());
	int flag;
	float maxv=0;
	for (int ia = 0;ia<= maxX;ia++){
		for (int ja = 0;ja<= maxY;ja++){
			Z2i::Point p(ia,ja);
			if (inputAlpha(p) != 0){
				flag = 0;
				for (int ib = 0;ib<= maxX;ib++){
					for (int jb = 0;jb<= maxY;jb++){
						Z2i::Point q(ib,jb);
						if (inputBeta(q) != 0 && p != q){
							int distance = (int)dist(p,q);
							if ( distance >= inputBeta(q) ){NULL;}
							else {
								if (distance + inputAlpha(p) <= inputBeta(q)){
									// la boule p est contenu dans la boule q
									flag = 1;
									break;
								}
							}
						}
					}
					if (flag) break;
				}
				if (!flag){
					result.setValue(p,inputAlpha(p) * inputAlpha(p));
					if (result(p) >= maxv) maxv = result(p); 
				}
			}
		}
	}

	Board2D board;
	board.clear();
	Display2DFactory::drawImage<Gray>(board, result , 0.0 ,maxv+0.1);
	board.saveSVG(outputfile);  
	return result;
}


/* Détermine l'axe médian par simplification des hyperboules par méthode des rayons caractéristiques ( voir page 19 )
   Entrée : Les deux transformées en distance 
   Sortie : l'axe médian 
*/
Grille imageHAM (Grille &dt1,Grille &dt2,char *outputfile){// axe médian avec hyperboules
	Grille result(dt1.domain());
	int maxX = dt1.domain().upperBound()[0];
	int maxY = dt1.domain().upperBound()[1];
	int flag;
	int rmp,rMp,rmq,rMq,maxv=0;
	float rcp,rcq,rcmax,taup,tauq;
	float distance;
	for (int i = 0 ; i <= maxX ; i++ ){
		for (int j = 0 ; j <= maxY ; j++ ){
			Z2i::Point p(i,j);
			if ((dt1(p) > 0.0) && (dt2(p) > 0.0)){//Ce point appartient au contour interieur
				flag = 0;
				rmp = min(dt1(p),dt2(p));
				rMp = max(dt1(p),dt2(p));
				rcmax = 0.0;
				for (int a = 0 ; a <= maxX ; a++ ){
					for (int b = 0 ; b <= maxY ; b++ ){
						Z2i::Point q(a,b);
						if ( (dt1(q) > 0.0) && (dt2(q) > 0.0) && (p != q)){
							rmq = min(dt1(q),dt2(q));
							rMq = max(dt1(q),dt2(q));
							distance = dist(p,q);
							if (distance + rMp <= rMq){
								flag = 1;
							}
 						}
					}
				}
				if (!flag) {  
					result.setValue(p,1);
					//if (result(p) > maxv) maxv = result(p);
				} 
			}
		}
	}
	Board2D board;
	board.clear();
	Display2DFactory::drawImage<Gray>(board, result , 0.0 ,1);
	board.saveSVG(outputfile);  
	cout << maxv << endl;
	return result;
}

float power(Z2i::Point centre,float r,Z2i::Point p){
	return distSquare(p,centre)-r*r;
}

/* Construction de la carte de couverture
   Entrée : ensemble d'hyperboule
   Sortie : carte de couverture
*/

Grille hyperballCovering(Grille &dtInt,Grille &dtExt){
	Grille result(dtInt.domain());
	int maxX = dtExt.domain().upperBound()[0];
	int maxY = dtExt.domain().upperBound()[1];
	int infX,infY,supX,supY;
	for (int i=0;i<=maxX;i++){
		for (int j=0;j<=maxY;j++){
			Z2i::Point p(i,j);
			infX = max(0,i-int(sqrt(dtExt(p)))-1);
			supX = min(result.domain().upperBound()[0],i+int(sqrt(dtExt(p)))+1);
			infY = max(0,j-int(sqrt(dtExt(p)))-1);
			supY = min(result.domain().upperBound()[1],j+int(sqrt(dtExt(p)))+1);
			for (int a = infX;a <= supX;a++){
				for (int b = infY;b <= supY;b++){
					Z2i::Point q(a,b);
					if (distSquare(p,q) < dtExt(p) ){
						result.setValue(q,result(q)+max(dtExt(p)-dtInt(p),-power(p,dtExt(p),q)));
					}
				}	
			}	
		}
	}
	int maxv = 0.0;
	for (int i=0;i<=maxX;i++){
		for (int j=0;j<=maxY;j++){
			Z2i::Point p(i,j);
			if (result(p) > maxv ) maxv = result(p);
		}
	}
	cout << maxv << endl;
	return result;
}

/*Algorithme de l'axe médian à partir de la notion de couverture
   Entrée : ensemble d'hyperboules et un paramètre alpha
   Sortie : l'axe médian
*/
Grille imageAMFromRecovering(Grille &dtInt,Grille &dtExt,Grille &recovering,char *outputfile){
	Grille result(recovering.domain());
	Grille AM(recovering.domain());
	Grille taken(recovering.domain());
	int maxX = dtExt.domain().upperBound()[0];
	int maxY = dtExt.domain().upperBound()[1];
	//Trier les boules 
	vector<Z2i::Point> boules;
	Z2i::Point ptMax;
	for (int a =0;a< maxX*maxY;a++){
		ptMax = Z2i::Point(0,0);
		for (int i=0;i<=maxX;i++){
			for (int j=0;j<=maxY;j++){
				Z2i::Point p(i,j);
				if (dtExt(ptMax)+dtInt(ptMax) < dtExt(p)+dtInt(p) && !(taken(p)) ) ptMax = p;
			}
		}
		if (!dtExt(ptMax)) break;
		boules.push_back(ptMax);
		taken.setValue(ptMax,1);
		//cout << ptMax << " " << dtExt(ptMax) << endl; 
	}
	// Recherche de l'axe médian
	int flag;
	int infX,infY,supX,supY;
	int maxv=0;
	for (int a = boules.size()-1;a>=0;a--){
	//for (int a=0;a<boules.size();a++){
		flag = 0;
		infX = max(0,boules[a][0]-int(sqrt(dtExt(boules[a])))-1);
		supX = min(result.domain().upperBound()[0],boules[a][0]+int(sqrt(dtExt(boules[a])))+1);
		infY = max(0,boules[a][1]-int(sqrt(dtExt(boules[a])))-1);
		supY = min(result.domain().upperBound()[1],boules[a][1]+int(sqrt(dtExt(boules[a])))+1);
		for (int i = infX ;i <=supX ;i++){
			for (int j = infY;j <= supY;j++){
				Z2i::Point p(i,j);
				// Si le point appartient à l'hyperboule
				if (distSquare(boules[a],p) < dtExt(boules[a]) ){
					// Si après enlèvement de l'hyperboule un point n'est plus couvert
					if (recovering(p) <= max(dtExt(p)-dtInt(p),-power(boules[a],dtExt(boules[a]),p)) ) {
						flag = 1;
						//cout << recovering(p) << " " << max(dtExt(p)-dtInt(p),-power(boules[a],dtExt(boules[a]),p)) << endl;				 
						break;
					}
				}
			}
		}
		if (flag){
			result.setValue(boules[a],dtInt(boules[a]));
			AM.setValue(boules[a],1);
			if (result(boules[a]) > maxv) maxv = result(boules[a]);
		}
		else{
			for (int i = infX ;i <=supX ;i++){
				for (int j = infY;j <= supY;j++){
					Z2i::Point p(i,j);
					if (distSquare(boules[a],p) < dtExt(boules[a]) ){
						//on retire l'hyperboule
						recovering.setValue(p,recovering(p) - max(dtExt(p)-dtInt(p),-power(boules[a],dtExt(boules[a]),p)) );
					}
				}
			}
		}
	}
	Board2D b;
	b.clear();
	Display2DFactory::drawImage<Gray>(b,AM,0,1);
	b.saveSVG(outputfile);
	cout << maxv << endl;
	return AM;
}

/* Algorithme de l'axe médian à partir de l'axe médian sur la forme intérieure
   Entrée : ensemble d'hyperboules et un paramètre alpha
   Sortie : l'axe médian 
*/

Grille imageAMfromMin(Grille &dtInt,Grille &dtExt,Grille &dtNeg,float alpha,char *outputfile){
	Grille result(dtInt.domain());
	int maxX = dtExt.domain().upperBound()[0];
	int maxY = dtExt.domain().upperBound()[1];
	Grille AMInt = imageAM(dtInt," ");
	Grille AMExt = imageAM(dtExt," ");
	int maxv = 0;
	for (int i = 0 ;i <=maxX ;i++){
		for (int j = 0 ;j <= maxY;j++){
			Z2i::Point p(i,j);
			if ((dtInt(p) == 0) && (AMExt(p) != 0)) {
				/*
				distance totale = d
				d = dtExt(p) + dtNeg(p)
				alpha = dtNeg(p) / d
				*/
				if ((dtNeg(p) > 0) && ((100*dtNeg(p)) / (dtExt(p)+dtNeg(p)) <= 100*alpha) ) 
					result.setValue(p,1);
			}
			else {
				if (AMInt(p) != 0) result.setValue(p,1);
			}
			//if (result(p) > maxv) maxv = result(p);
		}
	}
	Board2D b;
	b.clear();
	Display2DFactory::drawImage<Gray>(b,result,0,1);
	b.saveSVG(outputfile);

	return result;	
} 

/*
	Opère la restriction d'un ensemble de points sur un autre
	Entrée : un axe médian et une DT
	Sortie : l'axe médian réduit sur les point où la DT est strictement positive  
*/
Grille restriction(Grille &ma,Grille &dt,char * outputfile){
	Grille result(dt.domain());
	int maxX = dt.domain().upperBound()[0];
	int maxY = dt.domain().upperBound()[1];
	for (int i = 0 ;i <=maxX ;i++){
		for (int j = 0 ;j <= maxY;j++){
			Z2i::Point p(i,j);
			if ((ma(p) > 0) && (dt(p) > 0)) result.setValue(p,1);
		}
	}
	Board2D b;
	b.clear();
	Display2DFactory::drawImage<Gray>(b,result,0,1);
	b.saveSVG(outputfile);
	return result;
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
	
	for (int i =3;i<=27;i++){
		for (int j =3;j<=27;j++){
			test.setValue(Z2i::Point(i,j),128);
			
		}
	}
	
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
	Grille skeleton = homotopicThinning(imageN,"../../../skeleton.svg");
	fill(skeleton);
	Voronoi2D VM_ske = voronoiMap(skeleton,"",""); 
	Grille DT_ske = imageDT(skeleton,"../../../DT_skeleton.svg");
	Grille DTVM_ske = imageDT(VM_ske,"../../../DTVM_skeleton.svg");
	Grille AM_skel = imageAM(DTVM_ske,"../../../AM_skeleton.svg");
	

	//récupération des contours intérieurs et extérieurs
	
	vector< vector< Z2i::Point >  >  contours;

	getContours(contours,image,x,y,noiseLevel);
	Grille image1 ( Z2i::Domain(lower,upper));
	Grille image2 ( Z2i::Domain(lower,upper));
	Grille imageCB ( Z2i::Domain(lower,upper));
	
	constructImage(contours[0],image1);
	constructImage(contours[1],image2);

	Grille image1filled = image1;
	Grille image2filled = image2;
	fill(image1filled);
	fill(image2filled);
	
	

	Board2D boardCB;
	Display2DFactory::drawImage<Gray>(boardCB, imageCB, (unsigned int)0, (unsigned int)129);
	boardCB.saveSVG("../../../contourBoth.svg");
	
	Board2D board1,board2;
	Display2DFactory::drawImage<Gray>(board1, image1filled, (unsigned int)0, (unsigned int)129);
	Display2DFactory::drawImage<Gray>(board2, image2filled, (unsigned int)0, (unsigned int)129);
	board1.saveSVG("../../../contour1.svg");
	board2.saveSVG("../../../contour2.svg");


	/*Voronoi2D VM1 = voronoiMap(image1filled,"../../../VM_contour1.svg","../../../VMCells_contour1.svg"); 
	Voronoi2D VM2 = voronoiMap(image2filled,"../../../VM_contour2.svg","../../../VMCells_contour2.svg"); 

	//Calcul de la transformée en distance 
	
	Grille DT1 = imageDT(image1filled,"../../../DT_contour1.svg");
	Grille DT2 = imageDT(image2filled,"../../../DT_contour2.svg");

	Grille DTVM1 = imageDT(VM1,"../../../DT_VM_contour1.svg");
	Grille DTVM2 = imageDT(VM2,"../../../DT_VM_contour2.svg");


	Grille DTInt(DT1.domain());
	Grille DTExt(DT1.domain());
	//distingueIntExt(DT1,DT2,DTInt,DTExt);
	distingueIntExt(DTVM1,DTVM2,DTInt,DTExt);
	allFAlpha(DT1,DT2,0.05,"../../../AM_union_r=3.svg","",3);*/
	/*printImageColor(DTInt,"../../../DT_contourInt.svg");
	printImageColor(DTExt,"../../../DT_contourExt.svg");
	//Calcul de la moyenne et la différence des 2DT
	
	computeDTAverage(DT1,DT2,"../../../DT_contourAv.svg");
	computeDTDiff(DT1,DT2,"../../../DT_contourDiff.svg");

	// Calcul de l'axe médian
	Grille MA1 = imageAM(DTVM1,"../../../MA_contour1.svg");
	Grille MA2 = imageAM(DTVM2,"../../../MA_contour2.svg");

	// Reconstitution de la forme
	Grille RDT1 = imageRDT(MA1,"../../../RDT_MA1.svg");
	Grille RDT2 = imageRDT(MA2,"../../../RDT_MA2.svg");

	Grille rest = restriction(MA2,DTInt,"../../../MA2.svg");
	Grille r = imageRDT(rest,DTInt,"../../../RDT_MA_restreint_min.svg");
	Grille r2 = imageRDT(rest,DTExt,"../../../RDT_MA_restreint_max.svg");

	Grille alphaEqual = imageAlphaMinEqual(DTInt,DTExt,"../../../alpha_equal.svg");*/

	//preuve(image1,DT1,DT1,MA1);

	//Carte des angles
	
	//Grille angles = imageAngle(image1,image2,VM1,VM2,"../../../angles.svg"); 
	
	//Test Axe Médian du contour 

	//Grille DT_anneau = imageDT(imageN,"../../../DT_anneau.svg");
	//Grille MA_anneau = imageAM(DT_anneau,"../../../MA_anneau.svg");
	
	//Algos sur contours imprécis
		//Algo des hyperboules
	//Grille HMA = imageHAM (DT1,DT2,"../../../HAM_contour.svg");
	//Grille RDT_HMA = imageRDT(HMA,"../../../RDT_HMA.svg");
	//Grille RDT_DYM_HMA = RDTDynamiqueWithAM(HMA,DTVM1,DTVM2);
	//Grille RDT_alpha_HMA = imageAlphaMinWithAM(HMA,DTVM1,DTVM2,DTVM1,"../../../HAM_alphas.svg");
	
		//Algo de l'alpha axe médian
	//Grille alphaMA = imageAlphaAM(DT1,DT2,0.0,"../../../alphaAM_alpha=0.0.svg");
	//Grille alphaMA2 = imageAlphaAM(DT1,DT2,1.0,"../../../alphaAM_alpha=1.0.svg");
	//Grille alphaMA3 = imageAlphaAM(DT1,DT2,0.5,"../../../alphaAM_alpha=0.5.svg");
	//allFAlpha(DT1,DT2,0.1,"../../../allFAlphaAM.svg","../../../allFAlphaRDT.svg");

	
    	//printImageNB(DTIntNeg,"../../../Inverse.svg");
	/*Grille DTIntNeg = negatif(DTInt);
	Voronoi2D VMIntNeg = voronoiMap(DTIntNeg,"","");
    	Grille DTIntInv = imageDT(VMIntNeg,"../../../DT_inv_contourInt.svg");
    	Grille alphaMin = imageAlphaMin(DTInt,DTExt,DTIntInv,"../../../alphaMinEucl.svg");
	Grille alphaMinDyn = imageAlphaMinDynamique(DTInt,DTExt,"../../../alphaMinDynamique.svg");*/
	
	//Grille RDT_alphaMA = imageRDT(alphaMA,"../../../RDT_alphaAM_alpha=0.0.svg");
	//Grille RDT_alphaMA2 = imageRDT(alphaMA2,"../../../RDT_alphaAM_alpha=1.0.svg");
	//Grille RDT_alphaMA3 = imageRDT(alphaMA3,"../../../RDT_alphaAM_alpha=0.5.svg");
		//Algo de l'alpha-beta axe médian
	//Grille alphaBetaMA = imageAlphaBetaAM(DT1,DT2,0.0,0.0,"../../../alphabetaAM_alpha=0.0_beta=0.0.svg");
	//Grille alphaBetaMA2 = imageAlphaBetaAM(DT1,DT2,1.0,1.0,"../../../alphabetaAM_alpha=1.0_beta=1.0.svg");

	//Grille RDT_alphaBetaMA = imageRDT(alphaBetaMA,"../../../RDT_alphabetaAM_alpha=0.0_beta=0.0.svg");
	//Grille RDT_alphaBetaMA2 = imageRDT(alphaBetaMA2,"../../../RDT_alphabetaAM_alpha=1.0_beta=1.0.svg");
	
	//Grille RDTdym = RDTDynamique(DTInt,DTExt);
	//cout << nbComposantesConnexes(DTInt) << endl ;
	//Voronoi2D VMIntNeg = voronoiMap(DTIntNeg,"","");
    	//Grille DTIntInv = imageDT(VMIntNeg,"../../../DT_inv_contourInt.svg");
	/*Grille alphasGlo = alphaFromIntGlobal(DTExt,DTIntInv,"../../../alphasGlobal.svg");
	Grille alphasLoc = alphaFromIntLocal(DTExt,DTIntInv,"../../../alphasLocal.svg");
	vector<float> alphasLimLoc;
	alphas(alphasLimLoc,RDTdym,alphasLoc);
	for (int i =0;i< alphasLimLoc.size();i++) cout << alphasLimLoc[i] << " "  ;
	cout << endl;
	vector<float> alphasLimGlo;
	alphas(alphasLimGlo,RDTdym,alphasGlo);
	for (int i =0;i< alphasLimGlo.size();i++) cout << alphasLimGlo[i] << " "  ;
	cout << endl;
	Grille AMDym = MAfromRDTDym(RDTdym,"../../../AM_Dym");*/
	//Grille recouvrement = hyperballCovering(DTInt,DTExt);
	//Grille AMcovering = imageAMFromRecovering(DTInt,DTExt,recouvrement,"../../../AM_covering.svg");
	//Grille RDT_DYM_HMA = RDTDynamiqueWithAM(AMcovering,DTVM1,DTVM2);
	//Grille RDT_alpha_cover = imageAlphaMinWithAM(AMcovering,DTVM1,DTVM2,DTVM1,"../../../covering_alphas.svg");
	//Grille RDTfromRecovering = imageRDT(AMcovering,"../../../RDT_covering_Min.svg");
	//Grille AMFromAMMin = imageAMfromMin(DTInt,DTExt,DTIntInv,0.0,"../../../AM_AMMin_0.svg");
	//Grille AMFromAMMin2 = imageAMfromMin(DTInt,DTExt,DTIntInv,0.5,"../../../AM_AMMin_0,5.svg");
	//Grille AMFromAMMin3 = imageAMfromMin(DTInt,DTExt,DTIntInv,1.0,"../../../AM_AMMin_1,0.svg");
	//Grille RDT_DYM_add = RDTDynamiqueWithAM(AMFromAMMin,DTVM1,DTVM2,"../../../RDT_Dym_add_0,0.svg");
	//Grille RDT_DYM_add2 = RDTDynamiqueWithAM(AMFromAMMin2,DTVM1,DTVM2,"../../../RDT_Dym_add_0,5.svg");
	//Grille RDT_DYM_add3 = RDTDynamiqueWithAM(AMFromAMMin3,DTVM1,DTVM2,"../../../RDT_Dym_add_1,0.svg");
	//Grille RDT_alpha_add = imageAlphaMinWithAM(AMFromAMMin,DTVM1,DTVM2,DTVM1,"../../../adding_alphas_0,0.svg");
	//Grille RDT_alpha_add2 = imageAlphaMinWithAM(AMFromAMMin2,DTVM1,DTVM2,DTVM1,"../../../adding_alphas_0,5.svg");
	//Grille RDT_alpha_add3 = imageAlphaMinWithAM(AMFromAMMin3,DTVM1,DTVM2,DTVM1,"../../../adding_alphas_1,0.svg");
	//Grille RDT_FromMin = imageRDT(AMFromAMMin,"../../../RDT_FromMin_0.svg");
	//Grille RDT_FromMin2 = imageRDT(AMFromAMMin2,"../../../RDT_FromMin_0,5.svg");
	//Grille RDT_FromMin3 = imageRDT(AMFromAMMin3,"../../../RDT_FromMin_1,0.svg");
	//Création d'un contour avec boules réduites
	//contourBallReduced(x,y,noiseLevel,freeman,valMaxX+20,valMaxY+20); 
	
	return 0; 
}
