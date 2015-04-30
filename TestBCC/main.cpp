/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
/**
 * @file curvatureScaleSpace.cpp
 * @ingroup estimators
 *
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 *
 * @date 2014/04/02
 *
 * Output the image of curvature scale space
 * 
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <typeinfo>
#include "DGtal/base/Common.h"

#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/helpers/Surfaces.h"

//Grid curve
#include "DGtal/geometry/curves/FreemanChain.h"

//Estimators
#include "BinomialConvolver.h"


// Generation of resulting image
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/writers/PPMWriter.h"


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <vector>
#include <string>
#include <iomanip>
#define MEM_SIZE (128)
#define MAX_SOURCE_SIZE (0x100000)

using namespace DGtal;
using namespace std;
using namespace DGtal::Z2i;

typedef vector<Point>::const_iterator ConstIteratorOnPoints;
typedef vector<double>::const_iterator iter;


void displayVector(vector<double> data){
  int num = 0;
  for(iter i = data.begin(); i != data.end(); ++i){
    cout << "Vector n°" << num << " - " << *i << endl;
    num++;
  }

}

void computeCurvatureBCC(double h, const FreemanChain<Z2i::Integer> &fc, vector<double> &resCurvature,bool isClosed){
  vector<Z2i::Point> vectPoints;
  FreemanChain<Z2i::Integer>::getContourPoints( fc, vectPoints );
  typedef BinomialConvolver<vector<Z2i::Point>::const_iterator, double> myBC;
  typedef CurvatureFromBinomialConvolverFunctor< myBC, double >   myBCFCT;

  BinomialConvolverEstimator< myBC, myBCFCT> BCCurvatureEstimator;
  BCCurvatureEstimator.init( h, vectPoints.begin(), vectPoints.end(), isClosed );
  resCurvature.clear();
  resCurvature.resize(vectPoints.size());
  BCCurvatureEstimator.eval( vectPoints.begin(), vectPoints.end(), resCurvature.begin() );


/*
 *  General Variables
 */
//  unsigned int myN = BCCurvatureEstimator.getMyN();
//  double myH = BCCurvatureEstimator.getMyH();
//  Signal<double> myX = BCCurvatureEstimator.getMyX();
//  Signal<double> myY = BCCurvatureEstimator.getMyY();
//  Signal<double> myDX = BCCurvatureEstimator.getMyDX();
//  Signal<double> myDY = BCCurvatureEstimator.getMyDY();
//  Signal<double> myDDX = BCCurvatureEstimator.getMyDDX();
//  Signal<double> myDDY = BCCurvatureEstimator.getMyDDY();
//  ConstIteratorOnPoints myBegin;
//  ConstIteratorOnPoints myEnd;
//  cout << "myH is: " <<myH << endl;
//
//
//  std::vector<Z2i::Point>::const_iterator itb = vectPoints.begin();
//  std::vector<Z2i::Point>::const_iterator ite = vectPoints.end();
//
//
//  for ( std::vector<Z2i::Point>::const_iterator it = itb; it != ite; ++it ){
//    int i = 0;
//    typename std::map<ConstIteratorOnPoints,int>::const_iterator
//            map_it = myMapIt2Idx.find( it );
//    if ( map_it != myMapIt2Idx.end() )
//      i = map_it->second;
//    double denom = pow( myDX[ i ] * myDX[ i ] + myDY[ i ] * myDY[ i ], 1.5 );
//    double v = ( denom != double( 0.0 ) )? ( myDDX[ i ] * myDY[ i ] - myDDY[ i ] * myDX[ i ] ) / denom / myH : double( 0.0 );
//    std::cout << "Value: " << v << std::endl;
//    *resCurvature.begin()++ = v;
//  }
}




///////////////////////////////////////////////////////////////////////////////
namespace po = boost::program_options;

int main( int argc, char** argv )
{
 
  // parse command line ----------------------------------------------
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("FreemanChain,f", po::value<std::string>(), "FreemanChain file name")
    ("gridStepInit", po::value<double>()->default_value(1.0), "Grid step initial")
    ("gridStepIncrement", po::value<double>()->default_value(1.0), "Grid step increment ")
    ("output,o ", po::value<std::string>(), "set the output name ")
    ("gridStepFinal", po::value<double>()->default_value(1.0), "Grid step final")
    ("curvatureCutOff,c", po::value<double>()->default_value(10.0), "set the curvature limits to better display");
  

  typedef ImageContainerBySTLVector<Z2i::Domain, double > Image2D;


  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  }catch(const std::exception& ex){
    parseOK=false;
    trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
  }
  po::notify(vm);    
  if(!parseOK || vm.count("help")||argc<=1 || (!(vm.count("FreemanChain"))) )
    {
      trace.info()<< "Generate the Curvature Scale Space image using a binomial convolver based estimator." <<std::endl
                  << "The x axis is associated to the contour point and the y axis to the scale. The color represent the curvature values included between the cutoff values (set to 10 by default)."
                  << std::endl << "Basic usage: "<<std::endl
      << "\t main -f ../Research/images/contourS.fc --gridStepInit 0.001 --gridStepIncrement  0.0005 --gridStepFinal 0.05 -o result.ppm "<<std::endl
      << general_opt << "\n";
      return 0;
    }

  double h_initial = vm["gridStepInit"].as<double>();
  double h_increment = vm["gridStepIncrement"].as<double>();
  double h_final = vm["gridStepFinal"].as<double>();
  double curvatureCutOff = vm["curvatureCutOff"].as<double>();
  

  if(vm.count("FreemanChain")){
    string fileName = vm["FreemanChain"].as<std::string>();
    vector< DGtal::FreemanChain<Z2i::Integer>  > vectFcs =  PointListReader< Z2i::Point >:: getFreemanChainsFromFile<Z2i::Integer> (fileName);
    bool isClosed = vectFcs.at(0).isClosed(); 
    
    
    // Preparing resulting image:
    unsigned int height =  (int)((h_final-h_initial)/h_increment);
    // We add one point if the freemnchain is open.
    Z2i::Domain domain (Z2i::Point(0,0), Z2i::Point(vectFcs.at(0).size()+(isClosed? 0: 1), height-1));
    Image2D cssImage(domain);    
    HueShadeColorMap<double>  gradCurvature (-curvatureCutOff, curvatureCutOff);
    
   // trace.progressBar(0, height);
    double h= h_initial;
    for(double l= 0; l < height; l++ ){
      // Binomial estimator
      trace.progressBar(l, height);
      vector<double> curvaturesBCC;
      computeCurvatureBCC(h, vectFcs.at(0), curvaturesBCC, isClosed);
     // cout << "VectFcs.at(0)" << vectFcs.at(0) << endl;
      //displayVector(curvaturesBCC);
      // Output
      unsigned int j = 0;
        for ( std::vector<double>::const_iterator it = curvaturesBCC.begin(), it_end = curvaturesBCC.end(); it != it_end; ++it, ++j ){
          double c = *it;
          c = c<-curvatureCutOff? -curvatureCutOff: c;
          c = c>curvatureCutOff? curvatureCutOff: c;
          cssImage.setValue(Z2i::Point(j, l), c);
        }



      h=h+h_increment;
    }

    trace.progressBar(height, height);
    trace.info() <<std::endl;
    DGtal::GenericWriter<Image2D, 2, double, HueShadeColorMap<double> >::exportFile(vm["output"].as<std::string>(), cssImage, gradCurvature );
  }
  return 0;
}

