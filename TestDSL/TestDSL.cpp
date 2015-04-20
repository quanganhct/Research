#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
//#include "DGtal/io/readers/ITKReader.h"
#include "DGtal/io/readers/PGMReader.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"

#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/geometry/curves/FreemanChain.h"

#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/Color.h"
#include <iostream>

using namespace DGtal;

int main()
{
    
    typedef DGtal::ImageContainerBySTLVector< DGtal::Z2i::Domain, unsigned char> Image;
    std::string filename = "../Research/images/circle.pgm";
    Image image = DGtal::PGMReader<Image>::importPGM(filename);;
    //Image image = DGtal::ITKReader<Image>::importITK(filename);
    
    Z2i::KSpace ks;
    ks.init( image.domain().lowerBound(), image.domain().upperBound(), true );
    
    Z2i::DigitalSet set2d (image.domain());
    SetFromImage<Z2i::DigitalSet>::append<Image>(set2d, image, 1, 255);
    
    
    //Board2D aBoard;
    //aBoard << set2d;
    //aBoard << image.domain();
    
    SurfelAdjacency<2> sAdj( true );
    
    std::vector< std::vector< Z2i::Point >  >  m_vector;
    Surfaces<Z2i::KSpace>::extractAllPointContours4C( m_vector, ks, set2d, sAdj);
    
    
    
    //GradientColorMap<int> cmap_grad( 0, (const int)vectContoursBdryPointels.size() );
    //cmap_grad.addColor( Color( 50, 50, 255 ) );
    //cmap_grad.addColor( Color( 255, 0, 0 ) );
    //cmap_grad.addColor( Color( 255, 255, 10 ) );
    //cmap_grad.addColor( Color( 25, 255, 255 ) );
    //cmap_grad.addColor( Color( 255, 25, 255 ) );
    //cmap_grad.addColor( Color( 25, 25, 25 ) );
    
    for(unsigned int i=0; i<m_vector.size(); i++){
        //  Constructing and displaying FreemanChains from contours.
        //FreemanChain<Z2i::Integer> fc (vectContoursBdryPointels.at(i));
        //aBoard << SetMode( fc.className(), "InterGrid" );
        //aBoard<< CustomStyle( fc.className(),new CustomColors(  cmap_grad(i),  Color::None ) );
        //aBoard << fc;
        std::vector<Z2i::Point> vect = m_vector.at(i);
        std::cout << "Vect " << i << std::endl;
        for (unsigned int j=0; j<vect.size(); j++){
            std::cout << vect.at(j) << " x y : " << vect.at(j).myArray[0] << ":" << vect.at(j).myArray[1] << std::endl;
        }
    }
    
    //aBoard.saveEPS("freemanChainFromImage.eps");
    return 0;
}