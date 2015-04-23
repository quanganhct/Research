#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/BasicTypes.h"
//#include "DGtal/io/readers/ITKReader.h"
#include "DGtal/io/readers/PGMReader.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"

#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/geometry/curves/FreemanChain.h"
#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"

#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/Color.h"
#include <iostream>
#include <math.h>
#include <cmath>

using namespace DGtal;
using namespace std;

typedef Z2i::Z2::Point MyPoint;

typedef std::vector<Z2i::Z2::Point> Container;
      // Iterator on the container
typedef Container::const_iterator ConstIterator;
      // StandardDSS4 computer
typedef StandardDSS4Computer<ConstIterator> DSSComputer;  


pair<Container, ConstIterator> adjustVector(Container contour){
    DSSComputer front, back;
    Container tempt_f, tempt_b;

    front.init(contour.begin());
    back.init(contour.end()-1);

    cout << "BackBegin:" << *back.begin() << endl;
    cout << "Back A:" << back.primitive().a() << " B:" << back.primitive().b() << " U:" << back.primitive().mu() << " W:" << back.primitive().omega() << endl;

    ConstIterator current = contour.begin();

    while ((front.end() != contour.end()) && (front.extendFront())) {
        //cout << "Iter " << *current << endl;
        cout << "PointFront:" << front.front() << endl;
        //cout << "PointFront Back" << front.back() << endl;
        //current++;
        //cout << "Current " << *current << endl;
        tempt_f.push_back(front.front());
    }

    //back.init(current);
    //cout << "Back contour " << *(contour.end()) << " | " << *(contour.end()-1) << endl;


    while((back.begin() != contour.begin()) && (back.extendBack())) {

        cout << "PointBack:" << back.back() << endl;
        tempt_b.push_back(back.back());
        cout << "Back A:" << back.primitive().a() << " B:" << back.primitive().b() << " U:" << back.primitive().mu() << " W:" << back.primitive().omega() << endl;
    }

    ConstIterator iterator_front = tempt_f.begin();
    while (iterator_front != tempt_f.end()){
        contour.push_back(*iterator_front);
        iterator_front++;
    }

    ConstIterator iterator_back = tempt_b.end();

    Container newContour;

    while (iterator_back != tempt_b.begin()){
        iterator_back--;
        newContour.push_back(*iterator_back);
    }
    ConstIterator startPoint = newContour.end()-1;

    newContour.insert(newContour.end(), contour.begin(), contour.end());
    pair<Container, ConstIterator> p;
    p.first = newContour;
    p.second = startPoint;

    return p;
}

double calculCurvature(Container contour, ConstIterator current){
    DSSComputer front, back;
    ConstIterator it_f, it_b;

    front.init(current);
    back.init(current);

    it_f = current;
    it_b = current;

    while((front.end() != contour.end()) && front.extendFront()){
        it_f++;
    }
    if (front.end() == contour.end()){
        it_f--;
    }

    while((back.begin() != contour.begin()) && back.extendBack()){
        it_b++;
    }

    if (back.begin() == contour.begin()){
        it_b--;
    }

    cout << "BackBegin:" << *back.begin() << endl;
    cout << "Back A:" << back.primitive().a() << " B:" << back.primitive().b() << " U:" << back.primitive().mu() << " W:" << back.primitive().omega() << endl;

    MyPoint k = *current;
    MyPoint l = *it_b;
    MyPoint r = *it_f;

    double kl = sqrt(pow(k.myArray[0]-l.myArray[0], 2) + pow(k.myArray[1]-l.myArray[1], 2));
    double lr = sqrt(pow(l.myArray[0]-r.myArray[0], 2) + pow(l.myArray[1]-r.myArray[1], 2));
    double kr = sqrt(pow(k.myArray[0]-r.myArray[0], 2) + pow(k.myArray[1]-r.myArray[1], 2));

    double rayon = kl * lr * kr / sqrt((kl+lr+kr)*(kl+lr-kr)*(kl+kr-lr)*(kr+lr-kl));
    int det = (r.myArray[0]-k.myArray[0])*(l.myArray[1]-k.myArray[1]) - (l.myArray[0]-k.myArray[0])*(r.myArray[1]-k.myArray[1]);
    int sign = (det == 0) ? 0 : (det > 0) ? 1 : -1;
    cout << "Rayon:" << rayon << " K:" << *current << " R:" << *it_f << " L:" << *it_b << endl;
    cout << "KL:" << kl << " LR:" << lr << " KR:" << kr << endl;
    return sign/rayon;
}

int main()
{
    int width = 2;
    
    typedef DGtal::ImageContainerBySTLVector< DGtal::Z2i::Domain, unsigned char> Image;
    std::string filename = "../Research/images/QAcircle.pgm";
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
    
    std::vector<std::vector< Z2i::Point>>  m_vector;
    Surfaces<Z2i::KSpace>::extractAllPointContours4C( m_vector, ks, set2d, sAdj);
    
    
    //GradientColorMap<int> cmap_grad( 0, (const int)vectContoursBdryPointels.size() );
    //cmap_grad.addColor( Color( 50, 50, 255 ) );
    //cmap_grad.addColor( Color( 255, 0, 0 ) );
    //cmap_grad.addColor( Color( 255, 255, 10 ) );
    //cmap_grad.addColor( Color( 25, 255, 255 ) );
    //cmap_grad.addColor( Color( 255, 25, 255 ) );
    //cmap_grad.addColor( Color( 25, 25, 25 ) );
    
    // Construction of the computer
    DSSComputer theDSSComputer;  

    Container contour;  

    //add the whole round contour
    for (int j=0; j<m_vector.at(0).size(); j++){
        //cout << m_vector.at(0).at(j) << " x y : " << m_vector.at(0).at(j).myArray[0] << ":" << m_vector.at(0).at(j).myArray[1] << endl;
        contour.push_back(m_vector.at(0).at(j));
    }

    pair<Container, ConstIterator> p = adjustVector(contour);
    contour = p.first;
 
    ConstIterator startPoint = p.second;
    ConstIterator current = contour.begin();

    cout << "Start Point : " << *startPoint << endl;

    while(current != contour.end()){
        cout << *current << endl;
        current++;
    }

    vector<double> vectorCouvature;

    for (int i=0; i<10; i++){
        //cout << *startPoint << " " << calculCurvature(contour, startPoint) << endl;
        startPoint++;
    }
    /*
    theDSSComputer.init( current );
    while ( ( theDSSComputer.end() != contour.end() ) &&
          ( theDSSComputer.extendFront() ) ) {
        cout << "Iter " << *current << endl;
        cout << "PointFront " << theDSSComputer.front() << endl;
    }
    
      // Trace to the standard output
    cout << theDSSComputer << endl;

    DSSComputer::Primitive theDSS = theDSSComputer.primitive(); 
    cout << theDSS << endl; 
    */
    
    //aBoard.saveEPS("freemanChainFromImage.eps");
    return 0;
}