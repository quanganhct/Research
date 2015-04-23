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
using namespace Z2i;

typedef Z2i::Z2::Point MyPoint;

typedef std::vector<Z2i::Z2::Point> Container;
      // Iterator on the container
typedef Container::const_iterator ConstIterator;
      // StandardDSS4 computer
typedef StandardDSS4Computer<ConstIterator> DSSComputer;  

typedef DSSComputer::Primitive DSS;

int main()
{
    int width = 2;
    DSSComputer dss1, dss2;

    Container contour;
    contour.push_back(MyPoint(0,0));
    contour.push_back(MyPoint(1,0));
    contour.push_back(MyPoint(1,1));
    contour.push_back(MyPoint(2,1));
    contour.push_back(MyPoint(3,1));
    contour.push_back(MyPoint(3,2));
    contour.push_back(MyPoint(4,2));
    contour.push_back(MyPoint(5,2));
    contour.push_back(MyPoint(6,2));
    contour.push_back(MyPoint(6,3));
    contour.push_back(MyPoint(6,4));

    contour.push_back(MyPoint(6,4));
    contour.push_back(MyPoint(6,3));
    contour.push_back(MyPoint(6,2));
    contour.push_back(MyPoint(5,2));
    contour.push_back(MyPoint(4,2));
    contour.push_back(MyPoint(3,2));
    contour.push_back(MyPoint(3,1));
    contour.push_back(MyPoint(2,1)); 
    contour.push_back(MyPoint(1,1));    
    contour.push_back(MyPoint(1,0));
    contour.push_back(MyPoint(0,0));
    
    
    ConstIterator first = contour.begin();
    ConstIterator second = contour.end()-1;

    cout << *first << endl;
    cout << *second << endl;

    dss1.init( first );
    ConstIterator current = first;
    cout << "First round :" << endl;
    DSS d1 = dss1.primitive();
    while ((dss1.end() != contour.end()) && (dss1.extendFront())) {
        cout << *current << endl;
        d1 = dss1.primitive();
        cout << "A:" << d1.a() << " B:" << d1.b() << " U:" << d1.mu() << endl;
        current++;
    }

    dss2.init(second);
    current = second;
    cout << "Second round :" << endl;
    while ((dss2.begin() != contour.begin()) && (dss2.extendBack())){
        cout << *current << endl;
        d1 = dss2.primitive();
        cout << "A:" << d1.a() << " B:" << d1.b() << " U:" << d1.mu() << endl;
        current--;
    }

    //cout << dss1 << endl;
    DSSComputer::Primitive theDSS = dss1.primitive();  
   
    //cout << theDSS << endl;
    //aBoard.saveEPS("freemanChainFromImage.eps");
    return 0;
}