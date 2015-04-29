#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/BasicTypes.h"
//#include "DGtal/io/readers/ITKReader.h"
#include "DGtal/io/readers/PGMReader.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"

#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/geometry/curves/FreemanChain.h"
#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"
#include "DGTal/io/writers/GenericWriter.h"

#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/Color.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <OpenCL/opencl.h>
#include <stdlib.h>

#define MEM_SIZE (128)
#define MAX_SOURCE_SIZE (0x100000)

using namespace DGtal;
using namespace std;
using namespace DGtal::Z2i;

typedef Z2i::Z2::Point MyPoint;

typedef std::vector<Z2i::Z2::Point> Container;
      // Iterator on the container
typedef Container::const_iterator ConstIterator;
      // StandardDSS4 computer
typedef  StandardDSS4Computer<ConstIterator> DSSComputer;  


pair<Container, int> adjustVector(Container contour){
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
    int distance = 0;
    while (iterator_back != tempt_b.begin()){
        iterator_back--;
        newContour.push_back(*iterator_back);
        distance++;
    }

    newContour.insert(newContour.end(), contour.begin(), contour.end());
    pair<Container, int> p;
    p.first = newContour;
    p.second = distance - 1;
    //cout << "Start Point :" << *startPoint << endl;

    return p;
}

size_t shrRoundUp(size_t localWorkSize, size_t numItems) {
    size_t result = localWorkSize;
    while (result < numItems)
        result += localWorkSize;
    
    return result;
}

void displayABW(DSSComputer dss){
    double a = dss.primitive().a();
    double b = dss.primitive().b();
    double w = dss.primitive().omega();
    cout << " A:" << a << " B:" << b << " W:" << w << endl;
}

double computeWidth(DSSComputer dss){
    double a = abs(dss.primitive().a());
    double b = abs(dss.primitive().b());
    double w = dss.primitive().omega();
    return double(w-1)/max(a, b);
}

double calculCurvature(Container contour, ConstIterator current, int width){
    DSSComputer front, back;
    ConstIterator it_f, it_b;

    front.init(current);
    back.init(current);

    it_f = current;
    it_b = current;

    while((front.end() != contour.end()) && front.extendFront()){
        if(computeWidth(front) > width){
            front.retractFront();
            break;
        }
        it_f++;
    }

    if (front.end() == contour.end()){
        it_f--;
    }

    while((back.begin() != contour.begin()) && back.extendBack()){
        if(computeWidth(back) > width){
            back.retractBack();
            break;
        }
        it_b--;
    }
    //cout << *current << endl;
    //cout << " V front:" << abs(it_f-current) << " points " << computeWidth(front); displayABW(front);
    //cout << " V back:" << abs(it_b-current) << " points " << computeWidth(back); displayABW(back);


    if (back.begin() == contour.begin()){
        //it_b++;
    }

    //cout << "BackBegin:" << *back.begin() << endl;
    //cout << "Back A:" << back.primitive().a() << " B:" << back.primitive().b() << " U:" << back.primitive().mu() << " W:" << back.primitive().omega() << endl;

    MyPoint k = *current;
    MyPoint l = *it_b;
    MyPoint r = *it_f;

    double kl = sqrt(pow(k.myArray[0]-l.myArray[0], 2) + pow(k.myArray[1]-l.myArray[1], 2));
    double lr = sqrt(pow(l.myArray[0]-r.myArray[0], 2) + pow(l.myArray[1]-r.myArray[1], 2));
    double kr = sqrt(pow(k.myArray[0]-r.myArray[0], 2) + pow(k.myArray[1]-r.myArray[1], 2));

    double rayon = kl * lr * kr / sqrt((kl+lr+kr)*(kl+lr-kr)*(kl+kr-lr)*(kr+lr-kl));
    int det = (r.myArray[0]-k.myArray[0])*(l.myArray[1]-k.myArray[1]) - (l.myArray[0]-k.myArray[0])*(r.myArray[1]-k.myArray[1]);
    int sign = (det == 0) ? 0 : (det > 0) ? 1 : -1;
    //cout << "Rayon:" << rayon << " K:" << *current << " R:" << *it_f << " L:" << *it_b << endl;
    //cout << "KL:" << kl << " LR:" << lr << " KR:" << kr << endl;
    return sign/rayon;
}

void calculCurvatureForKernel(Container contour, ConstIterator current, int width, std::vector<int> *k_points
                                , std::vector<int> *l_points, std::vector<int> *r_points){
    DSSComputer front, back;
    ConstIterator it_f, it_b;

    front.init(current);
    back.init(current);

    it_f = current;
    it_b = current;

    while((front.end() != contour.end()) && front.extendFront()){
        if(computeWidth(front) > width){
            front.retractFront();
            break;
        }
        it_f++;
    }

    if (front.end() == contour.end()){
        it_f--;
    }

    while((back.begin() != contour.begin()) && back.extendBack()){
        if(computeWidth(back) > width){
            back.retractBack();
            break;
        }
        it_b--;
    }
    //cout << *current << endl;
    //cout << " V front:" << abs(it_f-current) << " points " << computeWidth(front); displayABW(front);
    //cout << " V back:" << abs(it_b-current) << " points " << computeWidth(back); displayABW(back);


    if (back.begin() == contour.begin()){
        //it_b++;
    }

    MyPoint k = *current;
    MyPoint l = *it_b;
    MyPoint r = *it_f;
    //cout << "K " << k.myArray[0] << " " << k.myArray[1] <<  endl;

    (*k_points).push_back(k.myArray[0]); (*k_points).push_back(k.myArray[1]);
    (*l_points).push_back(l.myArray[0]); (*l_points).push_back(l.myArray[1]);
    (*r_points).push_back(r.myArray[0]); (*r_points).push_back(r.myArray[1]);
}

void calculateVectorCurvature(int* l, int* k, int* r, float* curvature, int size){
    cl_device_id device_id = NULL;
    cl_context context = NULL;
    cl_command_queue command_queue = NULL;
    cl_program program = NULL;
    cl_kernel kernel = NULL;
    cl_platform_id platform_id = NULL;
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_int error;

    char string[MEM_SIZE];
 
    FILE *fp;
    char fileName[] = "./dsl.cl";
    char *source_str;
    size_t source_size;

    fp = fopen(fileName, "r");
    if (!fp) {
        fprintf(stderr, "Failed to load kernel.\n");
        exit(1);
    }
    source_str = (char*)malloc(MAX_SOURCE_SIZE);
    source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
    fclose(fp);

    error = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
    error = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, &ret_num_devices);
 
    /* Create OpenCL context */
    context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &error);
 
    /* Create Command Queue */
    command_queue = clCreateCommandQueue(context, device_id, 0, &error);
 
    const int mem_size = sizeof(int)*size*2;
    cl_mem k_points = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, mem_size, k, &error);
    cl_mem r_points = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, mem_size, r, &error);
    cl_mem l_points = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, mem_size, l, &error);

    cl_mem vectCourvature = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_float)*size, NULL, &error);

    program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &error);
    assert(error == CL_SUCCESS);

    error = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
    //error = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);

    size_t len;
    char buffer[9048];
    clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
    printf("%s\n", buffer);

    assert(error == CL_SUCCESS);

    kernel = clCreateKernel(program, "curvature", &error);
    assert(error == CL_SUCCESS);

    error = clSetKernelArg(kernel, 0, sizeof(cl_mem), &l_points);
    error |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &k_points);
    error |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &r_points);
    error |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &vectCourvature);
    error |= clSetKernelArg(kernel, 4, sizeof(int), &size);
    assert(error == CL_SUCCESS);

    const size_t local_ws = 64;
    const size_t global_ws = shrRoundUp(local_ws, size);

    cout << "groups :" << global_ws << endl;

    error = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global_ws, &local_ws, 0, NULL, NULL);

    std::vector<int> l_error;
    l_error.push_back(CL_INVALID_PROGRAM_EXECUTABLE);
    l_error.push_back(CL_INVALID_COMMAND_QUEUE);
    l_error.push_back(CL_INVALID_KERNEL);
    l_error.push_back(CL_INVALID_CONTEXT);
    l_error.push_back(CL_INVALID_KERNEL_ARGS);
    l_error.push_back(CL_INVALID_WORK_DIMENSION);
    l_error.push_back(CL_INVALID_WORK_GROUP_SIZE);
    l_error.push_back(CL_INVALID_WORK_ITEM_SIZE);
    l_error.push_back(CL_INVALID_GLOBAL_OFFSET);
    l_error.push_back(CL_OUT_OF_RESOURCES);
    l_error.push_back(CL_MEM_OBJECT_ALLOCATION_FAILURE);
    l_error.push_back(CL_INVALID_EVENT_WAIT_LIST);
    l_error.push_back(CL_OUT_OF_HOST_MEMORY);

    //cout << "error :" << error << endl;
    //cout << "CL_DEVICE_MAX_WORK_GROUP_SIZE " << CL_DEVICE_MAX_WORK_GROUP_SIZE << endl;
    /*
    for(int i=0; i<l_error.size(); i++){
        cout << l_error[i] << endl;
    }
    */
    //assert(error == CL_INVALID_WORK_GROUP_SIZE);
    assert(error == CL_SUCCESS);

    clEnqueueReadBuffer(command_queue, vectCourvature, CL_TRUE, 0, size, curvature, 0, NULL, NULL);

    clReleaseMemObject(k_points);
    clReleaseMemObject(l_points);
    clReleaseMemObject(r_points);
    clReleaseMemObject(vectCourvature);

    clReleaseKernel(kernel);
    clReleaseCommandQueue(command_queue);
    clReleaseContext(context);
}

int main()
{
    typedef DGtal::ImageContainerBySTLMap<DGtal::Z2i::Domain, unsigned char> Image2D;

    int width = 1;
    
    typedef DGtal::ImageContainerBySTLVector< DGtal::Z2i::Domain, unsigned char> Image;
    std::string filename = "../Research/images/rabbit.pgm";
    Image image = DGtal::PGMReader<Image>::importPGM(filename);;
     
    Z2i::KSpace ks;
    ks.init( image.domain().lowerBound(), image.domain().upperBound(), true );
    
    Z2i::DigitalSet set2d (image.domain());
    SetFromImage<Z2i::DigitalSet>::append<Image>(set2d, image, 1, 255);
    
    SurfelAdjacency<2> sAdj( true );
    
    std::vector<std::vector< Z2i::Point>>  m_vector;
    Surfaces<Z2i::KSpace>::extractAllPointContours4C( m_vector, ks, set2d, sAdj);
    
    
    // Construction of the computer
    DSSComputer theDSSComputer;  

    Container contour;  

    //add the whole round contour
    for (int j=0; j<m_vector.at(0).size(); j++){
        //cout << m_vector.at(0).at(j) << " x y : " << m_vector.at(0).at(j).myArray[0] << ":" << m_vector.at(0).at(j).myArray[1] << endl;
        contour.push_back(m_vector.at(0).at(j));
    }

    int size = contour.size();

    pair<Container, int> p = adjustVector(contour);
    contour = p.first;
 
    ConstIterator startPoint = contour.begin() + p.second;
    ConstIterator current = contour.begin();

    cout << "Start Point : " << *startPoint << endl;

    while(current != contour.end()){
        cout << *current << endl;
        current++;
    }

    vector<double> vectorCouvature;

    ConstIterator startPoint2 = startPoint;
    ConstIterator i = startPoint;
    
    
    for (int i=0; i<size; i++){
        double value = calculCurvature(contour, startPoint, width);
        vectorCouvature.push_back(value);
        cout << *startPoint << " curvature " << value << endl;
        startPoint++;
    }
    
    /*
    std::vector<int> k_points;
    std::vector<int> l_points;
    std::vector<int> r_points;
    std::vector<float> curvature(size);

    for (int i = 0; i < size; ++i)
    {
        calculCurvatureForKernel(contour, startPoint2, width, &k_points, &l_points, &r_points);
        startPoint2++;
    }

    int* k_array = &k_points[0];
    int* l_array = &l_points[0];
    int* r_array = &r_points[0];
    //float* c_array = &curvature[0];

    float* c_array = new float[size];
    
    for (int i=0; i<k_points.size(); i++){
        cout << "K " << k_points[i] << endl;
    }
    
    cout << "SIZE " << size << endl;

    calculateVectorCurvature(l_array, k_array, r_array, c_array, size);
    
    for (int i=0; i<size; i++){
        cout << c_array[i] << endl;
    }
    */
    
    cout << "size:" << vectorCouvature.size() << endl;
    Board2D board;
    Domain domain(MyPoint(0, -100), MyPoint(vectorCouvature.size()+1, 100));
    board << domain;
    for(int i=0; i<vectorCouvature.size(); i++){
        int c = int(vectorCouvature.at(i)*100.0);
        if (c>100){
            c=100;
        }
        if (c<-100){
            c=-100;
        }
      board << MyPoint(i, c);
    }

    board.saveEPS("result_rabbit.eps");
    



    return 0;
}