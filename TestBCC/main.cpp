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

#include "DGtal/base/Common.h"

#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/helpers/Surfaces.h"

//Grid curve
#include "DGtal/geometry/curves/FreemanChain.h"

//Estimators
#include "DGtal/geometry/curves/BinomialConvolver.h"


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
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#define MEM_SIZE (128)
#define MAX_SOURCE_SIZE (0x100000)
using namespace DGtal;


void
computeCurvatureBCC(double h, const FreemanChain<Z2i::Integer> &fc, std::vector<double> &resCurvature,bool isClosed){
 // std::cout<<"Hello World\n"; 
  std::vector<Z2i::Point> vectPoints;
  FreemanChain<Z2i::Integer>::getContourPoints( fc, vectPoints ); 
  
  typedef BinomialConvolver<std::vector<Z2i::Point>::const_iterator, double> MyBinomialConvolver;
  typedef CurvatureFromBinomialConvolverFunctor< MyBinomialConvolver, double >   CurvatureBCFct;

  BinomialConvolverEstimator< MyBinomialConvolver, CurvatureBCFct> BCCurvatureEstimator;
  BCCurvatureEstimator.init( h, vectPoints.begin(), vectPoints.end(), isClosed );
  resCurvature.clear();
  resCurvature.resize(vectPoints.size());
  BCCurvatureEstimator.eval( vectPoints.begin(), vectPoints.end(), resCurvature.begin() ); 
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
      << "\t main -f ../Research/cookbook/samples/contourS.fc --gridStepInit 0.001 --gridStepIncrement  0.0005 --gridStepFinal 0.05 -o result.ppm "<<std::endl
      << general_opt << "\n";
      return 0;
    }

  double h_initial = vm["gridStepInit"].as<double>();
  double h_increment = vm["gridStepIncrement"].as<double>();
  double h_final = vm["gridStepFinal"].as<double>();
  double curvatureCutOff = vm["curvatureCutOff"].as<double>();
  
/*

  OPEN CL PART


*/
  cl_device_id device_id = NULL;
  cl_context context = NULL;
  cl_command_queue command_queue = NULL;
  cl_mem memobj = NULL;
  cl_program program = NULL;
  cl_kernel kernel = NULL;
  cl_platform_id platform_id = NULL;
  cl_uint ret_num_devices;
  cl_uint ret_num_platforms;
  cl_int ret;
  char string[MEM_SIZE];

  FILE *fp;
  char fileName[] = "./hello.cl";
  char *source_str;
  size_t source_size;

/* Load the source code containing the kernel*/
  fp = fopen(fileName, "r");
  if (!fp) {
    fprintf(stderr, "Failed to load kernel.\n");
    exit(1);
  }
  source_str = (char*)malloc(MAX_SOURCE_SIZE);
  source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
  fclose(fp);

/* Get Platform and Device Info */
  ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
  ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, 1, &device_id, &ret_num_devices);

/* Create OpenCL context */
  context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);

/* Create Command Queue */
  command_queue = clCreateCommandQueue(context, device_id, 0, &ret);

/* Create Memory Buffer */
  memobj = clCreateBuffer(context, CL_MEM_READ_WRITE,MEM_SIZE * sizeof(char), NULL, &ret);

/* Create Kernel Program from the source */
  program = clCreateProgramWithSource(context, 1, (const char **)&source_str,
                                      (const size_t *)&source_size, &ret);

/* Build Kernel Program */
  ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);

/* Create OpenCL Kernel */
  kernel = clCreateKernel(program, "hello", &ret);

/* Set OpenCL Kernel Parameters */
  ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&memobj);

/* Execute OpenCL Kernel */
  ret = clEnqueueTask(command_queue, kernel, 0, NULL,NULL);

/* Copy results from the memory buffer */
  ret = clEnqueueReadBuffer(command_queue, memobj, CL_TRUE, 0,
                            MEM_SIZE * sizeof(char),string, 0, NULL, NULL);

/* Display Result */
  puts(string);

/* Finalization */
  ret = clFlush(command_queue);
  ret = clFinish(command_queue);
  ret = clReleaseKernel(kernel);
  ret = clReleaseProgram(program);
  ret = clReleaseMemObject(memobj);
  ret = clReleaseCommandQueue(command_queue);
  ret = clReleaseContext(context);

  free(source_str);




/*

  END OPEN CL PART

*/

  if(vm.count("FreemanChain")){
    std::string fileName = vm["FreemanChain"].as<std::string>();    
    std::vector< DGtal::FreemanChain<Z2i::Integer>  > vectFcs =  PointListReader< Z2i::Point >:: getFreemanChainsFromFile<Z2i::Integer> (fileName);     
    bool isClosed = vectFcs.at(0).isClosed(); 
    
    
    // Preparing resulting image:
    unsigned int height =  (int)((h_final-h_initial)/h_increment);
    // We add one point if the freemnchain is open.
    Z2i::Domain domain (Z2i::Point(0,0), Z2i::Point(vectFcs.at(0).size()+(isClosed? 0: 1), height-1));
    Image2D cssImage(domain);    
    HueShadeColorMap<double>  gradCurvature (-curvatureCutOff, curvatureCutOff);
    
    trace.progressBar(0, height);
    double h= h_initial;
    //bool flag = false;
    for(double l= 0; l < height; l++ ){
      // Binomial estimator
      trace.progressBar(l, height);
      std::vector<double> curvaturesBCC;
      computeCurvatureBCC(h, vectFcs.at(0), curvaturesBCC, isClosed);
      // Output
      unsigned int j = 0;
      //if(flag == false){
        for ( std::vector<double>::const_iterator it = curvaturesBCC.begin(), it_end = curvaturesBCC.end(); it != it_end; ++it, ++j ){
          double c = *it;
          c = c<-curvatureCutOff? -curvatureCutOff: c;
          c = c>curvatureCutOff? curvatureCutOff: c;
          cssImage.setValue(Z2i::Point(j, l), c); 
          //flag = true;
        }   
     // }
         

      h=h+h_increment;
    }
    trace.progressBar(height, height);
    trace.info() <<std::endl;
    DGtal::GenericWriter<Image2D, 2, double, HueShadeColorMap<double> >::exportFile(vm["output"].as<std::string>(), cssImage, gradCurvature );       
  }
  return 0;
}