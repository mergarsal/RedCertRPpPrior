#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>


#include "RotPrior.h"
#include "RotPriorUtils.h"
#include "../utils/generatePointCloud.h"


using namespace std;
using namespace Eigen;
using namespace RotPrior; 


int main(int argc, char** argv)
{
    
  std::cout << "Essential Matrix Estimation with rotation prior!\n";



  double N_iter = 1;
  //set experiment parameters
  double noise = 0.0;
  const size_t n_points = 100;
  double FoV = 100;  // in degrees
  double max_parallax = 2.;  // in meters
  double min_depth = 1.;     // in meters
  double max_depth = 8.;       // in meters
  double outlier_fraction = 0;
  double focal_length = 800; 

   std::srand(std::time(nullptr));



       Vector3 translation;
       Matrix3 rotation;
       
       Eigen::MatrixXd points_3D(3, n_points);
       std::vector<int> indices_outliers(1, n_points);
 
       // define struct with params
       UtilsRelPose::StrInputGenProblem str_in = UtilsRelPose::StrInputGenProblem(); 
       str_in.FoV = FoV;  
       str_in.min_depth = min_depth; 
       str_in.max_depth = max_depth; 
       str_in.focal_length = focal_length; 
       str_in.n_points = n_points; 
       str_in.noise = noise; 
       str_in.parallax = max_parallax; 
                                               
       // generate problem
       UtilsRelPose::StrOutputGenProblem str_out = UtilsRelPose::createSyntheticExperiment(str_in, 
                                                   UtilsRelPose::generateRandomTranslationDefault, 
                                                   UtilsRelPose::generatePerturbRotationY); 
       // extract data
       translation = str_out.translation; 
       rotation = str_out.rotation; 
       
       Eigen::Matrix<double, 3, n_points> p0,p1;
       Eigen::Matrix<double, 1, n_points> weights; 
                              
                                  
       for (int i=0; i < n_points; i++)
       {
                UtilsRelPose::CorrespondingFeatures correspondence = str_out.points_correspondences[i];
                
                p0.col(i) = correspondence.bearing_vector_0; 
                p1.col(i) = correspondence.bearing_vector_1; 
                weights[i] = 1.0;                                
                points_3D.col(i) = str_out.points_3D.col(i);        
        }
                                               
                                               


      /* DLT + PATTERN */
      
      Matrix2 rot_init_pattern; 
      Vector3 trans_init_pattern; 
      Matrix3 E_init_prior; 
      bool valid_init = computeInitPattern(p0,p1,weights, E_init_prior, rot_init_pattern, trans_init_pattern);  
      Matrix3 Rtinit; 
      Rtinit.setZero();
      Rtinit.block<2,2>(0,0)=rot_init_pattern;
      Rtinit.block<3,1>(0,2)=trans_init_pattern;
      Rtinit.block<2,2>(0,0) << rotation(0,0), rotation(0,2), rotation(2,0), rotation(2,2);
      Rtinit.block<3,1>(0,2) << 1, 0, 0;
      // std::cout << "Init pose:\n" << "Rot:\n" << rot_init_pattern << std::endl; 
      // std::cout << "Trans:\n" << trans_init_pattern << std::endl;
      
            
      EssentialEstimationOptions options;
      options.verbose = 1;
      options.estimation_verbose=1;
      options.stepsize_tol = 1e-07;
                   
      
      /* Run essential with prior estimation and DLT + pattern*/       
      RotPriorClass ess_prior_est_pattern(p0,p1,weights, options); 
      // run the actual estimation
      EssentialEstimationResult my_result_prior_pattern = ess_prior_est_pattern.getResults(Rtinit); 
         



   
      std::cout << "Ground truth rotation:\n" << rotation << std::endl;
      std::cout << "Ground truth translation Tgt:\n" << translation << std::endl;
         
      
      std::cout << "+++++++++++++++\n\nResults from E with prior estimation PATTERN\n"; 
      // ess_prior_est_pattern.printResult(my_result_prior_pattern);
      std::cout << "Rotation:\n" << my_result_prior_pattern.R_opt << std::endl; 
      std::cout << "Traslation:\n" << my_result_prior_pattern.t_opt << std::endl;    
            
      
     // std::cout << "Errores for mani:\n"; 
     // std::cout << "Error rotation: " << distR(rotation, my_result_prior_pattern.R_opt) << std::endl;
     // std::cout << "Error translation: " << distT(translation, my_result_prior_pattern.t_opt) << std::endl;
      
 

  return 0;

}
