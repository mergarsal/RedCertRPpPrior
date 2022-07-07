#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>


#include "RotPrior.h"
#include "RotPriorUtils.h"
#include "RotPriorConstraints.h"


// include problem generation
#include "../utils/generatePointCloud.h"

#include "IterCertAlg/SymmCert.h"
#include "IterCertAlg/SymmCert.cpp"
#include "IterCertAlg/RankYCert.h"
// #include "IterCertAlg/RankYCert.cpp"
#include "IterCertAlg/StaircaseCert.h"
#include "IterCertAlg/StaircaseCert.cpp"

#include <Eigen/Eigenvalues> 

using namespace std;
using namespace Eigen;
using namespace RotPrior; 

typedef Eigen::Matrix<double, 12, 1> Vector12;
typedef Eigen::Matrix<double, 13, 1> Vector13;
typedef Eigen::Matrix<double, 22, 1> Vector22;
typedef Eigen::Matrix<double, 8, 8>  Matrix8;
typedef Eigen::Matrix<double, 7, 7>  Matrix7;
typedef Eigen::Matrix<double, 6, 6>  Matrix6;          
typedef Eigen::Matrix<double, 5, 5>  Matrix5;         
        
int main(int argc, char** argv)
{
    
  std::cout << "Essential Matrix Estimation with rotation prior!\n";



  double N_iter = 1;
  //set experiment parameters
  double noise = 2.5;
  const size_t n_points = 50;
  double FoV = 100;  // in degrees
  double max_parallax = 2.;  // in meters
  double min_depth = 1.;     // in meters
  double max_depth = 8.;       // in meters
  double outlier_fraction = 0;
  double focal_length = 800; 

   // std::srand(std::time(nullptr));








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
       // str_in.use_custom_t = options.use_custom_t; 
       // str_in.max_rotation = options.max_rotation;                                      
       // str_in.custom_t = options.custom_t; 
                                               
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
      // Rtinit.block<2,2>(0,0) << rotation(0,0), rotation(0,2), rotation(2,0), rotation(2,2);
      // Rtinit.block<3,1>(0,2) << 1, 0, 0;
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
      
      /* 
        result.R_opt = R_opt;
        result.t_opt = t_opt;
        result.E_opt = E_opt;
      */
      Matrix3 E_opt = my_result_prior_pattern.E_opt; 
      Matrix2 R_opt = my_result_prior_pattern.R_opt; 
      Vector3 t_opt = my_result_prior_pattern.t_opt; 
      Matrix3 Ropt = Matrix3::Identity(); 
      Ropt(0,0) = R_opt(0,0); 
      Ropt(0,2) = R_opt(0,1); 
      Ropt(2,0) = R_opt(1,0); 
      Ropt(2,2) = R_opt(1,1);
      Vector3 q_opt = Ropt.transpose() * t_opt;
      
      Matrix6 Q_red;
      RotPrior::createDataCReduced(p0,p1,weights, Q_red); 
      Matrix9 C = Matrix9::Zero(); 
      C.block<6,6>(0,0) = Q_red;
      
      /* Checkin optimality */
      /** Left formulation **/
      std::cout << "########################\nLEFT  FORMULATION\n";
      Vector9 x_left; 
      x_left << E_opt(0,0), E_opt(1,0), E_opt(0,1), E_opt(2,1), E_opt(0,2), E_opt(1,2), t_opt;       
      Vector7 mult_left_init = Vector7::Zero();
       
      std::vector<Matrix9> Aleft; 
      RotPrior::createLeftEConstraints(Aleft);                              
   
          
      SymmCert::SymmCertOptions options_symm;

      
      SymmCert::SymmCertClass<Matrix9, Vector9, Vector7> cert_iter_symm(options_symm);
     
      SymmCert::SymmCertResult<Matrix9, Vector9, Vector7> res_left = cert_iter_symm.getResults(x_left, C, Aleft, mult_left_init, C);
      
      // cert_iter_symm.printResult(res_left); 
      
      Eigen::SelfAdjointEigenSolver<Matrix9> eig_left(res_left.opt_Hessian);
      std::cout << "Eigenvalues Hessian LEFT:\n" << eig_left.eigenvalues() << std::endl; 
      
      const int r_left = 4;
      Eigen::Matrix<double, 9, r_left> Yleft; 
      Matrix9 U = eig_left.eigenvectors();  
      Eigen::Matrix<double, 9, r_left> Ured = U.block<9,r_left>(0,5); 
      Eigen::Matrix4d Dleft = (eig_left.eigenvalues().block<r_left,1>(5,0)).cwiseAbs().cwiseSqrt().asDiagonal();
      Yleft = Ured * Dleft;
      
      StaircaseCert::StaircaseCertOptions options_stairs;
      options_stairs.use_euc_simplification = false;
      options_stairs.max_rank = 9;
      options_stairs.matrixU = U; 
      options_stairs.singularValues = eig_left.eigenvalues().cwiseAbs().cwiseSqrt();
      
                    
                          
                                  
                            
      
      StaircaseCert::StaircaseCertClass<Matrix9, Eigen::Matrix<double, 9, r_left>, Vector9, Vector7> cert_left_stair(options_stairs);
   
      StaircaseCert::StaircaseCertResult<Matrix9, Eigen::Matrix<double, 9, r_left>, Vector7> res_left_stair = cert_left_stair.getResults(x_left, C, Aleft, res_left.opt_mult, Yleft);
        
        
       Eigen::SelfAdjointEigenSolver<Matrix9> eig_left_sol(res_left_stair.opt_Hessian);
      std::cout << "Eigenvalues Hessian LEFT after staircase:\n" << eig_left_sol.eigenvalues() << std::endl; 
      
      
  
  
  
  
  
    
    
    
    
      
      /** Right formulation **/
      std::cout << "########################\nRIGHT FORMULATION\n";
      Vector9 x_right; 
      x_right << E_opt(0,0), E_opt(1,0), E_opt(0,1), E_opt(2,1), E_opt(0,2), E_opt(1,2), q_opt;       
      Vector7 mult_right_init = Vector7::Zero();
       
      std::vector<Matrix9> Aright; 
      RotPrior::createRightEConstraints(Aright);   
     
      SymmCert::SymmCertResult<Matrix9, Vector9, Vector7> res_right = cert_iter_symm.getResults(x_right, C, Aright, mult_right_init, C);
      
      // cert_iter_symm.printResult(res_right); 
      
      Eigen::SelfAdjointEigenSolver<Matrix9> eig_right(res_right.opt_Hessian);
      std::cout << "Eigenvalues Hessian RIGHT:\n" << eig_right.eigenvalues() << std::endl; 
      
      const int r_right = 5;
      Eigen::Matrix<double, 9, r_right> Yright; 
      Matrix9 U_right = eig_right.eigenvectors();  
      Eigen::Matrix<double, 9, r_right> Ured_right = U_right.block<9,r_right>(0,4); 
      Matrix5 Dright = (eig_right.eigenvalues().block<r_right,1>(4,0)).cwiseAbs().cwiseSqrt().asDiagonal();
      Yright = Ured_right * Dright;
      
      StaircaseCert::StaircaseCertOptions options_stairs_R;
      options_stairs_R.use_euc_simplification = true;
      options_stairs_R.max_rank = 9;
      options_stairs_R.matrixU = U_right; 
      options_stairs_R.singularValues = eig_right.eigenvalues().cwiseAbs().cwiseSqrt();
                                                                                                            
      StaircaseCert::StaircaseCertClass<Matrix9, Eigen::Matrix<double, 9, r_right>, Vector9, Vector7> cert_right_stair(options_stairs_R);
   
      StaircaseCert::StaircaseCertResult<Matrix9, Eigen::Matrix<double, 9, r_right>, Vector7> res_right_stair = cert_right_stair.getResults(x_right, 
                                                                                                          C, Aright, res_right.opt_mult, Yright);
      
      
      
      
      
      
      
      
      
      
      /** Both formulation **/
      std::cout << "########################\nBOTH  FORMULATION\n";
      Vector12 x_both; 
      x_both << E_opt(0,0), E_opt(1,0), E_opt(0,1), E_opt(2,1), E_opt(0,2), E_opt(1,2), t_opt, q_opt;       
      Vector13 mult_both_init = Vector13::Zero();
       
      std::vector<Matrix12> Aboth; 
      RotPrior::createBothEConstraints(Aboth);   
      Matrix12 C_ext = Matrix12::Zero(); 
      C_ext.block<6,6>(0,0) = Q_red; 
      
     
      options_symm.verbose = 1; 
      options.estimation_verbose = 1;
      SymmCert::SymmCertClass<Matrix12, Vector12, Vector13> cert_iter_symm_both(options_symm);
     
      SymmCert::SymmCertResult<Matrix12, Vector12, Vector13> res_both = cert_iter_symm_both.getResults(x_both, C_ext, Aboth, mult_both_init, C_ext);
      
      // cert_iter_symm_both.printResult(res_both); 
      
      Eigen::SelfAdjointEigenSolver<Matrix12> eig_both(res_both.opt_Hessian);
      std::cout << "Eigenvalues Hessian BOTH:\n" << eig_both.eigenvalues() << std::endl; 
      
      const int r_both = 7;
      Eigen::Matrix<double, 12, r_both> Yboth; 
      Matrix12 U_both = eig_both.eigenvectors();  
      Eigen::Matrix<double, 12, r_both> Ured_both = U_both.block<12,r_both>(0,5); 
      Matrix7 Dboth = (eig_both.eigenvalues().block<r_both,1>(5,0)).cwiseAbs().cwiseSqrt().asDiagonal();
      Yboth = Ured_both * Dboth;
       
      options_stairs.max_rank = 12;   
      options_stairs.matrixU = U_both; 
      options_stairs.singularValues = eig_both.eigenvalues().cwiseAbs().cwiseSqrt();                                                                                                  
      StaircaseCert::StaircaseCertClass<Matrix12, Eigen::Matrix<double, 12, r_both>, Vector12, Vector13> cert_both_stair(options_stairs);
   
      StaircaseCert::StaircaseCertResult<Matrix12, Eigen::Matrix<double, 12, r_both>, Vector13> res_both_stair = cert_both_stair.getResults(x_both, 
                                                                                                          C_ext, Aboth, res_both.opt_mult, Yboth);
                                                                                                          
                                    
                          
                          
                                    
                      
                      
                                                                                                                                                                                                                
      
      /** Adjugate formulation **/
      std::cout << "########################\nADJUGATE  FORMULATION\n";
      Vector12 x_adj; 
      x_adj << E_opt(0,0), E_opt(1,0), E_opt(0,1), E_opt(2,1), E_opt(0,2), E_opt(1,2), t_opt, q_opt;       
      Vector22 mult_adj_init = Vector22::Zero();
       
      std::vector<Matrix12> Aadj; 
      RotPrior::createAdjugateEConstraints(Aadj);   
     
      SymmCert::SymmCertClass<Matrix12, Vector12, Vector22> cert_iter_symm_adj(options_symm);
      SymmCert::SymmCertResult<Matrix12, Vector12, Vector22> res_adj = cert_iter_symm_adj.getResults(x_adj, C_ext, Aadj, mult_adj_init, C_ext);
      
      // cert_iter_symm_adj.printResult(res_adj); 
      
      Eigen::SelfAdjointEigenSolver<Matrix12> eig_adj(res_adj.opt_Hessian);
      std::cout << "Eigenvalues Hessian ADJ:\n" << eig_adj.eigenvalues() << std::endl;
      
      const int r_adj = 7;
      Eigen::Matrix<double, 12, r_adj> Yadj; 
      Matrix12 U_adj = eig_adj.eigenvectors();  
      Eigen::Matrix<double, 12, r_adj> Ured_adj = U_adj.block<12,r_adj>(0,5); 
      Matrix7 Dadj = (eig_adj.eigenvalues().block<r_adj,1>(5,0)).cwiseAbs().cwiseSqrt().asDiagonal();
      Yadj = Ured_adj * Dadj;
            
      options_stairs.max_rank = 12;                                                                                                      
      options_stairs.matrixU = U_adj; 
      options_stairs.singularValues = eig_adj.eigenvalues().cwiseAbs().cwiseSqrt();
      StaircaseCert::StaircaseCertClass<Matrix12, Eigen::Matrix<double, 12, r_adj>, Vector12, Vector22> cert_adj_stair(options_stairs);
   
      StaircaseCert::StaircaseCertResult<Matrix12, Eigen::Matrix<double, 12, r_adj>, Vector22> res_adj_stair = cert_adj_stair.getResults(x_adj, 
                                                                                                        C_ext, Aadj, res_adj.opt_mult, Yadj);
                                                                                                          
                                                                                                          
                                                                                                          
  return 0;

}
