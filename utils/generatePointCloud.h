/*
Original from Opengv.
Adaptation for the (central case) relative pose problem
*/

#pragma once

#include <stdlib.h>
#include <eigen3/Eigen/Dense>
#include <functional>



namespace UtilsRelPose
{

  /** A simple struct that contains the elements of a corresponding fatures as bearing (unit) vectors */
  struct CorrespondingFeatures {
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Eigen::Vector3d bearing_vector_0, bearing_vector_1;

    double weight_;
   /** Simple default constructor; does nothing */
    CorrespondingFeatures() {}

    CorrespondingFeatures(const Eigen::Vector3d & bearing_vector_0,
                          const Eigen::Vector3d & bearing_vector_1,
                          double weight_match = 1.0)
                          : weight_(weight_match), bearing_vector_0(bearing_vector_0),
                          bearing_vector_1(bearing_vector_1) {}

  }; // end of CorrespondingFeatures struct

 // Define types
 typedef std::vector<CorrespondingFeatures> bearing_vectors_t;
 typedef Eigen::VectorXd weights_t;
 
 
                 /* Generate random translation */
                using GenerateTranslation = std::function<Eigen::Vector3d(const double max_parallax,
                                                                          const Eigen::Vector3d& direction_parallax)>;
                /* Generate random rotation */
                using GenerateRotation = std::function<Eigen::Matrix3d(const double max_angle, 
                                                                       const Eigen::Vector3d& rot_angle)>;
        
                  
                /** A simple struct that contains the elements of a corresponding fatures as bearing (unit) vectors */
                struct StrInputGenProblem {
                   EIGEN_MAKE_ALIGNED_OPERATOR_NEW
                   
                   double FoV = 100;                    // in degrees
                   double min_depth = 1.0;              // in meters
                   double max_depth = 8.0;              // in meters 
                   bool allow_coplanar = false;         // allow to generate synthetic scenes with coplanar points
                   double max_X = 20.0;                 // in meters: this value is the absolute value. Max value for X-axis allowed 
                   double focal_length = 800;           // in pixels
                   size_t n_points = 100;
                   double noise = 0.1;                  // in pixels
                   double outlier_fraction = 0.0;
                   double min_epipolar_error_sq = 0.0;  // min epipolar error for the outliers
                   
                   // params for relpose
                   double parallax = 2.0;               // in meters
                   double max_rotation = 0.5;           // in degrees
                   
                   /* Not used */
                   bool use_custom_t = false;           // true if you want to provide t
                   bool use_custom_R = false;           // true if you want to provide R
                   Eigen::Vector3d custom_t; 
                   Eigen::Matrix3d custom_R;  
                   /* Until here */
                   
                   Eigen::Vector3d dir_trans; 
                   Eigen::Vector3d dir_rotation;
                   
                   // constructor
                   StrInputGenProblem(){}; 
                     
                }; // end of StrInputGenProblem struct
                  

                struct StrOutputGenProblem {
                   EIGEN_MAKE_ALIGNED_OPERATOR_NEW
                   
                    Eigen::Vector3d  translation;
                    Eigen::Matrix3d  rotation;
                    
                    bearing_vectors_t  points_correspondences;
                    Eigen::MatrixXd  points_3D;
                    std::vector<int>  indices_outliers;
                    
                    // constructor
                    StrOutputGenProblem(size_t n_points) {
                    translation.setZero(); 
                    rotation.setZero(); 
                    points_correspondences; 
                    points_3D = Eigen::MatrixXd(3, n_points); 
                
                    }; 
                    
                }; // end of StrOutputGenProblem struct
                
                // random rot
                
                Eigen::Matrix3d generatePerturbRotationY(const double max_angle_Y, 
                                                         const Eigen::Vector3d & dir_rotation); 

                Eigen::Vector3d generateRandomTranslationDefault(const double max_parallax, 
                                                                 const Eigen::Vector3d & dir_parallax); 
                                                                 
                Eigen::Matrix3d generateRandomRotationDefault(const double maxAngle, 
                                                              const Eigen::Vector3d & dir_rot); 

                StrOutputGenProblem createSyntheticExperiment(const StrInputGenProblem & in_param, 
                                                        const GenerateTranslation& generateTranslation = generateRandomTranslationDefault,
                                                        const GenerateRotation& generateRotation = generateRandomRotationDefault); 

                                              
                                                                                        

}  // end of namespace UtilsRelPose
