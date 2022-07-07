
#include "./generatePointCloud.h"

#include <vector>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <memory>

#include <Eigen/Core>

#include <functional>


#define PI 3.14159


namespace UtilsRelPose
{


/*void createSyntheticExperiment(
    const size_t n_points,
    double noise,
    double outlier_fraction,
    double FoV,
    double parallax,
    double min_depth,
    double max_depth,
    Eigen::Vector3d & translation,
    Eigen::Matrix3d & rotation,
    bearing_vectors_t & points_correspondences,
    Eigen::MatrixXd & gt,
    std::vector<int> & indices_outliers,
    bool allow_coplanar,   // allow to generate synthetic scenes with coplanar points
    double max_X,  // this value is the absolute value. Max value for X-axis allowed
    double min_epipolar_error_sq,   // min epipolar error for the outliers
    double focal_length    // focal length in pixelss
    )
    */
   

        Eigen::Matrix3d generateRandomRotationYAux(const double maxAngle, 
                                            const  Eigen::Vector3d & dir_rot)
        {
                double angle = maxAngle * (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
                
                Eigen::Matrix3d rot; 
                rot << cos(angle), 0, sin(angle), 
                       0, 1, 0,
                       -sin(angle), 0, cos(angle);
                       
                return rot;
        }


        Eigen::Matrix3d generateRandomRotationAux( double maxAngle )
        {
          // Create rotation as done in OPENGV
          Eigen::Vector3d rpy;
          rpy[0] = ((double) std::rand())/ ((double) RAND_MAX);
          rpy[1] = ((double) std::rand())/ ((double) RAND_MAX);
          rpy[2] = ((double) std::rand())/ ((double) RAND_MAX);

          rpy[0] = maxAngle*2.0*(rpy[0]-0.5);
          rpy[1] = maxAngle*2.0*(rpy[1]-0.5);
          rpy[2] = maxAngle*2.0*(rpy[2]-0.5);

          Eigen::Matrix3d R1;
          R1(0,0) = 1.0;
          R1(0,1) = 0.0;
          R1(0,2) = 0.0;
          R1(1,0) = 0.0;
          R1(1,1) = cos(rpy[0]);
          R1(1,2) = -sin(rpy[0]);
          R1(2,0) = 0.0;
          R1(2,1) = -R1(1,2);
          R1(2,2) = R1(1,1);

          Eigen::Matrix3d R2;
          R2(0,0) = cos(rpy[1]);
          R2(0,1) = 0.0;
          R2(0,2) = sin(rpy[1]);
          R2(1,0) = 0.0;
          R2(1,1) = 1.0;
          R2(1,2) = 0.0;
          R2(2,0) = -R2(0,2);
          R2(2,1) = 0.0;
          R2(2,2) = R2(0,0);

          Eigen::Matrix3d R3;
          R3(0,0) = cos(rpy[2]);
          R3(0,1) = -sin(rpy[2]);
          R3(0,2) = 0.0;
          R3(1,0) =-R3(0,1);
          R3(1,1) = R3(0,0);
          R3(1,2) = 0.0;
          R3(2,0) = 0.0;
          R3(2,1) = 0.0;
          R3(2,2) = 1.0;

          Eigen::Matrix3d rotation = R3 * R2 * R1;

          rotation.col(0) = rotation.col(0) / rotation.col(0).norm();
          rotation.col(2) = rotation.col(0).cross(rotation.col(1));
          rotation.col(2) = rotation.col(2) / rotation.col(2).norm();
          rotation.col(1) = rotation.col(2).cross(rotation.col(0));
          rotation.col(1) = rotation.col(1) / rotation.col(1).norm();
         
          return rotation;
        }



        // Create rotation with noisy IMU         
        Eigen::Matrix3d generatePerturbRotationY(const double max_angle_Y, 
                                        const Eigen::Vector3d & dir_rotation)
                                        
        {
                // NOTE: dir_rotation(0) has the max_angle for the perturbation
                Eigen::Matrix3d rot; 
                rot.setZero(); 
                                
                Eigen::Matrix3d rot_Y = generateRandomRotationYAux(max_angle_Y, dir_rotation);               
                
                Eigen::Matrix3d rot_pert = generateRandomRotationAux( dir_rotation(0) ); 
                
                rot = rot_pert * rot_Y;            
        
                return rot; 
        }
        
        Eigen::Vector3d generateRandomTranslation( double parallax )
        {
          //  std::cout << "Creating translation\n";
          Eigen::Vector3d translation;
          translation[0] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
          translation[1] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
          translation[2] = -(((double) rand())/ ((double) RAND_MAX));

          // std::cout << "Translation created!\n";
          return (parallax * translation.normalized());
        }
        
        
        Eigen::Vector3d generateRandomTranslation(const double parallax , const Eigen::Vector3d & dir_par)
        {

          // std::cout << "Translation created!\n";
          return (parallax * dir_par.normalized());
        }
        
        

        Eigen::Vector3d generateRandomTranslationDefault(const double max_parallax, 
                                                const Eigen::Vector3d & dir_parallax)
        {
                return generateRandomTranslation(max_parallax);
        
        }
        

        
                
        Eigen::Matrix3d generateRandomRotationDefault(const double maxAngle, 
                                                      const Eigen::Vector3d & dir_rot)
        {
                return generateRandomRotationAux( maxAngle );        
        }
        
        
         Eigen::Vector3d addNoiseGaussian(Eigen::Vector3d clean_point, double noise_level)
        {
                Eigen::Vector3d noisy_point = clean_point; 
                double last_entry = noisy_point(2); 
                noisy_point(0) /= last_entry; 
                noisy_point(1) /= last_entry; 
                noisy_point(2) = 1;
                
                
                double x_noise = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0 * noise_level;

                double y_noise = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0 * noise_level;
                
                
                noisy_point(0) += x_noise; 
                noisy_point(1) += y_noise;
                
                return noisy_point; 
                        

        }
        

                
                
                
         Eigen::Vector3d generateRandomPointTruncated(double FoV, double min_depth, 
                                        double max_depth, bool allow_coplanar, double max_X)
                         {

                  Eigen::Vector3d point3D;
                  // Note: tan(X) applies mod(X, pi) before computing the actual tan

                  // depth
                  point3D[2] = min_depth + (max_depth - min_depth) * ( ((double) std::rand() / (double) RAND_MAX));

             
                  double xmax = (tan(FoV * 0.5 * PI / 180)) * point3D[2];
                  if (xmax <= 0) xmax *= -1;  

                           
                  if (!allow_coplanar)
                    // update xmax
                    xmax = std::min(xmax, max_X);

                  point3D[0] = xmax * (((double) rand()) / ((double) RAND_MAX)-0.5)*2.0;
                  point3D[1] = xmax * (((double) rand()) / ((double) RAND_MAX)-0.5)*2.0;

                  return point3D;
                 
                }



StrOutputGenProblem createSyntheticExperiment(const StrInputGenProblem & in_param, const GenerateTranslation& generateTranslation, const GenerateRotation& generateRotation)
{
  // struct for the ouput 
  StrOutputGenProblem str_out(in_param.n_points);
  
  // generate a random pointcloud
  str_out.points_3D.setZero();

  // std::cout << "creating 3D points\n"; 
  
  for( size_t i = 0; i < in_param.n_points; i++ )
        str_out.points_3D.col(i) = generateRandomPointTruncated(in_param.FoV, in_param.min_depth,
                                in_param.max_depth, in_param.allow_coplanar, in_param.max_X );
                                                

  // 2. Create the relative rotation and translation
  // 2.1 Fix first frame
  Eigen::Matrix3d R1 = Eigen::Matrix3d::Identity();
  Eigen::Vector3d T1 = Eigen::Vector3d::Zero();

  // 2.2 Fix second frame
  Eigen::Matrix3d rotation = Eigen::Matrix3d::Zero();
  Eigen::Vector3d T2 = Eigen::Vector3d::Zero();

  // std::cout << " Creating second camera pose\n";
  bool valid_position = false;
  Eigen::MatrixXd points3D_frame2(3, in_param.n_points); 
  points3D_frame2.setZero();
  
  int max_n_iters = 200, n_iters = 0;

  
  
  do
  {
    // create random translation
    // note: we always create a backward movement

    T2 = generateTranslation(in_param.parallax, in_param.dir_trans);


    // do the sae with rotation 
    rotation = generateRotation(in_param.max_rotation, in_param.dir_rotation);


    // compute relative coordinates
    points3D_frame2.setZero();
    for( size_t i = 0; i < in_param.n_points; i++ )
    {
      Eigen::Vector3d temp = rotation.transpose() * (str_out.points_3D.col(i)- T2);
      points3D_frame2.col(i) = Eigen::Map<Eigen::Vector3d>(temp.data(), 3, 1);
    }

    
    // check condition for FoV
    Eigen::VectorXd ratio_X = (points3D_frame2.row(0)).cwiseQuotient((points3D_frame2.row(2)));

    double max_ratio_x = (ratio_X.array().abs()).maxCoeff();


    // check if any point was outside the FoV

    if (max_ratio_x > tan(in_param.FoV * 0.5 * PI / 180))
    {
        std::cout << "At least one point did not fulfill the condition\n"; 
        n_iters++;
        continue;
    }


    Eigen::VectorXd ratio_Y = (points3D_frame2.row(1)).cwiseQuotient((points3D_frame2.row(2)));
 
    double max_ratio_y = (ratio_Y.array().abs()).maxCoeff();
    // check if any point was outside the FoV

    if (max_ratio_y > tan(in_param.FoV * 0.5 * PI / 180))
    {
        std::cout << "At least one point did not fulfill the condition\n"; 
        n_iters++;
        continue;
    }

    // if we have arrive here, the position is valid
    valid_position = true;

  }while(!valid_position && (n_iters <= max_n_iters));

  // check invalid rotation 
  if ((!valid_position) && (n_iters > max_n_iters))
        std::cout << "[ERROR] Pose is not valid.\nSome points do not lie in the FOV\n";
  
  
  // save pose
  str_out.translation = T2; 
  str_out.rotation = rotation;


  // Generate the correspondences
  Eigen::Matrix3d K; 
  K.setZero();
  K(0, 0) = in_param.focal_length; 
  K(1, 1) = in_param.focal_length; 
  K(0, 2) = tan(in_param.FoV * 0.5 * PI / 180) * in_param.focal_length; 
  K(1, 2) = tan(in_param.FoV * 0.5 * PI / 180) * in_param.focal_length; 
  K(2, 2) = 1; 
  
  Eigen::Matrix3d invK = K.inverse(); 
  
  

  for( size_t i = 0; i < in_param.n_points; i++ )
  {
    Eigen::Vector3d obs1, obs2, o1, o2;
    obs1 = str_out.points_3D.col(i);
    obs2 = points3D_frame2.col(i);
        
    double last_o1 = obs1(2), last_o2 = obs2(2); 
    
    
    // observations in hom. coord.
    o1 << obs1(0) / last_o1, obs1(1) / last_o1, 1;  
    o2 << obs2(0) / last_o2, obs2(1) / last_o2, 1;  
      
    //add noise
    obs1 = K * o1;
    obs2 = K * o2;
    
    if(in_param.noise > 0.0 )
    {
      obs1 = addNoiseGaussian(obs1, in_param.noise);
      obs2 = addNoiseGaussian(obs2, in_param.noise);
    }

    CorrespondingFeatures correspondence;
    correspondence.bearing_vector_0 = (invK * obs1).normalized();
    correspondence.bearing_vector_1 = (invK * obs2).normalized(); 
    correspondence.weight_ = 1.0;

    // add it to the std::vector
    str_out.points_correspondences.push_back(correspondence);
  }

  return str_out; 

}


}  // end of namespace UtilsRelPose
