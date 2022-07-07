 #pragma once 
 
 
 #include "RotPriorTypes.h"
 
 
 namespace RotPrior{
   
      Matrix3 generatePerturbRotationY(const double max_angle_Y, 
                                        const Vector3 & dir_rotation);
      Matrix3 computeEfromRt(const Matrix3 & R, const Vector3 & t);
      
      void computeRtfromE( const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                              const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                              const   Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations, 
                              const Matrix3& E, Matrix3 & R, Vector3 & t );
                                                                 
      Matrix34 convertRt2ExpandRt(const Matrix3 & Rt, double pad_int = 1.);
       
      Matrix3 projectEtoPrior(  const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                                const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                                const   Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations, 
                                const   Matrix3& E, Matrix2 & rot, Vector3 & trans);
      
      
      Matrix2 projectRtoPrior(const Matrix3 & R); 
      
      Matrix3 computeEfromRtPrior(const Matrix2 & rot, const Vector3 trans); 
      
      Matrix3 computeEfromRtPrior(const Matrix3 & rt); 
      
      Matrix2 obtainRot(const Matrix2 & R); 
      
      
      void convertVector2Matrices(Matrix2 & rot, Vector3 & trans, const Matrix3 & x);
   
      Matrix3 SymProduct(const Matrix3 & Ydot, const Matrix3 & Y, const Matrix3 & nabla_Y);
       
       
      bool computeInitPattern(const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                              const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                              const   Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations, 
                              Matrix3& E, Matrix2 & rot, Vector3& trans);
      
        
      void createDataCReduced(const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                              const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                              const   Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations, 
                              Matrix6& Q_red);
      
      void createDataCReduced(Matrix9& C, Matrix6& Q_red);
      
      
      Matrix39 constructMixedDifTtr(const Matrix9 & B1, const Matrix9 & B2,
                const Matrix9 & B3,const Matrix9 & C , const Vector3 & t, Matrix3 & R);
                
      void constructBMatricesForMixedDiftR(Matrix9& B1, Matrix9 & B2, Matrix9 & B3);
      Matrix3 projectToEssentialManifold(const Matrix3 & E_hat);
                                 
      Matrix9 createMatrixR(const Matrix9 & C, const Vector3 & t);
      Matrix3 createMatrixT(const Matrix9 & C, const Matrix3 & R);
      Matrix9 createThat(const Vector3 & t);
      Matrix39 createRhat(const Matrix3& R);
      Eigen::MatrixXf skew(const Eigen::MatrixXf & M);
      Eigen::MatrixXf symm(const Eigen::MatrixXf & M);
      Matrix3 cross(const Vector3 & t);
      Vector9 vec(Matrix3 & M);
      
      Matrix3 computeEfromRt(const Matrix34& Rt);
      
      Matrix9 constructDataMatrix(const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                                   const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                                   const   Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations);
      
}  // end of namespace RotPrior
