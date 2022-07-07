#pragma once


/* RotPrior manifold related */
#include "RotPriorTypes.h"
#include "RotPriorManifold.h"
#include "RotPriorUtils.h"


#include <Eigen/Dense>



namespace RotPrior {


struct ProblemCachedMatrices{
EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
        Matrix3 NablaF_Y;
        Matrix3 Mt;
        Matrix9 Mr;
        /// DEFAULT CONSTRUCTOR with default values
        ProblemCachedMatrices(  const Matrix3 & Nabla = Matrix3::Identity(),
                                const Matrix3 &  mt = Matrix3::Identity(),
                                const Matrix9 & mr = Matrix9::Identity()) :
                                NablaF_Y(Nabla), Mt(mt), Mr(mr){}
};


class RotPriorProblem{
public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
        RotPriorProblem(){};  // Default
        
        RotPriorProblem(const Matrix9 & C);  

        // Constructor using two vectors of 3XN corresponding features
        RotPriorProblem(const Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                        const Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                        const Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations);

        // todo
        void setNumberPoints(const Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0) 
        {number_points_ = bearingVectors0.cols();}

                      
        ~RotPriorProblem(){};

         Matrix9 getDataMatrixC(void) {return data_matrix_C_;}
         

         void setPointCorrespondences(const Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                                      const Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                                      const Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations) 
                        {p0_ = bearingVectors0; p1_ = bearingVectors1; w_ = weights_observations; 
                        number_points_ = bearingVectors0.cols();}
   

         void setMatrixPrecon(Matrix3 & matrix_precon) {Matrix_precon_ = matrix_precon;}

         // Pseudo jacobi preconditioner based on the three largest eigenvalues of C
         Matrix3 computePseudoJacobiPrecon(void);


        /// ACCESSORS

        /** Get a const pointer to the SO(3) x S(2) product manifold */
         const RotPriorManifold& getRotPriorManifold() const {
            return domain_;
          }


          void setMr(const Matrix9 & M) {Mr_ = M;};
          void setMt(const Matrix3 & M) {Mt_ = M;};

          Matrix3 getMt(void) {return Mt_; };
          Matrix9 getMr(void) {return Mr_; };
          
          
           /// OPTIMIZATION AND GEOMETRY

          /** Given a matrix Y, this function computes and returns F(Y), the value of
   * the objective evaluated at X */
  double evaluate_objective(const Matrix3 &Y) const;

  double evaluate_objective(const Matrix3 &Y, ProblemCachedMatrices & problem_matrices) const;

    /** Given a matrix Y, this function computes and returns nabla F(Y), the
   * *Euclidean* gradient of F at Y. */

  Matrix3 Euclidean_gradient(const Matrix3 &Y, const ProblemCachedMatrices & problem_matrices) const;

  /** Given a matrix Y in the domain D of the SE-Sync optimization problem and
   * the *Euclidean* gradient nabla F(Y) at Y, this function computes and
   * returns the *Riemannian* gradient grad F(Y) of F at Y */

   Matrix3 Riemannian_gradient(const Matrix3 &Y, const Matrix3 &nablaF_Y) const;

   Matrix3 Riemannian_gradient(const Matrix3 &Y, const ProblemCachedMatrices & problem_matrices) const;


   // Matrix3 Riemannian_gradient(const Matrix3 &Y) const;
   


  /* Preconditioner */
  Matrix3 precondition(const Matrix3& X, const Matrix3 & Xdot) const;

  // Matrix3 preconditionRight(const Matrix3& X, const Matrix3 & Xdot) const;

  /** Given a matrix Y in the domain D of the SE-Sync optimization problem, the
   * *Euclidean* gradient nablaF_Y of F at Y, and a tangent vector dotY in
   * T_D(Y), the tangent space of the domain of the optimization problem at Y,
   * this function computes and returns Hess F(Y)[dotY], the action of the
   * Riemannian Hessian on dotY */

   Matrix3 Riemannian_Hessian_vector_product(const Matrix3 &Y,
                                                   const ProblemCachedMatrices & problem_matrices,
                                                   const Matrix3 &dotY) const;

 

    /** Given a matrix Y in the domain D of the SE-Sync optimization problem and a
  tangent vector dotY in T_Y(E), the tangent space of Y considered as a generic
  matrix, this function computes and returns the orthogonal projection of dotY
  onto T_D(Y), the tangent space of the domain D at Y*/
  Matrix3 tangent_space_projection(const Matrix3 &Y, const Matrix3 &dotY) const;

  /** Given a matrix Y in the domain D of the SE-Sync optimization problem and a
   * tangent vector dotY in T_D(Y), this function returns the point Yplus in D
   * obtained by retracting along dotY */
  Matrix3 retract(const Matrix3 &Y, const Matrix3 &dotY) const;


  // compute Residuals
  Eigen::Matrix<double, 1, Eigen::Dynamic> computeResiduals(const Matrix3 & Rt) const;

  // update points with weights
  void updateWeights(const Eigen::Matrix<double, 1, Eigen::Dynamic> & new_weights);  // This function modifies the weights



private:

  size_t number_points_;

  Eigen::Matrix<double, 3, Eigen::Dynamic> p0_, p1_;
  Eigen::Matrix<double, 1, Eigen::Dynamic> w_;


  /** The product manifold of SO(2) X S(2) ~ E(3) that is the domain of our method */
  RotPriorManifold domain_;

  Matrix9 data_matrix_C_;

  Matrix3 Mt_;

  Matrix9 Mr_;

  Matrix3 Matrix_precon_;

  // Matrices for the mixed term in the hessian
  Matrix9 B1_, B2_, B3_;

  // matrices for the preconditioner: B^T * C * B
  Matrix93 B_;


}; // end of RotPrior problem class




}  // end of RotPrior namespace


