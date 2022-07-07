#pragma once

#include <Eigen/Dense>

#include "RotPriorTypes.h"

/*Define the namespace*/
namespace RotPrior{

  class RotPriorManifold{
  

      public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
        /* Default constructor */
        RotPriorManifold(void){};

        /*Delete each component manifold*/
        ~RotPriorManifold(void){};
        
        /// GEOMETRY
        /** Given a generic matrix A in R^{3 x 3}, this function computes the
        * projection of A onto R (closest point in the Frobenius norm sense).  */
        Matrix3 project(const Matrix3 &A) const;
        
        
        /** Given an element Y in M and a tangent vector V in T_Y(M), this function
       * computes the retraction along V at Y using the QR-based retraction
       * specified in eq. (4.8) of Absil et al.'s  "Optimization Algorithms on
       * Matrix Manifolds").
       */
      Matrix3 retract(const Matrix3 &Y, const Matrix3 &V) const;
            
      /* Projections for each manifold */
      /// Sphere
      Vector3 ProjSphere(const Vector3 &t, const Vector3 &Vt) const;
      /// Rotation  
      Matrix2 ProjRotation(const Matrix2 & R, const Matrix2 & VR) const;
                                    
                               
   /** Given an element Y in M and a matrix V in T_X(R^{p x kn}) (that is, a (p
   * x kn)-dimensional matrix V considered as an element of the tangent space to
   * the *entire* ambient Euclidean space at X), this function computes and
   * returns the projection of V onto T_X(M), the tangent space of M at X (cf.
   * eq. (42) in the SE-Sync tech report).*/
  Matrix3 Proj(const Matrix3 &Y, const Matrix3 &V) const;
    
      };
} /*end of RotPrior namespace*/

