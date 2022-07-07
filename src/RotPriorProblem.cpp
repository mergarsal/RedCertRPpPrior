#include "RotPriorProblem.h"


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>


namespace RotPrior{
        RotPriorProblem::RotPriorProblem(const Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                                         const Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                                         const Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations) 
        {

                p0_ = bearingVectors0; 
                p1_ = bearingVectors1; 
                w_ = weights_observations; 
                
                data_matrix_C_ = constructDataMatrix(bearingVectors0, bearingVectors1, weights_observations);
                number_points_ = weights_observations.cols();

                Mt_.setZero();
                Mr_.setZero();
                constructBMatricesForMixedDiftR(B1_, B2_, B3_);


        }; //end of constructor
        RotPriorProblem::RotPriorProblem(const Matrix9& data_C)
        {
                data_matrix_C_ = data_C;

                Mt_.setZero();
                Mr_.setZero();
                constructBMatricesForMixedDiftR(B1_, B2_, B3_);

        }; //end of constructor

         
        
        Eigen::Matrix<double, 1, Eigen::Dynamic> RotPriorProblem::computeResiduals(const Matrix3 & Rt) const
        {
            Eigen::Matrix<double, 1, Eigen::Dynamic> residual;
            residual.resize(number_points_);

            Matrix3 E = computeEfromRt(convertRt2ExpandRt(Rt, 1.));

            for(size_t i = 0; i < number_points_; i++)
            {
                    Vector3 v1 = p1_.col(i);
                    Vector3 v0 = p0_.col(i);
                    residual(i) = (pow((v0.transpose() * E * v1), 2));
            }
    
            return residual;
        }        
        
        
        // update points with weights
        void RotPriorProblem::updateWeights(const Eigen::Matrix<double, 1, Eigen::Dynamic> & new_weights)
        {

            assert(new_weights.cols() == number_points_);
            for(size_t i = 0; i < number_points_; i++)   w_[i] = new_weights(i);


            // update data matrix
            data_matrix_C_ = constructDataMatrix(p0_, p1_, w_);
        }

        

        // Compute pseudo Jacobi preconditioner
        Matrix3 RotPriorProblem::computePseudoJacobiPrecon(void)
        {
                // If you use the 8 points algorithm, please
                // select the use_precon option

                Eigen::JacobiSVD<Matrix9> svd(data_matrix_C_, Eigen::ComputeFullU | Eigen::ComputeFullV);  //faster

                Vector9 svd_diag = svd.singularValues();
               
                return (Matrix3::Identity() / svd_diag(0));
        }

        

      
        
        // apply precon to the full X
        Matrix3 RotPriorProblem::precondition(const Matrix3& X, const Matrix3 & Xdot) const
        {
            return tangent_space_projection(X,  Matrix_precon_ * Xdot);
        }
       
        
       


        double RotPriorProblem::evaluate_objective(const Matrix3 &Y) const {
            // save this for latter (gradient & hessian)
            Vector3 t; 
            Matrix3 R;
            Matrix2 rot; 
            t = Y.block<3,1>(0,2);
            rot=Y.block<2,2>(0,0); 
            // expand matrix
            R.setZero();
            R(0,0)=rot(0,0); 
            R(0,2)=rot(0,1); 
            R(2,0)=rot(1,0);
            R(2,2)=rot(1,1);
            R(1,1) = 1;            
                        
            Matrix3 Mt = createMatrixT(data_matrix_C_, R);
            Matrix9 Mr = createMatrixR(data_matrix_C_, t);

            return (0.5 * (t.transpose() * Mt * t).trace());

        }

         double RotPriorProblem::evaluate_objective(const Matrix3 &Y, ProblemCachedMatrices & problem_matrices) const {

            Vector3 t; 
            Matrix3 R;
            Matrix2 rot; 
            t = Y.block<3,1>(0,2);
            rot=Y.block<2,2>(0,0); 
            // expand matrix
            R.setZero();
            R(0,0)=rot(0,0); 
            R(0,2)=rot(0,1); 
            R(2,0)=rot(1,0);
            R(2,2)=rot(1,1);
            R(1,1) = 1; 
            
        
            problem_matrices.Mt = createMatrixT(data_matrix_C_, R);
            problem_matrices.Mr = createMatrixR(data_matrix_C_, t);
            
                      
            
            return (0.5 * (t.transpose() * problem_matrices.Mt * t).trace());
        }



             Matrix3 RotPriorProblem::Euclidean_gradient(const Matrix3 &Y, const ProblemCachedMatrices & problem_matrices) const
            {
                   Vector3 t; 
                   Matrix3 R;
                   Matrix2 rot; 
                    t = Y.block<3,1>(0,2);
                    rot=Y.block<2,2>(0,0); 
                    // expand matrix
                    R.setZero();
                    R(0,0)=rot(0,0); 
                    R(0,2)=rot(0,1); 
                    R(2,0)=rot(1,0);
                    R(2,2)=rot(1,1);
                    R(1,1) = 1; 
                   
           
                Matrix3 G;
                G.setZero();
                
   
                
                Vector9 mr = problem_matrices.Mr * vec(R);
                
                Matrix3 Gr = Eigen::Map<Matrix3> (mr.data(), 3, 3);
                G(0,0)=Gr(0,0);
                G(0,1)=Gr(0,2);
                G(1,0)=Gr(2,0);
                G(1,1)=Gr(2,2);
                
                G.block<3,1>(0,2) = problem_matrices.Mt * t;
                return G;
            }

             Matrix3 RotPriorProblem::Riemannian_gradient(const Matrix3 &Y, const Matrix3 &nablaF_Y) const
             {
              return tangent_space_projection(Y, nablaF_Y);
             }
             

            Matrix3 RotPriorProblem::Riemannian_gradient(const Matrix3 &Y, const ProblemCachedMatrices & problem_matrices) const {
            
            Matrix3 Gg = Euclidean_gradient(Y, problem_matrices); 
            
              return tangent_space_projection(Y, Gg);
            }


       /** Given a matrix Y in the domain D of the SE-Sync optimization problem, and
           * a tangent vector dotY in T_D(Y), the tangent space of the domain of the
           * optimization problem at Y, this function computes and returns Hess
           * F(Y)[dotY], the action of the Riemannian Hessian on dotX */                                               

     Matrix3 RotPriorProblem::Riemannian_Hessian_vector_product(const Matrix3 &Y,
                                                   const ProblemCachedMatrices & problem_matrices,
                                                   const Matrix3 &dotY) const
                                                   {
                                                   // Euclidean Hessian-vector product
                                                    Matrix3 HessRiemannian;
                                                    HessRiemannian.setZero();


                                                    
                                                    Vector3 t = Y.block<3, 1>(0, 2); 
                                                    Matrix3 R;
                                                    Matrix2 rot = Y.block<2, 2>(0, 0); 

                                                    // expand matrix
                                                    R.setZero();
                                                    R(0,0)=rot(0,0); 
                                                    R(0,2)=rot(0,1); 
                                                    R(2,0)=rot(1,0);
                                                    R(2,2)=rot(1,1);
                                                    R(1,1) = 1; 
                                                    
                                                                                                        
                                                    Matrix34 dotYext = convertRt2ExpandRt(dotY, 0.0);
                                                    
                                                    Vector3 Vt = dotYext.block<3, 1>(0, 3);
                                                    Matrix3 VR = dotYext.block<3, 3>(0, 0);
                                                    
                                                    
                                                    Matrix39 mixed_terms_der = constructMixedDifTtr(B1_,
                                                                                                B2_, B3_, data_matrix_C_ , t, R);


                                                    // Compute Euclidean Hessian
                                                    Matrix9By12 coeff_R;
                                                    coeff_R.setZero(); 
                                                    coeff_R.block<9, 9>(0, 0) = problem_matrices.Mr;
                                                    coeff_R.block<9, 3>(0, 9) = 2 * mixed_terms_der.transpose();

                                                    Vector12 Vx;
                                                    Vx.setZero();
                                                    Vx.block<9, 1>(0, 0) = vec(VR);
                                                    Vx.block<3, 1>(9, 0) = Vt;

                                                    Vector9 mr = coeff_R * Vx;
                                                    Matrix3 HessR = Eigen::Map<Matrix3> (mr.data(), 3, 3);

                                                    Matrix3By12 coeff_T;
                                                    coeff_T.setZero(); 
                                                    coeff_T.block<3, 3>(0, 9) = problem_matrices.Mt;
                                                    coeff_T.block<3, 9>(0, 0) = 2 * mixed_terms_der;

                                                    Vector3 HessT = coeff_T * Vx;

                                                    // recover NabldaF(Y)
                                                    Vector3 nabla_Yt = (problem_matrices.NablaF_Y).block<3,1>(0,3);
                                                    Matrix3 nabla_YR = (problem_matrices.NablaF_Y).block<3,3>(0,0);

                                                    // clean
                                                    HessRiemannian.setZero();

                                                    // Riemannain Hessian for t (sphere)
                                                    HessRiemannian.block<3,1>(0, 2) = domain_.ProjSphere(t, HessT - (t.dot(nabla_Yt)) * Vt);

                                                    // Riemannain Hessian for R (rotation)
                                                                                        
                                                    // Riemannain Hessian for R (rotation)
                                                    Matrix3 HessRR = HessR - SymProduct(VR, R, nabla_YR);
                                                    
                                                    Matrix2 HessRred = Matrix2::Zero(); 
                                                    HessRred(0,0)=HessRR(0,0); 
                                                    HessRred(0,1)=HessRR(0,2);
                                                    HessRred(1,0)=HessRR(2,0);
                                                    HessRred(1,1)=HessRR(2,2);
                                                    HessRiemannian.block<2, 2>(0, 0) = domain_.ProjRotation(rot, HessRred);
                                                    return HessRiemannian;
                                                    
                                                  
                           }


     

          Matrix3 RotPriorProblem::tangent_space_projection(const Matrix3 &Y,
                                                            const Matrix3 &dotY) const 
                                                            { return domain_.Proj(Y, dotY); }


            Matrix3 RotPriorProblem::retract(const Matrix3 &Y, const Matrix3 &dotY) const
            {
                return domain_.retract(Y, dotY);

            }

      

} // end of Essential namespace
