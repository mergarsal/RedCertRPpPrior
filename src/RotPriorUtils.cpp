#include "RotPriorUtils.h"


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <Eigen/SVD> 
#include <Eigen/Core>
#include <Eigen/Geometry> 


#define SQRT2 1.41421356237
#define R180PI 57.2957795131

namespace RotPrior
{

        Matrix3 computeEfromRt(const Matrix3 & R, const Vector3 & t)
        {
            Matrix3 t_skew = cross(t);
            Matrix3 E = cross(t) * R;
            return E / (E.norm()) * SQRT2;
        }

        Matrix3 computeEfromRt(const Matrix34& Rt)
        {
            Vector3 t = Rt.block<3, 1>(0, 3);
            Matrix3 R = Rt.block<3, 3>(0, 0);
            return (computeEfromRt(R, t));
        }

         void computeRtfromE( const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                              const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                              const   Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations, 
                              const Matrix3& E, Matrix3 & R, Vector3 & t )
         {
            // Note: E should be already in \Me

            Eigen::JacobiSVD<Matrix3> svd(E, Eigen::ComputeFullU | Eigen::ComputeFullV);

            Matrix3 W;
            W.setZero();
            W(0, 1) = -1; W(1, 0) = 1; W(2, 2) = 1;

            Matrix3 R1, R2, V, U;
            U = svd.matrixU();
            V = svd.matrixV();

            // Rotation matrix
            R1 = U * W * V.transpose();

            R2 = U * W.transpose() * V.transpose();

            // for R1, R2 in SO(3)
            if (R1.determinant() < 0)       R1 = -R1;
            if (R2.determinant() < 0)       R2 = -R2;


            // Translation vector
            Vector3 t1;
            t1 = U.col(2);


            // Cheirality

            unsigned int N = weights_observations.cols();
            Eigen::MatrixXd X(9, N);
            Eigen::MatrixXd Y(6, N);
            X.setZero();
            Y.setZero();


            for (int i = 0; i < N; i++)
            {
                Vector9 temp;
                Vector3 v1 = bearingVectors1.col(i);
                Vector3 v0 = bearingVectors0.col(i);
                double weight_i = weights_observations[i];

                temp.setZero();
                for (int j = 0; j < 3; j++)    temp.block<3, 1>(j*3, 0) = weight_i * v1[j] * v0;
                X.col(i) = temp;


                Y.block<3, 1>(0, i) = v0 * v1.norm(); Y.block<3, 1>(3, i) = v1 * v0.norm();

            }

            Matrix3 ER1 = E * E.transpose() * R1;
            Matrix3 ER2 = E * E.transpose() * R2;

            Eigen::MatrixXd M11(N, 1);
            Eigen::MatrixXd M12(N, 1);
            M11 = X.transpose() * Eigen::Map<Eigen::Matrix<double, 9, 1> > (ER1.data(), 9, 1);
            M12 = X.transpose() * Eigen::Map<Eigen::Matrix<double, 9, 1> > (ER2.data(), 9, 1);
            double n_non_zeros_M11 = (M11.array() > 0).count();
            double n_non_zeros_M12 = (M12.array() > 0).count();


            if (n_non_zeros_M11 >= n_non_zeros_M12)      R = R1;
            else                                         R = R2;


            // Compute the right translation

            Eigen::MatrixXd M21(N, 1), M22(N, 1);
            Eigen::MatrixXd R_eye(6, 3); R_eye.block<3,3>(0,0) = - R.transpose(); R_eye.block<3,3>(3, 0) = Matrix3::Identity();


            M21 = Y.transpose() * R_eye * t1;

            double n_non_zeros_M21 = (M21.array() > 0).count();
            double n_non_zeros_M22 = (-M21.array() > 0).count();


            if (n_non_zeros_M21 >= n_non_zeros_M22)    t = t1;
            else                                       t = -t1;

        };





        Matrix34 convertRt2ExpandRt(const Matrix3 & Rt, double pad_int)
        {
                Matrix34 Rt_exp; 
                   
                
                Rt_exp.setZero(); 
                Rt_exp(0,0)=Rt(0,0); 
                Rt_exp(0,2)=Rt(0,1);
                Rt_exp(2,0)=Rt(1,0);
                Rt_exp(2,2)=Rt(1,1);                 
                Rt_exp(1, 1) = pad_int;
                Rt_exp.block<3, 1>(0, 3) = Rt.block<3, 1>(0, 2);
        
                
                return Rt_exp;
        }


        Matrix3 projectEtoPrior(const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                                const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                                const   Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations, 
                                const   Matrix3& E, Matrix2 & rot, Vector3 & trans)
        
        {
                // This function projects a general E 
                // to a Ez
                // 1. Get rotation and trans from E
                Matrix3 rot_gen; 
                // Vector3 trans;
                computeRtfromE(bearingVectors0, bearingVectors1, weights_observations, E, rot_gen, trans);
                
                // 2. Project rot to SO(2)
                rot = projectRtoPrior(rot_gen); 
                
                // 3. Recompute essential matrix
                Matrix3 E_p = computeEfromRtPrior(rot, trans);       
        
                return E_p; 
        }

        Matrix2 projectRtoPrior(const Matrix3 & R)
        {
                Eigen::Quaternion<double> eq = Eigen::Quaternion<double>(R); 
      
                      
                Vector4 init_rot = Vector4(); 
                init_rot << eq.w(), -eq.x(), -eq.y(), -eq.z();
                
                
                // let be a Z-rotation
                init_rot(1) = 0; 
                init_rot(3) = 0; 
                init_rot.normalize();  
                
                // NOTE: use conjugate 
                Eigen::Quaternion<double> rot_quat = Eigen::Quaternion<double>(init_rot[0], -init_rot[1], -init_rot[2], -init_rot[3]); 
                Matrix3 rot = rot_quat.toRotationMatrix(); 
                
                Matrix2 rot_r = Matrix2::Zero(); 
                rot_r(0,0)=rot(0,0);
                rot_r(0,1)=rot(0,2); 
                rot_r(1,0)=rot(2,0);
                rot_r(1,1)=rot(2,2);
                // convert back to rotation matrix                          
                return (rot_r); 
        }


        Matrix3 computeEfromRtPrior(const Matrix2 & rot, const Vector3 trans)
        {
                Matrix3 rt; 
                rt.block<2, 2>(0, 0) = rot; 
                rt.block<3, 1>(0, 2) = trans;
                return (computeEfromRtPrior(rt));
        }

        Matrix3 computeEfromRtPrior(const Matrix3 & rt)
        {
                return (computeEfromRt(convertRt2ExpandRt(rt, 1.0)));
        }
        
        

        Matrix2 obtainRot(const Matrix2 & R)
            {
                Matrix2 R_p;
                Eigen::JacobiSVD<Matrix2> svd(R, Eigen::ComputeFullU | Eigen::ComputeFullV);

                double detU = svd.matrixU().determinant();
                double detV = svd.matrixV().determinant();

                if (detU * detV > 0)
                    R_p =  svd.matrixU() * svd.matrixV().transpose();
                else
                {
                    Matrix2 Uprime = svd.matrixU();
                    Uprime.col(Uprime.cols() - 1) *= -1;
                    R_p = Uprime * svd.matrixV().transpose();
                }
            
                return R_p;
            
            }

        

         //template <typename M>
         Matrix3 SymProduct(const Matrix3 & Ydot, const Matrix3 & Y, const Matrix3 & nabla_Y) //const
         {
                Matrix3 P = Y.transpose() * nabla_Y;
                
                return (Ydot * 0.5 * (P + P.transpose()));
         }
         
        
         void createDataCReduced(const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                                 const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                                 const   Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations, 
                                 Matrix6& Q_red)
         {
               Matrix9 C = constructDataMatrix(bearingVectors0, bearingVectors1, weights_observations);
               // normalize C 
               C /= bearingVectors0.cols(); 
                
               // auxiliar matrix 
               Matrix9 Q = C; 
               
               // remove e00 = e11
               Q.col(0) += C.col(8);
               Q.row(0) += Q.row(8); 
               
               // remove e01 = -e10
               Q.col(6) -= Q.col(2); 
               Q.row(6) -= Q.row(2); 
               
               
               Q_red.setZero(); 
                
               int ki=0,kj=0;
               // we remove 2, 4 and 8-th
               for(int i=0;i<=7;i++)
               {
                if ((i==2) || (i==4))
                        continue;
                kj=ki;
                for(int j=i;j<=7;j++)
                {
                        if ((j==2) || (j==4))
                                continue;
                          
                        // else
                        Q_red(ki,kj) = Q(i,j);
                        Q_red(kj,ki) = Q(j,i);
                        kj++;
                }
                ki++;
               }
         
               Q_red /= Q_red.trace();
                
               return;      
         
         }
         
        Vector9 vec(Matrix3 & M)
        {    return Eigen::Map<Vector9> (M.data(), 9, 1);       }

        Matrix3 cross(const Vector3 & t)
        {
                Matrix3 t_cross;
                t_cross.setZero();
                t_cross(0, 1) = -t(2); t_cross(0, 2) = t(1);
                t_cross(1, 0) = t(2); t_cross(1, 2) = -t(0);
                t_cross(2, 0) = -t(1); t_cross(2, 1) = t(0);
                return t_cross;
        }
         
        Eigen::MatrixXf skew(const Eigen::MatrixXf & M)
        {
            return (0.5 * (M - M.transpose()));
        }

        Eigen::MatrixXf symm(const Eigen::MatrixXf & M)
        {
            return (0.5 * (M + M.transpose()));
        }
        
        Matrix39 createRhat(const Matrix3& R)
        {
            Matrix39 R_hat;
            R_hat.setZero();

            for (int i = 0; i < 3; i++)
            {
                Vector3 column_i = R.col(i);
                Matrix3 temp_cross = cross(-column_i);
                R_hat.block<3, 3>(0, i*3) = temp_cross.transpose();  

            }
            return R_hat;
        }



        Matrix9 createThat(const Vector3 & t)
        {
                Matrix9 that; that.setZero();
                Matrix3 temp_cross;
                temp_cross.setZero();
                temp_cross = cross(t);

                for (int i = 0; i < 3; i++) that.block<3, 3>(i*3, i*3) = temp_cross;
                return (that);
        }

        Matrix3 createMatrixT(const Matrix9 & C, const Matrix3 & R)
        {
                const Matrix39 R_hat = createRhat(R) ;
                const Matrix39 temp = R_hat * C;
                const Matrix3 temp_final = temp * (R_hat.transpose());

                return (0.5 * (temp_final + temp_final.transpose()));
        }

        Matrix9 createMatrixR(const Matrix9 & C, const Vector3 & t)
        {

                Matrix9 that = createThat(t);
                Matrix9 temp = (that.transpose() * (C * that));
                return 0.5 * (temp + temp.transpose());

        }
        

         Matrix9 constructDataMatrix(const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                                     const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                                     const   Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations)
        {

                Matrix9 C;
                Matrix9 temp;
                // clean output matrix
                C.setZero();

                for (int i = 0; i < bearingVectors0.cols(); i++)
                {
                        // clean data
                        temp.setZero();
                        const Vector3 v1 = bearingVectors1.col(i);
                        const Vector3 v0 = bearingVectors0.col(i);
                        const double weight = weights_observations[i];
                        for (int j = 0; j < 3; j++)    temp.block<3, 1>(j*3, 1) = v1[j] * v0;

                        C += weight * temp * temp.transpose();
                }
               
                return 0.5 * (C + C.transpose());
        }
        
        
         void createDataCReduced(Matrix9& C, Matrix6& Q_red)
         {
                
               // auxiliar matrix 

               Matrix9 Q = C / C.trace(); 
               
               // remove e00 = e11
               Q.col(0) += C.col(8);
               Q.row(0) += Q.row(8); 
               
               // remove e01 = -e10
               Q.col(6) -= Q.col(2); 
               Q.row(6) -= Q.row(2); 
               

               
                Q_red.setZero(); 
                
               int ki=0,kj=0;
               // we remove 2, 4 and 8-th
               for(int i=0;i<=7;i++)
               {
                if ((i==2) || (i==4))
                        continue;
                kj=ki;
                for(int j=i;j<=7;j++)
                {
                        if ((j==2) || (j==4))
                                continue;
                          

                        // else
                        Q_red(ki,kj) = Q(i,j);
                        Q_red(kj,ki) = Q(j,i);
                        kj++;
                }
                ki++;
               }
         

               Q_red /= Q_red.trace();
               
               return;      
         
         }
         
         /* initialization based on DLT + pattern */ 
    
         bool computeInitPattern(const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                                 const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                                 const   Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations, 
                                 Matrix3& E, Matrix2 & rot, Vector3& trans)
         {

               
              
               Matrix6 Q_red; 
               Q_red.setZero(); 
               createDataCReduced(bearingVectors0, bearingVectors1, weights_observations, Q_red);
               

               
               // compute init
               Eigen::JacobiSVD<Matrix6> svd(Q_red, Eigen::ComputeFullU | Eigen::ComputeFullV);

               Vector6 sigmas = svd.singularValues();
               
               Vector6 vecE = svd.matrixV().col(5);
               

               
               // check how well the problm instances is conditioned
               if ((sigmas(4) < 1e-06) || (sigmas(3) < 5e-04) )
                return false;
               

               // fill the other entries
               Vector9 exp_vecE; 
               exp_vecE.setZero(); 
               exp_vecE << vecE(0),vecE(1),-vecE(4),vecE(2),0,vecE(3),vecE(4),vecE(5),vecE(0);
               
            
               
               // project to essential set               
               Matrix3 E_initial = Eigen::Map<Matrix3>(exp_vecE.data(), 3, 3);

               
               
               E_initial = projectToEssentialManifold(E_initial);

               
               // project to prior 
               E = projectEtoPrior(bearingVectors0, bearingVectors1, weights_observations, E_initial, rot, trans); 

               
               return true;            
                
         }
            
         Matrix3 projectToEssentialManifold(const Matrix3 & E_hat)
         {
            const Eigen::JacobiSVD<Matrix3> svd(E_hat, Eigen::ComputeFullU | Eigen::ComputeFullV);
            
            const Eigen::DiagonalMatrix<double, 3> d(1, 1, 0);  // imply: ||E||_f ^2 = 2

            // check determinant
            Matrix3 Ur = svd.matrixU();
            Matrix3 Vrt = svd.matrixV().transpose();

            return (Ur * d * Vrt);
         }


        void constructBMatricesForMixedDiftR(Matrix9& B1, Matrix9 & B2, Matrix9 & B3)
        {
            // useful vars
            size_t r1 = 0, r4 =1, r7 = 2, r2 = 3, r5 = 4, r8 = 5, r3 = 6, r6 = 7, r9 = 8;

            // create first matrix
            B1.setZero();
            B1(1, r7) = -1; B1(2, r4) = 1; B1(4, r8) = -1; B1(5, r5) = 1; B1(7, r9) = -1; B1(8, r6) = 1;
            B1.transposeInPlace();

            // create second matrix
            B2.setZero();
            B2(0, r7) = 1; B2(2, r1) = -1; B2(3, r8) = 1; B2(5, r2) = -1; B2(6, r9) = 1; B2(8, r3) = -1;
            B2.transposeInPlace();

            // create third matrix
            B3.setZero();
            B3(0, r4) = -1; B3(1, r1) = 1; B3(3, r5) = -1; B3(4, r2) = 1; B3(6, r6) = -1; B3(7, r3) = 1;
            B3.transposeInPlace();

        }
        
        
         Matrix39 constructMixedDifTtr(const Matrix9 & B1, const Matrix9 & B2,
                const Matrix9 & B3,const Matrix9 & C , const Vector3 & t, Matrix3 & R)
        {
            Matrix39 mixed_t;
            Matrix9 temp_coeff;
            Vector9 term_ii;


            Matrix9 that = createThat(t);
            
            // NOTE: avoid the loop and define it directly
            // clean everything
            term_ii.setZero();
            mixed_t.setZero();
            temp_coeff.setZero();

            // first row
            temp_coeff = B1 * C * that;
            term_ii = ( 0.5 * (temp_coeff + temp_coeff.transpose())) * vec(R);
            mixed_t.block<1, 9>(0, 0) = term_ii.transpose();

            // second row
            term_ii.setZero();
            temp_coeff.setZero();
            temp_coeff = B2 * C * that;
            term_ii = ( 0.5 * (temp_coeff + temp_coeff.transpose())) * vec(R);
            mixed_t.block<1, 9>(1, 0) = term_ii.transpose();

            // third row
            term_ii.setZero();
            temp_coeff.setZero();
            temp_coeff = B3 * C * that;
            term_ii = ( 0.5 * (temp_coeff + temp_coeff.transpose())) * vec(R);
            mixed_t.block<1, 9>(2, 0) = term_ii.transpose();

           return mixed_t;
        }


        
        
}  // end of namespace RotPrior




