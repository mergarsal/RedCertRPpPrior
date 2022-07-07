#include "RotPriorManifold.h" 

#include "RotPriorUtils.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace RotPrior{


            Matrix3 RotPriorManifold::project(const Matrix3 & P) const
            {
                 Matrix3 A;
                 
                 A.setZero();
                 // 1. Project translation 
                 A.block<3, 1>(0, 2) = (P.block<3, 1>(0, 2)).normalized();
                 
                 // 2. Project rotation
                 
                 Matrix2 R = P.block<2,2>(0, 0);  

                // fill the output matrix
                // 
                A.block<2,2>(0,0) = obtainRot(R); 
                 

                return A;
            }
            
            
           
            
            Vector3 RotPriorManifold::ProjSphere(const Vector3 &t, const Vector3 &Vt) const
            { Vector3 tret = Vt - t.dot(Vt) * t;
            return  tret; }


            /* ROTATION */
            /** 
            Matrix2 RotPriorManifold::ProjRotation(const Matrix2& R, const Matrix2 & VR) const
            {
                Matrix2 RtVr = R.transpose() * VR;
                Matrix2 Vret = (0.5 * (RtVr - RtVr.transpose()));
                // Matrix2 symmRtVr = 0.5 * (RtVr + RtVr.transpose()); 
                // Matrix2 Vret = VR - R * symmRtVr;
                    return Vret;
            **/
            
            // NO se de donde sale esta
             /*
             Matrix2 RotPriorManifold::ProjRotation(const Matrix2& R, const Matrix2 & VR) const
            {
                Matrix2 RtVr = R.transpose() * VR;
                Matrix2 Vret = (R * 0.5 * (RtVr - RtVr.transpose()));
                    return Vret;
            }
            */
            
            /** Stiefel **/
            Matrix2 RotPriorManifold::ProjRotation(const Matrix2& R, const Matrix2 & VR) const
            {
                Matrix2 RtVr = R.transpose() * VR;
                Matrix2 symmRtVr = 0.5 * (RtVr + RtVr.transpose()); 
                Matrix2 Vret = VR - R * symmRtVr;
                    return Vret;
            }


            Matrix3 RotPriorManifold::retract(const Matrix3 &Y, const Matrix3 &V) const {

              // We use projection-based retraction, as described in "Projection-Like
              // Retractions on Matrix Manifolds" by Absil and Malick
              return project(Y + V); 
        }
        
        

         // We just call here to ProjSphere & ProjRotation2
         Matrix3 RotPriorManifold::Proj(const Matrix3 &Y, const Matrix3 & Vy) const{
         
            Matrix3 V_tan;
            V_tan.setZero();

            Matrix2 R, Vr;
            Vector3 t, Vt;
           
            t = Y.block<3,1>(0,2);
            R = Y.block<2,2>(0,0); 
           
            
            Vt = Vy.block<3,1>(0,2);
            Vr = Vy.block<2,2>(0,0); 
           
                        
            // for the sphere
            V_tan.block<3, 1>(0, 2) = ProjSphere(t, Vt);

            // for the rotation
            V_tan.block<2,2>(0, 0) = ProjRotation(R, Vr);
            
            
            return V_tan;
         }


} // end of essential namespace
