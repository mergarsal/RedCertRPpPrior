#include "RotPrior.h"
#include "RotPriorUtils.h"
#include "RotPriorProblem.h"




// R. optimizayion
#include <experimental/optional>
#include <Optimization/Riemannian/TNT.h>
#include "Optimization/Convex/Concepts.h"
#include "Optimization/Riemannian/GradientDescent.h"
#include "Optimization/Convex/ProximalGradient.h"

#include <functional>
#include <chrono>
#include <math.h>


#include <string>
#include <iomanip>
#include <algorithm>



using namespace std::chrono;
using namespace Optimization;
using namespace Optimization::Convex;




namespace RotPrior{
    EssentialEstimationResult RotPriorClass::getResults(Matrix3 & Rt)
    {

               Matrix9 C;

               // output
               EssentialEstimationResult result;

            
               Matrix3 E_8pts;
               Matrix3 matrix_precon = Matrix3::Identity(); // / p0_.cols();
               
               
                //////////////////////////////////////////////////
                if (options_.estimation_verbose)
                    std::cout << " Constructing data matrix C.\n";
                // Get starting timepoint
                auto start_time_init = high_resolution_clock::now();
                // create data matrix
                C = constructDataMatrix(p0_,p1_,w_);

                auto duration_construction_C = duration_cast<microseconds>(high_resolution_clock::now() - start_time_init);

                auto start_time_bounds = high_resolution_clock::now();
                //////////////////////////////////////////////////


               
               Matrix3 E_init = computeEfromRtPrior(Rt);
               
                               
                
               // Show initial objective
               double f_initial = 0.5 * vec(E_init).transpose() * C * vec(E_init);
               
     

               if (options_.estimation_verbose)
                   std::cout << "Initial objective value f_{init} = " << f_initial << std::endl;
            

	           // Define the problem
	            if (options_.estimation_verbose)         std::cout << "Creating problem " << std::endl;



               /// Define the problem
               RotPriorProblem problem(C);
               problem.setPointCorrespondences(p0_,p1_,w_);
               matrix_precon = problem.computePseudoJacobiPrecon();
               
               // compute precon if needed
               // set precon matrix
                problem.setMatrixPrecon(matrix_precon);

                // Cache these three matrices for gradient & hessian
                ProblemCachedMatrices problem_matrices;

                /// Function handles required by the TNT optimization algorithm
                // Preconditioning operator (optional)
                  std::experimental::optional<Optimization::Riemannian::LinearOperator<Matrix3, Matrix3, ProblemCachedMatrices>> precon;

                  
                    Optimization::Riemannian::LinearOperator<Matrix3, Matrix3, ProblemCachedMatrices> precon_op =
                        [&problem](const Matrix3 &Y, const Matrix3 &Ydot,
                                   const ProblemCachedMatrices & problem_matrices) {
                          return problem.precondition(Y, Ydot);
                        };
                    precon = precon_op;
                  



              // Objective
              Optimization::Objective<Matrix3, double, ProblemCachedMatrices> F =
                  [&problem](const Matrix3 &Y, ProblemCachedMatrices & problem_matrices){
                    return problem.evaluate_objective(Y, problem_matrices);
                  };


             /// Gradient
              Optimization::Riemannian::VectorField<Matrix3, Matrix3, ProblemCachedMatrices> grad_F = 
                                                [&problem](const Matrix3 &Y,
                                                ProblemCachedMatrices & problem_matrices) {
                                                
                // Compute and cache Euclidean gradient at the current iterate
                problem_matrices.NablaF_Y = problem.Euclidean_gradient(Y, problem_matrices);
                // Compute Riemannian gradient from Euclidean one
                return problem.Riemannian_gradient(Y, problem_matrices.NablaF_Y);
              };
              
              

              // Local quadratic model constructor
              Optimization::Riemannian::QuadraticModel<Matrix3, Matrix3, ProblemCachedMatrices> QM =
                  [&problem](
                      const Matrix3 &Y, Matrix3 &grad,
                      Optimization::Riemannian::LinearOperator<Matrix3, Matrix3, ProblemCachedMatrices> &HessOp,
                      ProblemCachedMatrices & problem_matrices) {
                    // Compute and cache Euclidean gradient at the current iterate

                    
                    problem_matrices.NablaF_Y  = problem.Euclidean_gradient(Y, problem_matrices);
                    // Compute Riemannian gradient from Euclidean gradient
                    
                    grad = problem.Riemannian_gradient(Y, problem_matrices.NablaF_Y );

                    
                    // Define linear operator for computing Riemannian Hessian-vector
                    // products
                    HessOp = [&problem](const Matrix3 &Y, const Matrix3 &Ydot,
                                        const ProblemCachedMatrices & problem_matrices) {
                                        Matrix3 Hss = problem.Riemannian_Hessian_vector_product(Y, problem_matrices, Ydot);

                                        return Hss;
                    };
                  };

                  // Riemannian metric

                  // We consider a realization of the product of Stiefel manifolds as an
                  // embedded submanifold of R^{r x dn}; consequently, the induced Riemannian
                  // metric is simply the usual Euclidean inner product
                  Optimization::Riemannian::RiemannianMetric<Matrix3, Matrix3, double, ProblemCachedMatrices>
                      metric = [&problem](const Matrix3 &Y, const Matrix3 &V1, const Matrix3 &V2,
                                          const ProblemCachedMatrices & problem_matrices) {

                        Matrix2 R1, R2;
                        Vector3 t1, t2;
                        R1 = V1.block<2,2>(0,0);
                        R2 = V2.block<2,2>(0,0);

                        t1 = V1.block<3,1>(0,2);
                        t2 = V2.block<3,1>(0,2);

                        return ((R1 * R2.transpose()).trace() + t1.dot(t2));
                      };

              // Retraction operator
              Optimization::Riemannian::Retraction<Matrix3, Matrix3, ProblemCachedMatrices> retraction =
                  [&problem](const Matrix3 &Y, const Matrix3 &Ydot, const ProblemCachedMatrices & problem_matrices) {
                    return problem.retract(Y, Ydot);
                  };





              // Stop timer
              auto stop_init = high_resolution_clock::now();

              auto duration_init = duration_cast<microseconds>( stop_init- start_time_init);





              // set up params for solver
              Optimization::Riemannian::TNTParams<double> params;
              params.gradient_tolerance = options_.tol_grad_norm;
              params.relative_decrease_tolerance = options_.rel_func_decrease_tol;
              params.max_iterations = options_.max_RTR_iterations;
              params.max_TPCG_iterations = options_.max_tCG_iterations;
              params.preconditioned_gradient_tolerance = options_.preconditioned_grad_norm_tol;

              params.stepsize_tolerance = options_.stepsize_tol;
              // params.Delta0 = 1.3;

              params.verbose = options_.verbose;
              /** An optional user-supplied function that can be used to instrument/monitor
               * the performance of the internal Riemannian truncated-Newton trust-region
               * optimization algorithm as it runs. */
              // std::experimental::optional<EssentialTNTUserFunction> user_fcn;

              /// Run optimization!
              auto start_opt = high_resolution_clock::now();

              Optimization::Riemannian::TNTResult<Matrix3, double> TNTResults = Optimization::Riemannian::TNT<Matrix3, Matrix3, double, ProblemCachedMatrices>(
                        F, QM, metric, retraction, Rt, problem_matrices, precon, params);
 

                auto duration_opt = duration_cast<microseconds>(high_resolution_clock::now() - start_opt);
               auto duration_total = duration_cast<microseconds>(high_resolution_clock::now() - start_time_init); 
                

               // Extract the results
               if (options_.estimation_verbose)
               {
                std::cout << "Final point SO(3) x S(2): " << TNTResults.x << std::endl;
                std::cout << "Final objective: " << TNTResults.f << std::endl;
               }

	           Matrix3 Rs_opt = TNTResults.x;
	           Matrix2 R_opt = Rs_opt.block<2,2>(0,0);
	           Vector3 t_opt = Rs_opt.block<3,1>(0,2);

	           Matrix3 E_opt = computeEfromRtPrior(R_opt, t_opt);
	                      

                    

	               //  assign all the results to this struct
	               result.f_hat = TNTResults.f / C.trace();
	               	               
                       result.R_opt = R_opt;
                       result.t_opt = t_opt;
                       result.E_opt = E_opt;
                       // Save time
                       result.elapsed_init_time = duration_init.count();
                       result.elapsed_C_time = duration_construction_C.count();  // in microsecs
                       result.elapsed_iterative_time = duration_opt.count();  // in microsecs
                       result.elapsed_estimation_time = duration_total.count();

           

               return result;

    }; // end of fcn getResult


void RotPriorClass::printResult(EssentialEstimationResult & my_result)
{
    // print data

      std::cout << "Data from estimation\n--------------------\n";

      std::cout << "######## Times [in microseconds]:\n";
      std::cout << "Data Matric C construction: " << my_result.elapsed_C_time << std::endl;
      std::cout << "-----\n Total init: " << my_result.elapsed_init_time << std::endl;

      std::cout << "Iterative Method: " << my_result.elapsed_iterative_time << std::endl;

      std::cout << "---------------------\n";
      std::cout << "Total time: " << my_result.elapsed_estimation_time << std::endl << std::endl;

      std::cout << "\n Recovered R:\n" << my_result.R_opt << std::endl;
      std::cout << "\n Recovered t:\n" << my_result.t_opt << std::endl;

    

};  //end of print fcn

}  // end of essential namespace
