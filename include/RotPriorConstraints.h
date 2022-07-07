#pragma once 


#include <Eigen/Core>

#include "RotPriorTypes.h"



namespace RotPrior{
    
        void createLeftEConstraints(std::vector<Matrix9> & As);
        
        void createRightEConstraints(std::vector<Matrix9> & As);

        void createBothEConstraints(std::vector<Matrix12> & As);
        
        void createAdjugateEConstraints(std::vector<Matrix12> & As);


} // end of essential namespace
