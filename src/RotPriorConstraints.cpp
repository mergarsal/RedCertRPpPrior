#include "RotPriorConstraints.h"


namespace RotPrior
{

    void createLeftEConstraints(std::vector<Matrix9> & As)
    {
        
        size_t e1=0, e4=1, e2=2, e8=3, e3=4, e6=5; 
        size_t t1=6, t2=7, t3=8; 
        
        /** Left formulation **/
        Matrix9 Ai; 
        Ai.setZero(); 
        Ai(t1,t1)=1;
        Ai(t2,t2)=1;
        Ai(t3,t3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e4,e4)=1;
        Ai(e6,e6)=1;
        Ai(t1,t1)=-1;
        Ai(t3,t3)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e3,e3)=1;
        Ai(e8,e8)=1;
        Ai(e1,e1)=1;
        Ai(t1,t1)=-1;
        Ai(t2,t2)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e1,e4)=1;
        Ai(e3,e6)=1;
        Ai(t1,t2)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e2,e8)=1;
        Ai(t1,t3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e3,e4)=-1;
        Ai(e6,e1)=1;
        Ai(t2,t3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
       
        Ai.setZero();
        Ai(e1,e1)=1;
        Ai(e2,e2)=1;
        Ai(e3,e3)=1;
        Ai(t3,t3)=-1;
        Ai(t2,t2)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
    
    }


    void createRightEConstraints(std::vector<Matrix9> & As)
    {
        
        size_t e1=0, e4=1, e2=2, e8=3, e3=4, e6=5; 
        size_t q1=6, q2=7, q3=8; 
        
        /** Right formulation **/
        Matrix9 Ai; 
        Ai.setZero();
        Ai(q1,q1)=1;
        Ai(q2,q2)=1;
        Ai(q3,q3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        Ai.setZero();
        Ai(e3,e3)=1;
        Ai(e6,e6)=1;
        Ai(e1,e1)=1;
        Ai(q1,q1)=-1;
        Ai(q2,q2)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        Ai.setZero();
        Ai(e2,e2)=1;
        Ai(e8,e8)=1;
        Ai(q1,q1)=-1;
        Ai(q3,q3)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        Ai.setZero();
        Ai(e1,e1)=1;
        Ai(e4,e4)=1;
        Ai(e3,e3)=1;
        Ai(q2,q2)=-1;
        Ai(q3,q3)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        Ai.setZero();
        Ai(e1,e2)=1;
        Ai(e3,e8)=-1;
        Ai(q1,q2)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        Ai.setZero();
        Ai(e4,e6)=1;
        Ai(q1,q3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        Ai.setZero();
        Ai(e2,e3)=1;
        Ai(e8,e1)=1;
        Ai(q2,q3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
    }

    
    void createBothEConstraints(std::vector<Matrix12> & As)
    {
        
        size_t e1=0, e4=1, e2=2, e8=3, e3=4, e6=5; 
        size_t t1=6, t2=7, t3=8, q1=9, q2=10, q3=11; 
        
        /** Both formulation **/
        Matrix12 Ai; 
        Ai.setZero(); 
        Ai(t1,t1)=1;
        Ai(t2,t2)=1;
        Ai(t3,t3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e4,e4)=1;
        Ai(e6,e6)=1;
        Ai(t1,t1)=-1;
        Ai(t3,t3)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e3,e3)=1;
        Ai(e8,e8)=1;
        Ai(e1,e1)=1;
        Ai(t1,t1)=-1;
        Ai(t2,t2)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e1,e4)=1;
        Ai(e3,e6)=1;
        Ai(t1,t2)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e2,e8)=1;
        Ai(t1,t3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e3,e4)=-1;
        Ai(e6,e1)=1;
        Ai(t2,t3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
       
        Ai.setZero();
        Ai(e1,e1)=1;
        Ai(e2,e2)=1;
        Ai(e3,e3)=1;
        Ai(t3,t3)=-1;
        Ai(t2,t2)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        

        Ai.setZero();
        Ai(q1,q1)=1;
        Ai(q2,q2)=1;
        Ai(q3,q3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        Ai.setZero();
        Ai(e3,e3)=1;
        Ai(e6,e6)=1;
        Ai(e1,e1)=1;
        Ai(q1,q1)=-1;
        Ai(q2,q2)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        Ai.setZero();
        Ai(e2,e2)=1;
        Ai(e8,e8)=1;
        Ai(q1,q1)=-1;
        Ai(q3,q3)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        /*
        Ai.setZero();
        Ai(e1,e1)=1;
        Ai(e4,e4)=1;
        Ai(e3,e3)=1;
        Ai(q2,q2)=-1;
        Ai(q3,q3)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        */
        
        Ai.setZero();
        Ai(e1,e2)=1;
        Ai(e3,e8)=-1;
        Ai(q1,q2)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        Ai.setZero();
        Ai(e4,e6)=1;
        Ai(q1,q3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        Ai.setZero();
        Ai(e2,e3)=1;
        Ai(e8,e1)=1;
        Ai(q2,q3);
        As.push_back(0.5 * (Ai + Ai.transpose()));
    
    }
    
    void createAdjugateEConstraints(std::vector<Matrix12> & As)
    {
        
        size_t e1=0, e4=1, e2=2, e8=3, e3=4, e6=5; 
        size_t t1=6, t2=7, t3=8, q1=9, q2=10, q3=11; 
        
        /** Left formulation **/
        Matrix12 Ai; 
        Ai.setZero(); 
        Ai(t1,t1)=1;
        Ai(t2,t2)=1;
        Ai(t3,t3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e4,e4)=1;
        Ai(e6,e6)=1;
        Ai(t1,t1)=-1;
        Ai(t3,t3)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e3,e3)=1;
        Ai(e8,e8)=1;
        Ai(e1,e1)=1;
        Ai(t1,t1)=-1;
        Ai(t2,t2)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e1,e4)=1;
        Ai(e3,e6)=1;
        Ai(t1,t2)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e2,e8)=1;
        Ai(t1,t3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e3,e4)=-1;
        Ai(e6,e1)=1;
        Ai(t2,t3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
       
        /*
        Ai.setZero();
        Ai(e1,e1)=1;
        Ai(e2,e2)=1;
        Ai(e3,e3)=1;
        Ai(t3,t3)=-1;
        Ai(t2,t2)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        */
        
        /** Right formulation **/
        

        Ai.setZero();
        Ai(q1,q1)=1;
        Ai(q2,q2)=1;
        Ai(q3,q3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        /*
        Ai.setZero();
        Ai(e3,e3)=1;
        Ai(e6,e6)=1;
        Ai(e1,e1)=1;
        Ai(q1,q1)=-1;
        Ai(q2,q2)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        */
        
        Ai.setZero();
        Ai(e2,e2)=1;
        Ai(e8,e8)=1;
        Ai(q1,q1)=-1;
        Ai(q3,q3)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        Ai.setZero();
        Ai(e1,e1)=1;
        Ai(e4,e4)=1;
        Ai(e3,e3)=1;
        Ai(q2,q2)=-1;
        Ai(q3,q3)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        Ai.setZero();
        Ai(e1,e2)=1;
        Ai(e3,e8)=-1;
        Ai(q1,q2)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        Ai.setZero();
        Ai(e4,e6)=1;
        Ai(q1,q3)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));

        Ai.setZero();
        Ai(e2,e3)=1;
        Ai(e8,e1)=1;
        Ai(q2,q3);
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        /** Norm E **/
        Ai.setZero();
        Ai(e1,e1)=2;
        Ai(e2,e2)=1;
        Ai(e3,e3)=2;
        Ai(e4,e4)=1;
        Ai(e6,e6)=1;
        Ai(e8,e8)=1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        /** Adjugate of E **/
        Ai.setZero();
        Ai(e6,e8)=-1;
        Ai(t1,q1)=-1;        
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e3,e8)=1;
        Ai(e1,e2)=-1;
        Ai(t2,q1)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e2,e6)=1;
        Ai(t3,q1)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e1,e4)=-1;
        Ai(e6,e3)=-1;
        Ai(t1,q2)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e3,e3)=1;
        Ai(e1,e1)=1;
        Ai(t2,q2)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e3,e4)=1;
        Ai(e1,e6)=-1;
        Ai(t3,q2)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e4,e8)=1;
        Ai(t1,q3)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e2,e3)=-1;
        Ai(e1,e8)=-1;
        Ai(t2,q3)=-1;
        As.push_back(0.5 * (Ai + Ai.transpose()));
        
        Ai.setZero();
        Ai(e2,e4)=-1;
        Ai(t3,q3)=-1,
        As.push_back(0.5 * (Ai + Ai.transpose()));
    }


}  // end of namespace
