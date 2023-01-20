#include <Eigen/Core>
#include <iostream>
#include <LBFGS.h>
#include <cmath>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace LBFGSpp;

double foo(const VectorXd& x, VectorXd& grad)
{
    
    std::cout << "input: " <<grad[1] << std::endl;
    const int n = x.size();
    VectorXd d(n);
    for(int i = 0; i < n; i++)
        d[i] = i;

    double f = (x - d).squaredNorm();
    grad.noalias() = 2.0 * (x - d);
    std::cout << "output: " <<grad[1] << std::endl;
    return f;
}

double lorentzianE(const VectorXd& x)
{
    double E = 0;
    for(int i=0; i<x.size();i++)
    {
        for(int j=0; j<x.size();j++)
        {
            if(i<j)
            {
                double denorm = pow((x[i] - x[j] -3),2)+5;
                E += -1/denorm*1000;
            }
        }
    }
    std::cout << "E: " << E << std::endl;
    return E;
}

VectorXd lorentzianGrad(const VectorXd& x)
{
    VectorXd grad(x.size());
    const int size = x.size();
    for(int i=0;i<size;i++)
    {
        double sumDer =0;
        for(int j=0;j<size;j++)
        {
            if(i<j)
            {
                double denorm = pow((x[i] - x[j] -3),2)+5;
                sumDer += 2*(x[i]-x[j]-3)/pow(denorm, 2);
               
            }
            else if(i>j)
            {
                double denorm = pow((x[j] - x[i] -3),2)+5;
                sumDer += -2*(x[j]-x[i]-3)/pow(denorm, 2);
                
            }
        };
        grad[i] = sumDer*1000;   
    }
    std::cout << grad.transpose() << std::endl;
    return grad;
}



double lorentzian2D(const VectorXd& x, VectorXd& grad)
{    
    double E = lorentzianE(x);
    grad = lorentzianGrad(x);
    return E;
};




int main()
{
    const int n = 2;
    LBFGSParam<double> param;
    
    param.max_iterations = 3;
    param.ftol = 1e-7;
    param.epsilon = 1e-6;
    std::cout << param.epsilon << std::endl;
    LBFGSSolver<double> solver(param);

    VectorXd x(n);

    x[0] = 10;
    x[1] = 40;
    std::cout << "initial x = \n" << x.transpose() << std::endl;

    


    
    double fx;
    int niter = solver.minimize(lorentzian2D, x, fx);

    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = " << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;

    return 0;
}
