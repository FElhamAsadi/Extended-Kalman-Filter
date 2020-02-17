#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0.0, 0.0, 0.0, 0.0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()
     || estimations.size() == 0) {
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }

  // accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];

    // coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse = rmse/estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;
}


MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  
  // Calculate a Jacobian here.
  float mag_x = sqrt(pow(x_state(0),2.0) + pow(x_state(1),2.0));
  float cross_p = x_state(2)*x_state(1) - x_state(3)*x_state(0);
     
  MatrixXd Jacob = MatrixXd(3,4);

  if (fabs(mag_x) < 0.00001){
    cout << "CalculateJacobian () - Error - Division by Zero" << endl;
    return Jacob;
    
  } 
  
  Jacob << x_state(0)/mag_x, x_state(1)/mag_x, 0.0 , 0.0,
         -x_state(1)/pow(mag_x,2.0), x_state(0)/pow(mag_x,2.0), 0.0, 0.0,
         x_state(1)*cross_p/pow(mag_x,3.0),  -x_state(0)*cross_p/pow(mag_x,3.0), x_state(0)/mag_x, x_state(1)/mag_x;
  
  return Jacob;
  
}

VectorXd Tools::CartesianToPolar(const VectorXd& x_state) {
    
    VectorXd Polar = VectorXd(3);
  
    float mag_x = sqrt(pow(x_state(0), 2.0) + pow(x_state(1), 2.0));
	if (mag_x < 0.0001) {
      cout << "CartesianToPolar () - Error - Division by Zero" << endl;
      return Polar;

	}
  
	else {
      // normalize the angle between -pi to pi
      float Phi = atan2(x_state(1), x_state(0));
      while (Phi > M_PI){
      	Phi  -= 2 * M_PI;
      }
  	  while (Phi < -M_PI){
    	Phi  += 2 * M_PI;
      }
      
      Polar << mag_x,
			  Phi,
			  (x_state(0) * x_state(2) + x_state(1) * x_state(3)) / mag_x;
	}
     return Polar;
}