#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225,  0.0,
               0.0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09,  0.0,  0.0,
              0.0, 0.0009,  0.0,
              0.0,  0.0, 0.09;

  H_laser_ << 1.0, 0.0, 0.0, 0.0,
              0.0, 1.0, 0.0, 0.0;
  
  // the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1.0,  0.0, 0.0,  0.0,
             0.0,  1.0,  0.0, 0.0,
             0.0,  0.0, 1.0,  0.0,
             0.0,  0.0,  0.0, 1.0;

  // set the acceleration noise components
  noise_ax = 9.0;
  noise_ay = 9.0;

  //ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ = MatrixXd::Zero(4, 4);
  
 
   // state covariance matrix P
   ekf_.P_ = MatrixXd(4, 4);
   ekf_.P_ << 1.0, 0.0, 0.0, 0.0,
              0.0,  1.0, 0.0,  0.0,
              0.0,  0.0, 10.0,  0.0,
              0.0,  0.0,  0.0, 10.0;
  
   // Jacobian matrix
   ekf_.Hj_ = MatrixXd(3,4);
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() { }

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * Initialize the state ekf_.x_ with the first measurement.
     * Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    cout << "EKF Initialization: " <<endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1,1,5.2,0;
    
    // first measurement  
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates  and initialize state.
      ekf_.R_ = R_radar_; 
      ekf_.x_(0) =  measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]);
      ekf_.x_(1) =  measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]);

      ekf_.H_ = VectorXd(3);
	  ekf_.H_ << tools.CartesianToPolar(ekf_.x_);
      ekf_.Hj_<< tools.CalculateJacobian(ekf_.x_);
    }
    
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

      // Initialize state.
      ekf_.R_ = R_laser_; 
      ekf_.H_ = H_laser_;
      ekf_.x_(0) = measurement_pack.raw_measurements_[0];
      ekf_.x_(1) = measurement_pack.raw_measurements_[1];
      	
    }
    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
  
    return;
  }

  /**
   * Prediction: 
   * Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_ (0,2) = dt;
  ekf_.F_ (1,3) = dt;
  
  ekf_.Q_ << pow(dt,4.0)/4.0*noise_ax,  0.0, pow(dt,3.0)/2.0*noise_ax,  0.0,
              0.0, pow(dt,4.0)/4.0*noise_ay,  0.0, pow(dt,3.0)/2.0*noise_ay,
             pow(dt,3.0)/2.0*noise_ax,  0.0, pow(dt,2.0)*noise_ax,  0.0,
              0.0, pow(dt,3.0)/2.0*noise_ay,  0.0, pow(dt,2.0)*noise_ay; 

            
  ekf_.Predict();

  
  /**
   *  Update :
   *  Use the sensor type to perform the update step.
   *  Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      
    // Radar updates
    ekf_.R_ = R_radar_;
    ekf_.H_ = VectorXd(3);
	ekf_.H_ << tools.CartesianToPolar(ekf_.x_);
    ekf_.Hj_<< tools.CalculateJacobian(ekf_.x_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);


  } else {
    
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "------------------------------- " << endl;
  cout << "x_ = " <<endl;
  cout << ekf_.x_ << endl;
  cout << "   " <<endl;
  cout << "P_ = " << endl;
  cout << ekf_.P_ << endl;
}
