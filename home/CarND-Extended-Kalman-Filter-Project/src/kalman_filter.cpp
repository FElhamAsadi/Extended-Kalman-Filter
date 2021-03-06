#include "kalman_filter.h"
using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in, MatrixXd &Hj_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
  Hj_ = Hj_in;
}

void KalmanFilter::Predict() {
  /**
   * predict the state
   */
  x_ = F_ * x_;
  P_ = (F_ * P_) * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  
  // update the state by using Kalman Filter equations
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = (P_ * H_.transpose() ) * S.inverse();

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  
  //update the state by using Extended Kalman Filter equations 
  VectorXd y = z - H_;
  
  // normalize the angle between -pi to pi
  while(y(1) > M_PI){
    y(1) -= 2 * M_PI;
  }
  while(y(1) < -M_PI){
    y(1) += 2 * M_PI;
  }
  
  MatrixXd S = (Hj_ * P_) * Hj_.transpose() + R_;
  MatrixXd K = (P_ * Hj_.transpose() ) * S.inverse();

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj_) * P_;

}
