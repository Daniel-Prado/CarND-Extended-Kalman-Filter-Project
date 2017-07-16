#include "kalman_filter.h"


using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}


void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  
  // <DPR>
  x_ = F_ * x_; // as seen in section 8 of Lesson 5.
  MatrixXd Ft = F_.transpose(); // as seen in section 9 of Lesson 5.
  P_ = F_ * P_ * Ft + Q_;
  // </DPR>

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  // <DPR>  As seen in sections 7 and 13 of lesson 5.
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
  
  // </DPR>

}

// <DPR>
const float KalmanFilter::rewrapAngleRestricted(const float angle) {

  const float TWO_PI = 2 * M_PI;
  
  if(angle > M_PI)
    return angle - TWO_PI;
  else if(angle < -M_PI)
    return angle + TWO_PI;
  
  return angle;
}
// </DPR>


void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  
  // <DPR>
  float rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  float phi = atan2(x_(1), x_(0));
  float rho_dot;
  if (fabs(rho) < 0.0001) {
    rho_dot = 0;
  } else {
    rho_dot = (x_(0)*x_(2) + x_(1)*x_(3))/rho;
  }
  VectorXd z_pred(3);
  z_pred << rho, phi, rho_dot;
  VectorXd y = z - z_pred;
  // As required in the project, we should ensure that the resulting phi angle
  // of the y vector calculated above is between -pi and pi. This is not guaranteed yet because
  // for example, z could be close to pi, and z_pred close to -pi and hence the substraction
  // would be close to 2*pi.
  y(1) = rewrapAngleRestricted(y(1));
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
  
  // </DPR>
  
}
