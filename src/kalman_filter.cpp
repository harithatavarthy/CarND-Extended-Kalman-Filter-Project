#include "kalman_filter.h"
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
    
    x_= F_ * x_ ;
    MatrixXd Ft_ = F_.transpose();
    P_ = F_ * P_ * Ft_ + Q_;
    
    /*
    std::cout << "Inside Predict function";
    std::cout << "F_:"<< F_ << std::endl;
    std::cout << "Ft_:"<< Ft_ << std::endl;
    std::cout << "Q_:"<< Q_ << std::endl;
    
    std::cout << "x_:"<< x_ << std::endl;
    std::cout << "P_:"<< P_ << std::endl;
     */
    
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
    
    //std::cout << "Inside KalmanFilter update function";
    //std::cout << "H_:"<< H_ << std::endl;
    //std::cout << "Predicted x_:"<< x_ << std::endl;
    //std::cout << "Predicted P_:"<< P_ << std::endl;
    VectorXd z_pred = H_ * x_;
    //std::cout << "Z_PRED:"<<z_pred<< std::endl;
    
    VectorXd y_ = z - z_pred;
    //std::cout << "y_:" << y_ << std::endl;
    MatrixXd Ht_ = H_.transpose();
    //std::cout << "Ht_:"<< Ht_ << std::endl;
    MatrixXd S_ = H_ * P_ * Ht_ + R_;
    //std::cout << "S_:"<< S_ << std::endl;
    MatrixXd Si_ = S_.inverse();
    //std::cout << "Si_:"<< Si_ << std::endl;
    MatrixXd K_ =  P_ * Ht_ * Si_;
    //std::cout << "K_:"<< K_ << std::endl;
    
    //new state
    x_ = x_ + (K_ * y_);
    //std::cout << "New x_:"<< x_ << std::endl;
    long x_size = x_.size();
    //std::cout << "size:"  << x_size;
    MatrixXd I_ = MatrixXd::Identity(x_size, x_size);
    P_ = (I_ - K_ * H_) * P_;
    //std::cout << "New P_:"<< P_ << std::endl;
     
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
    // Convert the predicted measurement from cartesian to polar form
    
    VectorXd z_pred_radar =  Eigen::VectorXd(3);
    float px = x_[0];
    float py = x_[1];
    float vx = x_[2];
    float vy = x_[3];
    
    float c1 = px*px+py*py;
    float c2 = sqrt(c1);
    z_pred_radar[0] = c2;
    
    if(px!= 0)
    {
    //z_pred_radar[1] = atan(py/px);
        z_pred_radar[1] = atan2(py,px);
    }
    else{ z_pred_radar[1] = 0;}
    //std::cout << "Theta predicted : " << z_pred_radar[1];
    
    z_pred_radar[2] = (px*vx + py*vy)/c2;
    
    
    VectorXd y_radar = z - z_pred_radar;
    //std::cout << "Theta error : " << y_radar[1];
    
    if(y_radar[1] > M_PI){y_radar[1] = y_radar[1] - 2*M_PI;}
    else if (y_radar[1] < (-1*(M_PI))){y_radar[1] = y_radar[1] + 2*M_PI;}

    //std::cout << "Theta error after adjustment : " << y_radar[1];
    
    MatrixXd Hjt_ = Hj_.transpose();
    MatrixXd S_ = Hj_ * P_ * Hjt_ + R_;
    
    MatrixXd Si_ = S_.inverse();
    MatrixXd K_ =  P_ * Hjt_ * Si_;
    
    //new state
    x_ = x_ + (K_ * y_radar);
    long x_size = x_.size();
    //std::cout << "size:"  << x_size;
    MatrixXd I_ = MatrixXd::Identity(x_size, x_size);
    P_ = (I_ - K_ * Hj_) * P_;
    
    
}
