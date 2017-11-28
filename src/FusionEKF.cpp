#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
//#include <fstream>


#define EPS 0.0001 // A very small number

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
    noise_ax = 9;
    noise_ay = 9;
    

    
    ofs.open("../data/output.txt", std::ofstream::out | std::ofstream::app);
    
    

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
    
    
    //  Initialize the transition matrix F_
    ekf_.F_ = MatrixXd(4,4);
    ekf_.Q_ = MatrixXd(4,4);
    ekf_.F_ << 1, 0, 1, 0,
			  0, 1, 0, 1,
			  0, 0, 1, 0,
			  0, 0, 0, 1;
    
    
    // Initialize Measurement matrix for Laser measurement update
    
    H_laser_ << 1,0,0,0,
                0,1,0,0;
    
    // Initialize state covariance matrix  P
    ekf_.P_ = MatrixXd(4,4);
    ekf_.P_ <<  1,0,0,0,
                0,1,0,0,
                0,0,1000,0,
                0,0,0,1000;
    


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {ofs.close();}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack, const VectorXd &gt_values) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
        float rho = measurement_pack.raw_measurements_[0];
        float theta = measurement_pack.raw_measurements_[1];
        float rho_dot = measurement_pack.raw_measurements_[2];
        float px = rho * cos(theta);
        float py = rho * sin(theta);
        float vx = rho_dot * cos(theta);
        float vy = rho_dot * sin(theta);
        
        ekf_.x_ << px,py,vx,vy;
        //std::cout << " RADAR DATA initialized:";
      
        
        
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
        ekf_.x_ << measurement_pack.raw_measurements_[0],measurement_pack.raw_measurements_[1],0,0;
    }

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
      // Deal with the special case initialisation problems
      if (fabs(ekf_.x_(0)) < EPS and fabs(ekf_.x_(1)) < EPS){
          ekf_.x_(0) = EPS;
          ekf_.x_(1) = EPS;
      }
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

    
    //compute the elapsed time between the current and previous measurements
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;
    float dt_2  = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;
    if(dt_4 > MAXFLOAT){dt_4 = MAXFLOAT;}
    //float noise_ax = 9;
    //float noise_ay = 9;
    //std::cout << "dt_:"<< dt << std::endl;
    //std::cout << "dt_2:"<< dt_2 << std::endl;
    //std::cout << "dt_3:"<< dt_3 << std::endl;
    //std::cout << "dt_4:"<< dt_4 << std::endl;
    
    //std::cout << "dt_4/4*noise_ax:" << dt_4/(4*noise_ax);
    //std::cout << "dt_3/2*noise_ax:" << dt_4/(4*noise_ax);
    
    float maxxQ;
    float maxyQ;
    
    if (((dt_4/4)*noise_ax) > MAXFLOAT){maxxQ = MAXFLOAT;}
    else { maxxQ = (dt_4/4)*noise_ax ;}
    
    if (((dt_4/4)*noise_ay) > MAXFLOAT){maxyQ = MAXFLOAT;}
    else { maxyQ = dt_4/4*noise_ay ;}
    
    
    
    // Modify the F Matrix (State Transition Matrix)  so that the time is integrated
    ekf_.F_(0,2) = dt;
    ekf_.F_(1,3) = dt;
    
    //set the process covariance Matrix Q
    /*
    ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
                0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
                dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
                0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
     */
    
    ekf_.Q_ <<  maxxQ, 0, (dt_3/2)*noise_ax, 0,
                0, maxyQ, 0, (dt_3/2)*noise_ay,
                (dt_3/2)*noise_ax, 0, dt_2*noise_ax, 0,
                0, (dt_3/2)*noise_ay, 0, dt_2*noise_ay;
    
    //std::cout<< "calling prediction:";
  ekf_.Predict();
    //std::cout<< "Done prediction:";


  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
      ekf_.R_ = Eigen::MatrixXd(3,3);
      ekf_.R_ = R_radar_;
      
      Hj_ = tools.CalculateJacobian(ekf_.x_);
      
      ekf_.Hj_ = Hj_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
      
  } else {
    // Laser updates
      ekf_.R_ = Eigen::MatrixXd(2,2);
      ekf_.R_ = R_laser_;
      ekf_.H_ = H_laser_;
      ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
    
    
    // output the estimation
    ofs << ekf_.x_[0] << "\t";
    ofs << ekf_.x_(1) << "\t";
    ofs << ekf_.x_(2) << "\t";
    ofs << ekf_.x_(3) << "\t";
    
    
    
    // output the measurements
    if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        // output the estimation
        ofs << measurement_pack.raw_measurements_[0] << "\t";
        ofs << measurement_pack.raw_measurements_[1] << "\t";
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // output the estimation in the cartesian coordinates
        float ro = measurement_pack.raw_measurements_[0];
        float phi = measurement_pack.raw_measurements_[1];
        ofs << ro * cos(phi) << "\t"; // p1_meas
        ofs << ro * sin(phi) << "\t"; // ps_meas
    }
    
    
    
    // output the ground truth packages
    ofs << gt_values(0) << "\t";
    ofs << gt_values(1) << "\t";
    ofs << gt_values(2) << "\t";
    ofs << gt_values(3) << "\n";
    //ofs << "\n";
    //cout << "ground_truth:" << ground_truth[1] << "\n";
    
   }
