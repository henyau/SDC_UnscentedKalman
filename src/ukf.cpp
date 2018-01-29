#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

//Helper func to return normalized (wrapped -pi to pi) angle in radians
double AngleWrap(double inAng)
{
	double inAng_wrapped = fmod(inAng+M_PI, 2*M_PI);
	if (inAng_wrapped<0)
		inAng_wrapped += 2*M_PI;
	return (inAng_wrapped - M_PI);
}

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 5.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.8;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  
  n_x_ = 5; // state vector length
  n_aug_ = n_x_ + 2; // augmented state vector length, includes noise vector 
  const int  n_sig = n_aug_*2+1; //number of sigma points to evaluate at

  lambda_ = 3-n_aug_; //design parameter

  weights_ = VectorXd(n_sig);//set weights
  weights_(0) = lambda_/(lambda_+n_aug_);

  for(int i = 1; i<n_sig; i++)
  	weights_(i) = 0.5f/(lambda_+n_aug_);

  is_initialized_ = false;

  Xsig_pred_ = MatrixXd::Zero(n_x_, n_sig);
  time_us_ = 0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  process incoming measurements depending on sensor type
  */
  //if not yet initialized, use first measurement to set ICs
  if(!is_initialized_){
	if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
		VectorXd xvec = VectorXd(meas_package.raw_measurements_);
		x_<< xvec[0] * cos(xvec[1]), xvec[0] * sin(xvec[1]), 0,0,0;	
  	}
	else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
		VectorXd xvec = VectorXd(meas_package.raw_measurements_);
		x_<< xvec[0], xvec[1],0,0,0;		
	}
	//TODO, for now use scaled I, but could use knowledge of noise distribution for more accurate start
	P_ = 0.05f*MatrixXd::Identity(5,5);	
  	is_initialized_ = true;
	time_us_ = meas_package.timestamp_;	
  }
  else{
	double delta_t = (meas_package.timestamp_-time_us_)*1E-6;
        if(delta_t > 1.0f)//if reset, the start time is wrong
	{
		delta_t = 0.05f;
		is_initialized_ = false;
	}
	time_us_ = meas_package.timestamp_;
//        cout<<"delta_t: "<< delta_t<<endl;
   	Prediction(delta_t);
	if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
		UpdateRadar(meas_package);		
	else if(meas_package.sensor_type_ == MeasurementPackage::LASER)
		UpdateLidar(meas_package);
  }


}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
   Predict sigma points, the state, and the state covariance matrix.
  */
  
  const int  n_sig = n_aug_*2+1; //number of sigma points to evaluate at
 
  
 //1.) generate sigma points
  //create matrix of augmented state vectors evaluated at sigma points  
  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, n_sig);
  VectorXd x_aug = VectorXd::Zero(n_aug_);//augmented state vector
 
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_+1) = 0;

  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug.block(0,0,n_x_,n_x_) = P_;
 
  //Q block
  P_aug(n_x_,n_x_) = std_a_*std_a_;
  P_aug(n_x_+1,n_x_+1) = std_yawdd_*std_yawdd_;

  MatrixXd P_sqrt = P_aug.llt().matrixL();


  for(int i = 0; i< n_sig; i++)
	Xsig_aug.col(i) = x_aug;
  

  Xsig_aug.block(0, 1, n_aug_, n_aug_) += (sqrt(lambda_+n_aug_))*P_sqrt;
  Xsig_aug.block(0, n_aug_+1, n_aug_, n_aug_) -= (sqrt(lambda_+n_aug_))*P_sqrt;
  
  //2.) Predict state at sigma points
  for (int i = 0; i<n_sig; i++){
    //Xsig_pred_.col(i)(3) = AngleWrap(Xsig_pred_.col(i)(3));

      //use labels for clarity
    float vel = Xsig_aug(2, i);
    float yaw = AngleWrap(Xsig_aug(3, i));
    float yawd = Xsig_aug(4, i);
    float nu_a = Xsig_aug(5, i);
    float nu_yawdd = Xsig_aug(6, i);
    
    VectorXd uncert_v = VectorXd(n_x_);
    VectorXd update_v = VectorXd(n_x_);
    uncert_v<<0.5*delta_t*delta_t*cos(yaw)*nu_a,0.5*delta_t*delta_t*sin(yaw)*nu_a, delta_t*nu_a,0.5*delta_t*delta_t*nu_yawdd,delta_t*nu_yawdd;
   

    if(fabs(yawd) <=1E-4)
    	update_v <<vel*cos(yaw)*delta_t, vel*sin(yaw)*delta_t,0,yawd*delta_t,0;
    else
    	update_v <<(vel/yawd)*(sin(yaw+yawd*delta_t)-sin(yaw)), (vel/yawd)*(-cos(yaw+yawd*delta_t)+cos(yaw)),0,yawd*delta_t,0;
    

    Xsig_pred_.col(i) = Xsig_aug.col(i).head(n_x_)+update_v+uncert_v;
  }
 
  //3.) compute mean x and covariance 
  //predict state mean
  x_.fill(0.0f);//clear
  P_.fill(0.0f);//clear
  
  for (int i=0; i<n_sig; i++){
      x_ += weights_(i)*Xsig_pred_.col(i);
  }
  //predict state covariance matrix
  x_(3) = AngleWrap(x_(3)); 

  VectorXd tempV;
  for (int i=0; i<n_sig; i++){
  	tempV = Xsig_pred_.col(i)-x_;
        tempV(3) = AngleWrap(tempV(3));
        P_ += weights_(i)*tempV*tempV.transpose();
  }
//  cout<<"P: "<<P_<<endl;
//  cout<<"x: "<<x_<<endl;
  
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Modify the state vector, x_, and covariance, P_.

  Calculate the lidar NIS.
  */

  VectorXd zvec = VectorXd(meas_package.raw_measurements_);
//  cout<<"start Update Lidar"<<endl;
//  cout<<"x_start: "<< x_ <<endl;
  int n_z = 2; //lidar measurement dimension
  int n_sig = 2*n_aug_+1;
  MatrixXd Zsig = MatrixXd::Zero(n_z, n_sig);

   //mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(n_z);  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z,n_z);
 
 //create matrix of predicted measurements (at sigma points)
  Zsig = Xsig_pred_.block(0,0,2, n_sig);//for lidar data, no conversion in pos needed

  //calculate mean predicted measurement
  for(int i = 0; i<n_sig; i++){
      z_pred += Zsig.col(i)*weights_(i);
  }
 
  VectorXd zdiff; 
  //calculate innovation covariance matrix S
  for(int i = 0; i<n_sig; i++){
	zdiff = Zsig.col(i) - z_pred;
	zdiff(1) = AngleWrap(zdiff(1));
	S += weights_(i)*(zdiff*zdiff.transpose());
  }
  S(0, 0) +=std_laspx_*std_laspx_;
  S(1, 1) +=std_laspy_*std_laspy_;
 
 //Compute cross correlation matrix T 
  MatrixXd T = MatrixXd::Zero(n_x_,n_z);
  VectorXd xdiff;
  for(int i = 0; i< n_sig; i++){
	xdiff = Xsig_pred_.col(i) - x_;
	zdiff = Zsig.col(i) - z_pred;
	//normalize angles
	zdiff(1) = AngleWrap(zdiff(1));
	xdiff(3) = AngleWrap(xdiff(3));

  	T += weights_(i)*xdiff*zdiff.transpose();
  }

  //Kalman gain matrix
  MatrixXd K = T*S.inverse();
  //update dates and covariance matrix!
  zdiff = zvec - z_pred;
	//normalize angles
  zdiff(1) = AngleWrap(zdiff(1));

  x_ += K*zdiff;
  P_ -= K*S*K.transpose();

  x_(3) = AngleWrap(x_(3));
//  cout<<"P: "<<P_<<endl;
//  cout<<"x: "<<x_<<endl;
//  cout<<"K "<< K <<endl;
//  cout<<"Xsig_pred: "<<Xsig_pred_<<endl;

//  cout<<"Updated Lidar: "<<endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  
  VectorXd zvec = VectorXd(meas_package.raw_measurements_);
  //no available info on orientation or velocities
  //x_<< xvec[0] * cos(xvec[1]), xvec[0] * sin(xvec[1]), 0,0,0;		
  cout<<"start Update Radar"<<endl;
  cout<<"x_start: "<< x_ <<endl;
  int n_z = 3; //radar measurement dimension
  int n_sig = 2*n_aug_+1;
  MatrixXd Zsig = MatrixXd::Zero(n_z, n_sig);
 // Zsig.fill(0.0f);

   //mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(n_z);  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z,n_z);
//  S.fill(0.0f);

  double px;
  double py;
  double v;
  double sai;
//  float said;
 
  double rho;
  double phi;
  double phid;
  
//create matrix of predicted measurements (at sigma points)
  for(int i = 0; i<n_sig; i++){
      px = Xsig_pred_.col(i)(0);
      py = Xsig_pred_.col(i)(1);
      v = Xsig_pred_.col(i)(2);
      sai = Xsig_pred_.col(i)(3);
//      said = Xsig_pred_.col(i)(4);
      
      rho = sqrt(px*px+py*py);
      phi = atan2(py,px);
      phid = (px*v*cos(sai)+py*v*sin(sai))/rho;
      Zsig.col(i)(0) = rho;
      Zsig.col(i)(1) = phi;
      Zsig.col(i)(2) = phid;
  }


  //calculate mean predicted measurement
  //z_pred.fill(0.0f);
  for(int i = 0; i<n_sig; i++){
      z_pred += Zsig.col(i)*weights_(i);
  }
 

  VectorXd zdiff; 
  //calculate innovation covariance matrix S
  for(int i = 0; i<n_sig; i++){
	zdiff = Zsig.col(i) - z_pred;
	zdiff(1) = AngleWrap(zdiff(1));
	S += weights_(i)*(zdiff*zdiff.transpose());
  }
  S(0, 0) +=std_radr_*std_radr_;
  S(1, 1) +=std_radphi_*std_radphi_;
  S(2, 2) +=std_radrd_*std_radrd_;


 //Compute cross correlation matrix T
  //compute Kalman gain matrix
  
  MatrixXd T = MatrixXd::Zero(n_x_,n_z);
  VectorXd xdiff;
//  VectorXd zdiff;
  for(int i = 0; i< n_sig; i++){
	xdiff = Xsig_pred_.col(i) - x_;
	zdiff = Zsig.col(i) - z_pred;
	//normalize angles
	zdiff(1) = AngleWrap(zdiff(1));
	xdiff(3) = AngleWrap(xdiff(3));

  	T += weights_(i)*xdiff*zdiff.transpose();
  }

  //Kalman gain matrix
  MatrixXd K = T*S.inverse();
  //update dates and covariance matrix!
  zdiff = zvec - z_pred;
	//normalize angles
  zdiff(1) = AngleWrap(zdiff(1));

  x_ += K*zdiff;
  P_ -= K*S*K.transpose();

  x_(3) = AngleWrap(x_(3));
//  cout<<"P: "<<P_<<endl;
//  cout<<"x: "<<x_<<endl;
//  cout<<"K "<< K <<endl;
//  cout<<"Xsig_pred: "<<Xsig_pred_<<endl;
  double NIS = zdiff.transpose()*S.inverse()*zdiff;
  cout<<"Radar NIS: "<< NIS<<endl;


}
