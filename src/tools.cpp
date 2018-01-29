#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE here.
  */
	
	VectorXd ret_RMSE = VectorXd::Zero(4);
	if(estimations.size() != ground_truth.size() || estimations.size() == 0 ){
		std:cout<<"Error in estimation or ground_truth vector length"<<std::endl;
		return ret_RMSE;
	}
		
		
	for(int i = 0; i< estimations.size(); i++)
		ret_RMSE +=  VectorXd((estimations[i] - ground_truth[i]).array()*(estimations[i] - ground_truth[i]).array());

	ret_RMSE /= estimations.size();
	ret_RMSE = ret_RMSE.array().sqrt();
	return ret_RMSE;
	
}
