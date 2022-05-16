#pragma once
#include "ImplicitScheme.h"

class ConservativeScheme : public ImplicitScheme
{
public:
	ConservativeScheme(
		double T,
		unsigned numberOfStepsX,
		double alpha_0,
		double alpha_1,
		double beta_0,
		double beta_1,
		double k_0,
		std::function<double(double)> psi_0,
		std::function<double(double)> psi_1,
		std::function<double(double)> phi,
		std::function<double(double, double)> f,
		std::function<double(double, double)> solution,
		std::function<double(double)> F,
		std::function<double(double)> k,
		const std::string& path,
		int gridPitchT = -1.);
	//~ConservativeScheme();


	double virtual calculate(const bool isSaveGraph = false, const bool isSave = false) override;

private:
	std::function<double(double)> F;
	std::function<double(double)> k;


	double function(double t, double x, double u);
};
#pragma once
