#pragma once

#include "ImplicitScheme.h"

class FullyImplicitScheme : public ImplicitScheme
{
public:
	FullyImplicitScheme(
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
		const std::string& path);
	//~FullyImplicitSchema();

	virtual double calculate(const bool isSaveGraph = false, const bool isSaveError = false) override;

private:

};
