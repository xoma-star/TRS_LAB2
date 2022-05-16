#include "ConservativeScheme.h"

ConservativeScheme::ConservativeScheme(
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
	int gridPitchT) : ImplicitScheme(
		T,
		numberOfStepsX,
		alpha_0,
		alpha_1,
		beta_0, beta_1,
		k_0,
		psi_0,
		psi_1,
		phi,
		f,
		solution,
		path,
		gridPitchT
	)
{
	this->F = F;
	this->k = k;
}

double ConservativeScheme::calculate(const bool isSaveGraph, const bool isSaveError)
{
	double* elementsMatrix = new double[(numberOfStepsX + 1) * 4]; // 0 - верх диаг.; 1 - ср. диаг.; 2 - нижн. диаг.;
	double k_e;

	double coeff = 1.;

	for (unsigned i = 0; i <= numberOfStepsX; ++i) {
		layers[0 * (numberOfStepsT + 1) + i] = phi(i * gridPitchX);
	}

	for (unsigned t = 0; t < numberOfStepsT; ++t) {

		//верхн.
		elementsMatrix[0 * (numberOfStepsX + 1) + 0] = beta_0;

		//ср.
		elementsMatrix[1 * (numberOfStepsX + 1) + 0] = gridPitchX * alpha_0 - beta_0;
		elementsMatrix[1 * (numberOfStepsX + 1) + numberOfStepsX] = alpha_1 * gridPitchX + beta_1;

		//нижн; 
		elementsMatrix[2 * (numberOfStepsX + 1) + numberOfStepsX - 1] = -beta_1;

		//прав;
		elementsMatrix[3 * (numberOfStepsX + 1) + 0] = gridPitchX * psi_0(gridPitchT * (t + 1));
		elementsMatrix[3 * (numberOfStepsX + 1) + numberOfStepsX] = gridPitchX * psi_1(gridPitchT * (t + 1));



		for (unsigned x = 1; x < numberOfStepsX; ++x) {
			double V_plus = k((layers[t * (numberOfStepsX + 1) + x] + layers[t * (numberOfStepsX + 1) + x + 1]) / 2.);
			double V_minus = k((layers[t * (numberOfStepsX + 1) + x] + layers[t * (numberOfStepsX + 1) + x - 1]) / 2.);

			//double V_plus = ( k(layers[t * (numberOfStepsX + 1) + x]) + k(layers[t * (numberOfStepsX + 1) + x + 1]))  / 2.;
			//double V_minus = ( k(layers[t * (numberOfStepsX + 1) + x]) + k(layers[t * (numberOfStepsX + 1) + x - 1]) ) / 2.;


			//верхн.
			elementsMatrix[0 * (numberOfStepsX + 1) + x] = gridPitchT * coeff / (gridPitchX * gridPitchX) * V_plus;


			//ср.
			elementsMatrix[1 * (numberOfStepsX + 1) + x] = -1. - gridPitchT * coeff / (gridPitchX * gridPitchX) * (V_plus + V_minus);


			//нижн;
			elementsMatrix[2 * (numberOfStepsX + 1) + x - 1] = gridPitchT * coeff / (gridPitchX * gridPitchX) * V_minus;


			//прав;
			elementsMatrix[3 * (numberOfStepsX + 1) + x] = -layers[t * (numberOfStepsX + 1) + x] -
				gridPitchT * (1 - coeff) / gridPitchX * (V_plus * (layers[t * (numberOfStepsX + 1) + x + 1] - layers[t * (numberOfStepsX + 1) + x]) / gridPitchX -
					V_minus * (layers[t * (numberOfStepsX + 1) + x] - layers[t * (numberOfStepsX + 1) + x - 1]) / gridPitchX) -
				gridPitchT * function(gridPitchT * t, gridPitchX * x, layers[t * (numberOfStepsX + 1) + x + 1]);
		}

		progonka(elementsMatrix, t + 1);
	}

	delete[] elementsMatrix;

	if (isSaveGraph) {
		saveToFile(fileName + "_" + std::to_string(numberOfStepsX) + "_" + std::to_string(numberOfStepsT) + "Solution.txt", layers);
	}

	if (isSaveError) {
		return error(fileName);
	}
	else {
		return error();
	}
}

double ConservativeScheme::function(double t, double x, double u)
{
	return F(u) * f(t, x);
}
