#include "FullyImplicitScheme.h"

FullyImplicitScheme::FullyImplicitScheme(
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
	const std::string& path) : ImplicitScheme(
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
		path
	)
{
}


double FullyImplicitScheme::calculate(const bool isSaveGraph, const bool isSaveError)
{
	auto start = std::chrono::steady_clock::now();

	double* elementsMatrix = new double[(numberOfStepsX + 1) * 4]; // 0 - верх диаг.; 1 - ср. диаг.; 2 - нижн. диаг.;

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
			//верхн.
			elementsMatrix[0 * (numberOfStepsX + 1) + x] = k_0 * gridPitchT / (gridPitchX * gridPitchX);

			//ср.
			elementsMatrix[1 * (numberOfStepsX + 1) + x] = -2. * gridPitchT * k_0 / (gridPitchX * gridPitchX) - 1;

			//нижн;
			elementsMatrix[2 * (numberOfStepsX + 1) + x - 1] = k_0 * gridPitchT / (gridPitchX * gridPitchX);

			//прав;
			elementsMatrix[3 * (numberOfStepsX + 1) + x] = -layers[t * (numberOfStepsX + 1) + x] - gridPitchT * f(gridPitchT * (t + 1), gridPitchX * (x));
		}

		progonka(elementsMatrix, t + 1);
	}

	delete[] elementsMatrix;

	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

	if (isSaveError) {

		return error(fileName);
	}
	else {
		return error();
	}
}
