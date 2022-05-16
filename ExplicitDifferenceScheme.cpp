#include "ExplicitDifferenceScheme.h"

ExplicitDifferenceScheme::ExplicitDifferenceScheme(
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
	const std::string& path) : DifferenceScheme(
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

double ExplicitDifferenceScheme::calculate(const bool isSaveGraph, const bool isSaveError)
{
	/*std::cout << gridPitchX << "\t" << gridPitchT << "\n";*/
	auto start = std::chrono::steady_clock::now();
	for (unsigned i = 0; i <= numberOfStepsX; ++i) {
		layers[0 * (numberOfStepsT + 1) + i] = phi(i * gridPitchX);
	}

	for (unsigned t = 0; t < numberOfStepsT; ++t) {
		for (unsigned x = 1; x < numberOfStepsX; ++x) {
			layers[(t + 1) * (numberOfStepsX + 1) + x] = layers[t * (numberOfStepsX + 1) + x] +
				k_0 * gridPitchT / (gridPitchX * gridPitchX) *
				(layers[(t) * (numberOfStepsX + 1) + x - 1] - 2. * layers[t * (numberOfStepsX + 1) + x] +
					layers[t * (numberOfStepsX + 1) + x + 1]) +
				gridPitchT * f(t * gridPitchT, x * gridPitchX);



			//std::cout << "1 " << layers[(t + 1) * (numberOfStepsX + 1) + 0] << "\n";
			//std::cout << "2 " << layers[(t + 1) * (numberOfStepsX + 1) + 1] << "\n";

			layers[(t + 1) * (numberOfStepsX + 1) + 0] = (psi_0((t + 1) * gridPitchT) * gridPitchX - beta_0 * layers[(t + 1) * (numberOfStepsX + 1) + 1]) /
				(alpha_0 * gridPitchX - beta_0);

			layers[(t + 1) * (numberOfStepsX + 1) + numberOfStepsX] = (psi_1((t + 1) * gridPitchT) * gridPitchX +
				beta_1 * layers[(t + 1) * (numberOfStepsX + 1) + numberOfStepsX - 1]) / (alpha_1 * gridPitchX + beta_1);


		}
	}

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
