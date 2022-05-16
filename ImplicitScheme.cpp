#include "ImplicitScheme.h"

ImplicitScheme::ImplicitScheme(
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
	const std::string& path,
	int numberOfStepsT
) : DifferenceScheme(
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
	numberOfStepsT
)
{

}

void ImplicitScheme::progonka(double* matrix, unsigned time)
{
	for (unsigned x = 1; x <= numberOfStepsX; x++)
	{
		matrix[1 * (numberOfStepsX + 1) + x] = matrix[1 * (numberOfStepsX + 1) + x] -
			matrix[0 * (numberOfStepsX + 1) + x - 1] *
			(matrix[2 * (numberOfStepsX + 1) + x - 1] /
				matrix[1 * (numberOfStepsX + 1) + x - 1]);

		matrix[3 * (numberOfStepsX + 1) + x] = matrix[3 * (numberOfStepsX + 1) + x] -
			matrix[3 * (numberOfStepsX + 1) + x - 1] *
			(matrix[2 * (numberOfStepsX + 1) + x - 1] /
				matrix[1 * (numberOfStepsX + 1) + x - 1]);

		matrix[2 * (numberOfStepsX + 1) + x - 1] = 0;
	}

	layers[(numberOfStepsX + 1) * time + numberOfStepsX] = matrix[3 * (numberOfStepsX + 1) + numberOfStepsX] / matrix[1 * (numberOfStepsX + 1) + numberOfStepsX];
	for (int x = numberOfStepsX - 1; x >= 0; --x)
		layers[(numberOfStepsX + 1) * time + x] = (matrix[3 * (numberOfStepsX + 1) + x] -
			matrix[0 * (numberOfStepsX + 1) + x] *
			layers[(numberOfStepsX + 1) * time + x + 1]) /
		matrix[1 * (numberOfStepsX + 1) + x];
}
