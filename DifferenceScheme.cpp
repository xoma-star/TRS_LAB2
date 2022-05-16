#include "DifferenceScheme.h"

DifferenceScheme::DifferenceScheme(
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
	int numberOfStepsT) : alpha_0{ alpha_0 },
	T{ T },
	numberOfStepsX{ numberOfStepsX },
	alpha_1{ alpha_1 },
	beta_0{ beta_0 },
	beta_1{ beta_1 },
	k_0{ k_0 },
	psi_0{ psi_0 },
	psi_1{ psi_1 },
	phi{ phi },
	f{ f },
	fileName{ path },
	solution{ solution }
{
	if (numberOfStepsT == -1) {
		this->numberOfStepsT = 2. * numberOfStepsX * numberOfStepsX * T;
	}
	else {
		this->numberOfStepsT = numberOfStepsT;
	}
	gridPitchX = 1. / numberOfStepsX;
	gridPitchT = 1. / this->numberOfStepsT;
	layers = new double[(this->numberOfStepsT + 1) * (numberOfStepsX + 1)];

	//std::cout << "X " << numberOfStepsX << "\n\n";
	//std::cout << "T " << numberOfStepsT << "\n\n";
	std::cout << "h = " << gridPitchX << ";  t = " << gridPitchT << "\n";
	//std::cout << (numberOfStepsT + 1) * (numberOfStepsX + 1) << "\n\n";
}

void DifferenceScheme::setNumberOfStepsSpace(unsigned numberOfStepsX)
{
	this->numberOfStepsX = numberOfStepsX;

	numberOfStepsT = 2. * numberOfStepsX * numberOfStepsX;
	gridPitchX = 1. / numberOfStepsX;
	gridPitchT = 1. / numberOfStepsT;

	delete[] layers;
	std::cout << "SET  h = " << gridPitchX << ";  t = " << gridPitchT << "\n";
	layers = new double[(numberOfStepsT + 1) * (numberOfStepsX + 1)];
}

void DifferenceScheme::setNumberOfStepsTimeAndSpace(unsigned numberOfStepsT, double numberOfStepsX)
{
	this->numberOfStepsX = numberOfStepsX;
	this->numberOfStepsT = numberOfStepsT;
	gridPitchT = 1. / numberOfStepsT;
	gridPitchX = 1. / numberOfStepsX;

	delete[] layers;
	std::cout << "SET  h = " << gridPitchX << ";  t = " << gridPitchT << "\n";
	layers = new double[(numberOfStepsT + 1) * (numberOfStepsX + 1)];

}

DifferenceScheme::~DifferenceScheme()
{
	delete[] layers;
}

double DifferenceScheme::error(const std::string& fileName)
{

	std::unique_ptr<double[]> error{ new double[(numberOfStepsT + 1) * (numberOfStepsX + 1)] };

	double max = -1.;


	for (unsigned t = 0; t <= numberOfStepsT; ++t) {
		for (unsigned x = 0; x <= numberOfStepsX; ++x) {
			error[t * (numberOfStepsX + 1) + x] = fabs(layers[t * (numberOfStepsX + 1) + x] - solution(gridPitchT * t, gridPitchX * x));
			if (error[t * (numberOfStepsX + 1) + x] >= max) {
				max = fabs(layers[t * (numberOfStepsX + 1) + x] - solution(gridPitchT * t, gridPitchX * x));
			}

		}
	}
	if (!fileName.empty()) {
		saveToFile(fileName + "_" + std::to_string(numberOfStepsX) + "_" + std::to_string(numberOfStepsT) + "Error.txt", error.get());
	}

	return max;
}

void DifferenceScheme::drawGraph(std::vector<std::pair<std::string, std::string>>&& pathTofiles,
	std::string&& graphName, const std::string path, std::string&& labelX, std::string&& labelY, std::string&& labelZ)
{
	std::string command = "";
	std::string lineUsing = "\" u 1:2:3";
	std::string lineName = " title ";
	std::string lineWith = " with lines";

	command += "set ticslevel 0\n";
	command += "splot";
	command += "\"" + pathTofiles[0].first + lineUsing + lineName + "\"n = " + pathTofiles[0].second + "\"" + lineWith;

	for (unsigned i = 1; i < pathTofiles.size(); ++i) {
		command += ", ";

		command += "\"" + pathTofiles[i].first + lineUsing + lineName + "\"n = " + pathTofiles[i].second + "\"" + lineWith;

	}
	command += ";set term wxt title \"" + graphName + "\"";
	command += ";set xlabel \"" + labelX + "\"";
	command += ";set ylabel \"" + labelY + "\"";
	command += ";set zlabel \"" + labelZ + "\"";

	std::ofstream graphic("file");
	graphic << "cd \"" + path + "\"" << "\n";
	graphic << command << "; pause mouse keypress" << "\n";
	graphic.close();
	std::system("gnuplot -persist file");

	graphic.close();
}


void DifferenceScheme::drawGraphNorm(std::vector<std::pair<std::string, std::string>>&& pathTofiles,
	std::string&& graphName, const std::string path, std::string&& labelX, std::string&& labelY)
{
	std::string lineUsing = "\" using ";
	std::string lineName = " title ";
	std::string lineWith = " with lines ls ";
	std::string ls = " lw 3 ";

	int numLine = 1;
	std::string command = "";

	std::string counter;
	command += "\"" + pathTofiles[0].first + lineUsing + "1:2" + lineName + "\"n = " + pathTofiles[0].second + "\"" + lineWith + "1" + ls;
	for (unsigned i = 1; i < pathTofiles.size(); ++i) {
		command += ", ";
		counter = std::to_string(numLine + 1);
		command += "\"" + pathTofiles[i].first + lineUsing + "1:2" + lineName + "\"n = " + pathTofiles[i].second + "\"" + lineWith + counter + ls;

		numLine++;
	}
	command += ";set term wxt title \"" + graphName + "\"";
	command += ";set xlabel \"" + labelX + "\"";
	command += ";set ylabel \"" + labelY + "\"";

	std::ofstream graphic("file");
	graphic << "cd \"" + path + "\"" << "\n";
	graphic << "plot " << command << "; pause mouse keypress" << "\n";
	graphic.close();
	std::system("gnuplot -persist file");

	graphic.close();
}

//const std::unique_ptr<double[]> array
void DifferenceScheme::saveToFile(const std::string& fileName, double* array)
{
	std::ofstream fileToWrite(fileName);

	for (unsigned t = 0; t <= numberOfStepsT; ++t) {
		for (unsigned x = 0; x <= numberOfStepsX; ++x) {
			fileToWrite << x * gridPitchX << "\t" << t * gridPitchT << "\t" << array[t * (numberOfStepsX + 1) + x] << "\n";
		}
	}
	fileToWrite.close();


}
