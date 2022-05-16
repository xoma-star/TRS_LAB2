#pragma once

#include <iostream>
#include <vector>
#include <functional>
#include <string>
#include <fstream>
#include <memory>
#include <chrono>

class DifferenceScheme
{
public:
	DifferenceScheme(
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
		int numberOfStepsT = -1
	);

	~DifferenceScheme();

	void drawGraph(std::vector<std::pair<std::string, std::string>>&& pathTofiles,
		std::string&& graphName, const std::string path, std::string&& labelX, std::string&& labelY, std::string&& labelZ);

	void drawGraphNorm(std::vector<std::pair<std::string, std::string>>&& pathTofiles,
		std::string&& graphName, const std::string path, std::string&& labelX, std::string&& labelY);


	void setNumberOfStepsSpace(unsigned numberOfStepsX);
	void setNumberOfStepsTimeAndSpace(unsigned numberOfStepsT, double numberOfStepsX);



protected:

	virtual double calculate(const bool isSaveGraph = false, const bool isSaveError = false) { return -1; };
	double error(const std::string& fileName = "");


	double* layers;


	std::function<double(double)> psi_0;
	std::function<double(double)> psi_1;
	std::function<double(double)> phi;
	std::function<double(double, double)> f;

	std::function<double(double, double)> solution;



	double alpha_0;
	double alpha_1;
	double beta_0;
	double beta_1;
	double k_0;
	double gridPitchX;
	double gridPitchT;

	unsigned numberOfStepsX;
	unsigned numberOfStepsT;

	std::string fileName;

	double T;

	void saveToFile(const std::string& fileName, double* array);

private:



};
