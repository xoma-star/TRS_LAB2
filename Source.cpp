#include <iostream>

# define M_PI           3.14159265358979323846

#include "ExplicitDifferenceScheme.h"
#include "FullyImplicitScheme.h"
#include "CrankNicholsonScheme.h"
#include "ConservativeScheme.h"

constexpr int k_0 = 1. / (4 * M_PI * M_PI);
constexpr int alpha_0 = 1;
constexpr int alpha_1 = 1;
constexpr int beta_0 = 0;
constexpr int beta_1 = 1;
constexpr double T = 1.;


constexpr unsigned numberOfStepsX = 100;

double f(double t, double x);
double phi(double x);
double psi_0(double t);
double psi_1(double t);


double K(double u);
double F(double u);



double solution(double t, double x);

#define pathToFiles "C:/Project/USATU_Lab/TRS_LAB/Lab_1/Lab_2"
//#define pathToFiles "C:/Users/hp/OneDrive/Документы/УГАТУ/3 курс/ТРС_2/Lab_2"

void TASK1(bool isDraw = false);
void TASK2(bool isDraw = false);
void TASK3(bool isDraw = false);

int main() {

	TASK1(false);
	TASK2(false);
	TASK3(false);


	return 0;
}

double f(double t, double x)
{
	return sin(2. * M_PI * x) - 1. / (2. * M_PI * M_PI);
}

double phi(double x)
{
	return x * x;
}

double psi_0(double t)
{
	return 2. * M_PI * (1 - exp(-t));
}

double psi_1(double t)
{
	return 1.;
}


double solution(double t, double x)
{
	return (1 - exp(-t)) * sin(2 * M_PI * x) + x * x;
}


void TASK1(bool isDraw)
{

	ExplicitDifferenceScheme explicitScheme(T, numberOfStepsX, alpha_0, alpha_1,
		beta_0, beta_1, k_0, psi_0, psi_1, phi, f, solution, "explicitSchem");

	std::cout << explicitScheme.calculate() << "\n\n";

	explicitScheme.setNumberOfStepsSpace(50);
	std::cout << explicitScheme.calculate() << "\n\n";

	explicitScheme.setNumberOfStepsSpace(100);
	std::cout << explicitScheme.calculate() << "\n\n";

	//std::ofstream fileToWrite("ExplicitSchemeError_X.txt");


	double error;

	std::cout << "\nSTART T = 20000, X = 10, 50, 100;\n";

	explicitScheme.setNumberOfStepsTimeAndSpace(20000, 10);
	std::cout << explicitScheme.calculate(true, true) << "\n\n";
	//error = explicitScheme.calculate();
	//std::cout << error << "\n\n";
	//fileToWrite << "10\t" << error << "\n";

	explicitScheme.setNumberOfStepsTimeAndSpace(20000, 50);
	std::cout << explicitScheme.calculate(true, true) << "\n\n";
	//error = explicitScheme.calculate();
	//std::cout << error << "\n\n";
	//fileToWrite << "50\t" << error << "\n";

	explicitScheme.setNumberOfStepsTimeAndSpace(20000, 100);
	std::cout << explicitScheme.calculate() << "\n\n";
	//error = explicitScheme.calculate();
	//std::cout << error << "\n\n";
	//fileToWrite << "100\t" << error << "\n";

	//fileToWrite.close();

	std::cout << "\nEND T = 20000, X = 10, 50, 100;\n";



	//fileToWrite.open("ExplicitSchemeError_T.txt");

	std::cout << "\nSTART T = 200, 5000, 20000, X = 10\n";

	explicitScheme.setNumberOfStepsTimeAndSpace(200, 10);
	std::cout << explicitScheme.calculate(true, true) << "\n\n";
	//error = explicitScheme.calculate();
	//std::cout << error << "\n\n";
	//fileToWrite << "200\t" << error << "\n";

	explicitScheme.setNumberOfStepsTimeAndSpace(5000, 10);
	std::cout << explicitScheme.calculate(true, true) << "\n\n";
	//error = explicitScheme.calculate(true, true);
	//std::cout << error << "\n\n";
	//fileToWrite << "5000\t" << error << "\n";

	explicitScheme.setNumberOfStepsTimeAndSpace(20000, 10);
	std::cout << explicitScheme.calculate() << "\n\n";
	//error = explicitScheme.calculate();
	//std::cout << error << "\n\n";
	//fileToWrite << "20000\t" << error << "\n";

	//fileToWrite.close();


	std::cout << "\nEND T = 200, 5000, 20000, X = 10\n";



	explicitScheme.setNumberOfStepsTimeAndSpace(10, 10);
	std::cout << explicitScheme.calculate() << "\n\n";

	if (isDraw) {
		explicitScheme.drawGraph({ {"explicitSchem_10_20000Error.txt", "10 steps - 20000 times"} }, "temp",
			pathToFiles, "X", "T", "error");
		explicitScheme.drawGraph({ {"explicitSchem_50_20000Error.txt", "50 steps - 20000 times"} }, "temp",
			pathToFiles, "X", "T", "error");

		explicitScheme.drawGraph({ {"explicitSchem_10_200Error.txt", "10 steps - 200 times"} }, "temp",
			pathToFiles, "X", "T", "error");
		explicitScheme.drawGraph({ {"explicitSchem_10_5000Error.txt", "10 steps - 5000 times"} }, "temp",
			pathToFiles, "X", "T", "error");

		//explicitScheme.drawGraph({ {"explicitSchem_10_5000Error.txt", "10 steps - 5000 times"} }, "temp",
		//	pathToFiles, "X", "T", "error");

		//explicitScheme.drawGraphNorm({ {"ExplicitSchemeError_X.txt", "Change X; Error"} }, "t", pathToFiles, "T", "Error");
		//explicitScheme.drawGraphNorm({ {"ExplicitSchemeError_T.txt", "Change T; Error"} }, "t", pathToFiles, "T", "Error");

	}
}

void TASK2(bool isDraw)
{
	std::cout << "\n\n\----------------------\n\n\n";
	std::cout << "fullyImplicitSchema\n\n";

	FullyImplicitScheme fullyImplicitScheme(T, numberOfStepsX, alpha_0, alpha_1,
		beta_0, beta_1, k_0, psi_0, psi_1, phi, f, solution, "fullyImplicitSchema");

	std::cout << fullyImplicitScheme.calculate() << "\n\n";

	fullyImplicitScheme.setNumberOfStepsSpace(50);
	std::cout << fullyImplicitScheme.calculate(true, true) << "\n\n";

	fullyImplicitScheme.setNumberOfStepsSpace(100);
	std::cout << fullyImplicitScheme.calculate() << "\n\n";


	std::cout << "\n\n -------------------------------- \n";

	std::cout << "crankNicholsonScheme\n\n";
	CrankNicholsonScheme crankNicholsonScheme(T, numberOfStepsX, alpha_0, alpha_1,
		beta_0, beta_1, k_0, psi_0, psi_1, phi, f, solution, "crankNicholsonScheme");

	std::cout << crankNicholsonScheme.calculate() << "\n\n";


	crankNicholsonScheme.setNumberOfStepsSpace(50);
	std::cout << crankNicholsonScheme.calculate(true, true) << "\n\n";

	crankNicholsonScheme.setNumberOfStepsSpace(100);
	std::cout << crankNicholsonScheme.calculate() << "\n\n";

	if (isDraw) {
		fullyImplicitScheme.drawGraph({ {"fullyImplicitSchema_50_5000Error.txt", "FI 50 steps - 5000 times"}
			}, "temp",
			pathToFiles, "X", "T", "error");

		crankNicholsonScheme.drawGraph({ {"crankNicholsonScheme_50_5000Error.txt", "CN 50 steps - 5000 times"}
			}, "temp",
			pathToFiles, "X", "T", "error");
	}

}

void TASK3(bool isDraw)
{
	std::cout << "\n\n -------------------------------- \n";

	std::cout << "CHECK SCHEME ON TASK_1 START\n\n";

	ConservativeScheme conservativeSchemeCHECK(T, numberOfStepsX, alpha_0, alpha_1,
		beta_0, beta_1, k_0, psi_0, psi_1, phi, f, solution, [](double u) {return 1; }, [](double u) {return 1; }, "conservativeSchemeCHECK");

	std::cout << conservativeSchemeCHECK.calculate() << "\n";

	conservativeSchemeCHECK.setNumberOfStepsSpace(50);
	std::cout << conservativeSchemeCHECK.calculate() << "\n";

	conservativeSchemeCHECK.setNumberOfStepsSpace(100);
	std::cout << conservativeSchemeCHECK.calculate() << "\n";

	std::cout << "\nCHECK SCHEME ON TASK_1 END\n\n\n";

	std::cout << "START TASK_3. Conservative Scheme\n\n";

	ConservativeScheme conservativeScheme(T, numberOfStepsX, alpha_0, alpha_1,
		beta_0, beta_1, k_0, psi_0, psi_1, phi, f, solution, F, K, "conservativeScheme");

	std::cout << conservativeScheme.calculate() << "\n";

	conservativeScheme.setNumberOfStepsSpace(50);
	std::cout << conservativeScheme.calculate(true) << "\n";

	conservativeScheme.setNumberOfStepsSpace(100);
	std::cout << conservativeScheme.calculate() << "\n";

	std::cout << "END TASK_3. Conservative Scheme\n\n\n";

	if (isDraw) {
		conservativeScheme.drawGraph({ {"conservativeScheme_50_5000Solution.txt", "conservative 50 steps - 5000 times"}
			}, "temp",
			pathToFiles, "X", "T", "u");
	}

	//conservativeScheme.drawGraph({ {"fullyImplicitSchema_10_200Error.txt", "FI"}, {"crankNicholsonScheme_10_200Error.txt", "CN"} }, "temp",
		//pathToFiles, "X", "T", "error");

}

double K(double u)
{
	return u;
	//return 1.;

}

double F(double u)
{
	return u * u * u;
	//return 1.;

}
