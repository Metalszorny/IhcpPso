// IhcpPso1dCPlusPlusMultiThreaded.cpp : Defines the entry point for the console application.
//

#include "stdafx.h";
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <algorithm>
#include <time.h>
#include <iomanip>
#include <iterator>
#include <thread>

#pragma region Fields

double _initialTemperature;
int _horizontalSplitting;
int _verticalSplitting;
int _monitoredIndex;
double _timeDifference;
int _dimensionNumber;
double* _htcValues;
double* _specificHeats;
int _specificHeatsSizeRows;
double* _heatConductivities;
int _heatConductivitiesSizeRows;
int _rangeMin;
int _rangeMax;
double _simulationTime;
double _height;
double _radius;
double _density;
double _dX;
double _dY;
double* _referenceCooldownCurve;
int _referenceCooldownCurveSizeRows;
double _tK;
double* _particleInitialVelocity;
int _particleInitialVelocitySize;
double* _particleInitialPosition;
int _particleInitialPositionSize;
int _particleInformerNumber;
double _particleConstant1;
double _particleConstant2;
double _particleEpsilon;
double* _particlePosition;
double* _particleVelocity;
double* _particleBestPosition;
double* _particleFitness;
double* _particleBestFitness;
double* _globalBestFitness;
double* _globalBestPosition;
int _globalBestSize;
int _particleOptimalisationType;
int _particleSwarmSize;
int _maxEpochs;
int _maxStaticEpochs;
double _weight;
int* _particleInformers;
std::string _exitReason;
int _epoch;

#pragma endregion

#pragma region Methods

void ReadDataFromFile(std::string dataFilePath);
void CalculateReferenceCooldownCurve();
double* CalculateCooldownCurve1D(bool isParticle, int particleIndex, double* currentTemperature, double* previousTemperature, double* g);
double GetHeatConductivity(double temperature);
double GetSpecificHeat(double temperature);
double GetAlpha(double heatConductivity, double specificHeat);
double GetHeatTransferCoefficient(double temperature, bool isParticle, int particleIndex);
void WriteReferenceCooldownLogToFile();
void WritePsoGlobalBestLogToFile();
void WriteExitResultAndTimeLogToFile(long elaspedMilliseconds);
void InitializeParticle(int particleIndex);
double ObjectiveFunction(int particleIndex);
void Optimize();
void UpdateRing();
void OptimizePosition(int particleIndex, double mbest, double beta);
void UpdateClercPosition(int particleIndex);
void UpdateErrorValues(int particleIndex);
void UpdatePosition(int particleIndex);
void UpdateQuantumPosition(int particleIndex, double mbest, double beta);
void UpdateVelocity(int particleIndex, double* bestLocalPosition);
bool CheckIsInRange(double* position);
double* GetBestLocalPositions(int particleIndex);
void GetInitialPosition(int particalIndex);
void GetInitialVelocity(int particalIndex);
int* GetIntegerRange(int startIndex, int count);
int* Shuffle(int* index, int indexSize);
double GetSwarmAverageBestPosition();
std::vector<std::string> explode(std::string originalString, char delimeter);
template<typename Out> void split(std::string s, char delim, Out result);
double* SetCurrentTemperature();
double* SetPreviousTemperature();
double* SetG();

int main()
{
	_exitReason = "";
	ReadDataFromFile("ConfigurationIn.txt");
	CalculateReferenceCooldownCurve();
	WriteReferenceCooldownLogToFile();

	// Get particle swarm initial values.
	_particlePosition = new double[_particleSwarmSize * _dimensionNumber];
	_particleVelocity = new double[_particleSwarmSize * _dimensionNumber];
	_particleBestPosition = new double[_particleSwarmSize * _dimensionNumber];
	_particleFitness = new double[_particleSwarmSize];
	_particleBestFitness = new double[_particleSwarmSize];
	_globalBestSize = 1;
	_globalBestFitness = new double[_globalBestSize];
	_globalBestFitness[_globalBestSize - 1] = DBL_MAX;
	_globalBestPosition = new double[_globalBestSize * _dimensionNumber];
	std::vector<std::thread> taskArray;

	// Set particle swarm initial values.
	for (int i = 0; i < _particleSwarmSize; i += 1)
	{
		taskArray.push_back(std::thread(InitializeParticle, i));
	}

	for (int i = 0; i < _particleSwarmSize; i += 1)
	{
		taskArray[i].join();
	}

	// Set informers.
	UpdateRing();

	// Set initial best values.
	for (int i = 0; i < _particleSwarmSize; i += 1)
	{
		if (_particleBestFitness[i] < _globalBestFitness[_globalBestSize - 1])
		{
			_globalBestFitness[_globalBestSize - 1] = _particleBestFitness[i];

			for (int j = 0; j < _dimensionNumber; j += 1)
			{
				_globalBestPosition[(_globalBestSize - 1) + j] = _particlePosition[(i * _dimensionNumber) + j];
			}
		}
	}

	clock_t start = clock();
	Optimize();
	clock_t finish = clock();
	WritePsoGlobalBestLogToFile();
	WriteExitResultAndTimeLogToFile((long)(finish - start));

	return 0;
}

void ReadDataFromFile(std::string dataFilePath)
{
	std::ifstream fileStream(dataFilePath);

	if (fileStream.is_open())
	{
		std::string currentLine = "";

		while (std::getline(fileStream, currentLine))
		{
			if (currentLine.find("#") != std::string::npos)
			{
				std::transform(currentLine.begin(), currentLine.end(), currentLine.begin(), ::tolower);

				if (currentLine.find("initial temperature") != std::string::npos)
				{
					std::getline(fileStream, currentLine);
					std::remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_initialTemperature = atof(currentLine.c_str());
				}
				else if (currentLine.find("horizontal splitting") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_horizontalSplitting = atoi(currentLine.c_str());
				}
				else if (currentLine.find("vertical splitting") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_verticalSplitting = atoi(currentLine.c_str());
				}
				else if (currentLine.find("monitored index") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_monitoredIndex = atoi(currentLine.c_str());
				}
				else if (currentLine.find("range min") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_rangeMin = atoi(currentLine.c_str());
				}
				else if (currentLine.find("range max") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_rangeMax = atoi(currentLine.c_str());
				}
				else if (currentLine.find("simulation time") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_simulationTime = atof(currentLine.c_str());
				}
				else if (currentLine.find("height") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_height = atof(currentLine.c_str());
				}
				else if (currentLine.find("htcin") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					dataFilePath = currentLine;
				}
				else if (currentLine.find("time difference") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_timeDifference = atof(currentLine.c_str());
				}
				else if (currentLine.find("pso epsilon") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_particleEpsilon = atof(currentLine.c_str());
				}
				else if (currentLine.find("radius") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_radius = atof(currentLine.c_str());
				}
				else if (currentLine.find("density") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_density = atof(currentLine.c_str());
				}
				else if (currentLine.find("optimization type") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');

					if (currentLine == "clerc")
					{
						_particleOptimalisationType = 1;
					}
					else if (currentLine == "quantum")
					{
						_particleOptimalisationType = 2;
					}
					else
					{
						_particleOptimalisationType = 1;
					}
				}
				else if (currentLine.find("swarm size") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_particleSwarmSize = atoi(currentLine.c_str());
				}
				else if (currentLine.find("max epochs") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_maxEpochs = atoi(currentLine.c_str());
				}
				else if (currentLine.find("max static epochs") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_maxStaticEpochs = atoi(currentLine.c_str());
				}
				else if (currentLine.find("weight") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_weight = atof(currentLine.c_str());
				}
				else if (currentLine.find("constant1") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_particleConstant1 = atof(currentLine.c_str());
				}
				else if (currentLine.find("constant2") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_particleConstant2 = atof(currentLine.c_str());
				}
				else if (currentLine.find("informer number") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_particleInformerNumber = atoi(currentLine.c_str());
				}
				else if (currentLine.find("initial position") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_particleInitialPositionSize = atoi(currentLine.c_str());

					if (_particleInitialPositionSize > 0)
					{
						_particleInitialPosition = new double[_particleInitialPositionSize];

						for (int i = 0; i < _particleInitialPositionSize; i += 1)
						{
							std::getline(fileStream, currentLine);
							remove(currentLine.begin(), currentLine.end(), ' ');
							replace(currentLine.begin(), currentLine.end(), ',', '.');
							_particleInitialPosition[i] = atof(currentLine.c_str());
						}
					}
				}
				else if (currentLine.find("initial velocity") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_particleInitialVelocitySize = atoi(currentLine.c_str());

					if (_particleInitialVelocitySize > 0)
					{
						_particleInitialVelocity = new double[_particleInitialVelocitySize];

						for (int i = 0; i < _particleInitialVelocitySize; i += 1)
						{
							std::getline(fileStream, currentLine);
							remove(currentLine.begin(), currentLine.end(), ' ');
							replace(currentLine.begin(), currentLine.end(), ',', '.');
							_particleInitialVelocity[i] = atof(currentLine.c_str());
						}
					}
				}
				else if (currentLine.find("heat conductivity") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_heatConductivitiesSizeRows = atoi(currentLine.c_str());

					if (_heatConductivitiesSizeRows > 0)
					{
						_heatConductivities = new double[_heatConductivitiesSizeRows * 2];

						for (int i = 0; i < (_heatConductivitiesSizeRows * 2); i += 2)
						{
							std::getline(fileStream, currentLine);
							replace(currentLine.begin(), currentLine.end(), ',', '.');
							replace(currentLine.begin(), currentLine.end(), '\t', ' ');
							std::vector<std::string> result = explode(currentLine, ' ');
							_heatConductivities[i + 0] = atof(result[0].c_str());
							_heatConductivities[i + 1] = atof(result[1].c_str());
						}
					}
				}
				else if (currentLine.find("specific heat") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_specificHeatsSizeRows = atoi(currentLine.c_str());

					if (_specificHeatsSizeRows > 0)
					{
						_specificHeats = new double[_specificHeatsSizeRows * 2];

						for (int i = 0; i < (_specificHeatsSizeRows * 2); i += 2)
						{
							std::getline(fileStream, currentLine);
							replace(currentLine.begin(), currentLine.end(), ',', '.');
							replace(currentLine.begin(), currentLine.end(), '\t', ' ');
							std::vector<std::string> result = explode(currentLine, ' ');
							_specificHeats[i + 0] = atof(result[0].c_str());
							_specificHeats[i + 1] = atof(result[1].c_str());
						}
					}
				}
			}
		}

		fileStream.close();
	}

	std::ifstream fileStream2(dataFilePath);

	if (fileStream2.is_open())
	{
		std::string currentLine = "";

		while (std::getline(fileStream2, currentLine))
		{
			if (currentLine.find("#") != std::string::npos)
			{
				std::transform(currentLine.begin(), currentLine.end(), currentLine.begin(), ::tolower);

				if (currentLine.find("htc values") != std::string::npos)
				{
					std::getline(fileStream2, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_dimensionNumber = atoi(currentLine.c_str());

					if (_dimensionNumber > 0)
					{
						_htcValues = new double[_dimensionNumber * 2];

						for (int i = 0; i < (_dimensionNumber * 2); i += 2)
						{
							std::getline(fileStream2, currentLine);
							replace(currentLine.begin(), currentLine.end(), ',', '.');
							std::vector<std::string> result = explode(currentLine, ' ');
							_htcValues[i + 0] = atof(result[0].c_str());
							_htcValues[i + 1] = atof(result[1].c_str());
						}
					}
				}
			}
		}

		fileStream2.close();
	}

	_dX = (_radius / 1000) / _horizontalSplitting;
	_dY = (_height / 1000) / _verticalSplitting;
	_tK = 0;
}

double* SetCurrentTemperature()
{
	double* returnValue = new double[_horizontalSplitting];

	for (int i = 0; i < _horizontalSplitting; i += 1)
	{
		returnValue[i] = _initialTemperature;
	}

	return returnValue;
}

double* SetPreviousTemperature()
{
	double* returnValue = new double[_horizontalSplitting];

	for (int i = 0; i < _horizontalSplitting; i += 1)
	{
		returnValue[i] = _initialTemperature;
	}

	return returnValue;
}

double* SetG()
{
	double* returnValue = new double[_horizontalSplitting];

	for (int i = 0; i < _horizontalSplitting; i += 1)
	{
		returnValue[i] = 0;
	}

	return returnValue;
}

void CalculateReferenceCooldownCurve()
{
	_referenceCooldownCurveSizeRows = (int)(_simulationTime / _timeDifference) + 1;
	_referenceCooldownCurve = new double[_referenceCooldownCurveSizeRows * 2];
	double* currentTemperature = SetCurrentTemperature();
	double* previousTemperature = SetPreviousTemperature();
	double* g = SetG();

	for (int i = 0; i < (_referenceCooldownCurveSizeRows * 2); i += 2)
	{
		currentTemperature = CalculateCooldownCurve1D(false, -1, currentTemperature, previousTemperature, g);
		_referenceCooldownCurve[i + 0] = ((i / 2) * _timeDifference);
		_referenceCooldownCurve[i + 1] = currentTemperature[_monitoredIndex];

		for (int j = 0; j < _horizontalSplitting; j += 1)
		{
			previousTemperature[j] = currentTemperature[j];
		}
	}

	delete[] currentTemperature;
	delete[] previousTemperature;
	delete[] g;
}

double* CalculateCooldownCurve1D(bool isParticle, int particleIndex, double* currentTemperature, double* previousTemperature, double* g)
{
	double heatConductivity = GetHeatConductivity(currentTemperature[_horizontalSplitting - 1]);
	double alpha = GetAlpha(heatConductivity, GetSpecificHeat(currentTemperature[_horizontalSplitting - 1]));
	// tempTi[0] = tempTi-1[0] + dt * aplha * (1 / (dX * dX) * 2 * (tempTi-1[1] - tempTi-1[0]) + g[0] / HC)
	currentTemperature[0] = previousTemperature[0] + _timeDifference * alpha * (1 / (_dX * _dX) * 2 *
		(previousTemperature[1] - previousTemperature[0]) + g[0] / heatConductivity);
	double heatTransferCoefficient = GetHeatTransferCoefficient(currentTemperature[_horizontalSplitting - 1], isParticle, particleIndex);
	// tempTi[n] = tempTi-1[n] + dt * alpha * (1 / (dX * dX) * 2 * (tempTi-1[n - 1] - tempTi-1[n] - dX / HC * (HTC * (tempTi-1[n] - tK))) + 1 / (horizontalSplitting * dX) * (-1 / HC) * (HTC * (tempTi-1[n] - tK)) + g[n] / HC)
	currentTemperature[_horizontalSplitting - 1] = previousTemperature[_horizontalSplitting - 1] +
		_timeDifference * alpha * (1 / (_dX * _dX) * 2 * (previousTemperature[_horizontalSplitting - 2] -
			previousTemperature[_horizontalSplitting - 1] - _dX / heatConductivity * (heatTransferCoefficient *
			(previousTemperature[_horizontalSplitting - 1] - _tK))) + 1 / (_horizontalSplitting * _dX) *
				(-1 / heatConductivity) * (heatTransferCoefficient * (previousTemperature[_horizontalSplitting - 1] -
					_tK)) + g[_horizontalSplitting - 1] / heatConductivity);

	for (int i = 1; i < _horizontalSplitting - 1; i += 1)
	{
		// tempTi[i] = tempTi-1[i] + dt * alpha * (1 / (dX * dX) * (tempTi-1[i - 1] + tempTi-1[i + 1] - 2 * tempTi-1[i]) + 1 / (i * dX) * 1 / (2 * dX) * (tempTi-1[i + 1] - tempTi-1[i - 1]) + g[i] / HC)
		currentTemperature[i] = previousTemperature[i] + _timeDifference * alpha * (1 / (_dX * _dX) *
			(previousTemperature[i - 1] + previousTemperature[i + 1] - 2 * previousTemperature[i]) +
			1 / (i * _dX) * 1 / (2 * _dX) * (previousTemperature[i + 1] - previousTemperature[i - 1]) +
			g[i] / heatConductivity);
	}

	return currentTemperature;
}

double GetHeatConductivity(double temperature)
{
	if (_heatConductivitiesSizeRows > 0)
	{
		double heatConductivity0 = 0;
		double heatConductivity1 = 0;
		double temperature0 = 0;
		double temperature1 = 0;
		int i = 0;

		while ((i < (_heatConductivitiesSizeRows * 2)) &&
			(_heatConductivities[i + 0] <= temperature))
		{
			heatConductivity0 = _heatConductivities[i + 1];
			temperature0 = _heatConductivities[i + 0];
			i += 2;
		}

		if (i < (_heatConductivitiesSizeRows * 2))
		{
			heatConductivity1 = _heatConductivities[i + 1];
			temperature1 = _heatConductivities[i + 0];

			if (i == 0)
			{
				temperature0 = 0;
				heatConductivity0 = 0;
			}

			return ((heatConductivity1 - heatConductivity0) / (temperature1 - temperature0) * temperature -
				((heatConductivity1 - heatConductivity0) / (temperature1 - temperature0) * temperature0 - heatConductivity0));
		}
		else
		{
			return _heatConductivities[(_heatConductivitiesSizeRows * 2) - 1];
		}
	}

	return 0;
}

double GetSpecificHeat(double temperature)
{
	if (_specificHeatsSizeRows > 0)
	{
		double specificHeat0 = 0;
		double specificHeat1 = 0;
		double temperature0 = 0;
		double temperature1 = 0;
		int i = 0;

		while ((i < (_specificHeatsSizeRows * 2)) &&
			(_specificHeats[i + 0] <= temperature))
		{
			specificHeat0 = _specificHeats[i + 1];
			temperature0 = _specificHeats[i + 0];
			i += 2;
		}

		if (i < (_specificHeatsSizeRows * 2))
		{
			specificHeat1 = _specificHeats[i + 1];
			temperature1 = _specificHeats[i + 0];

			if (i == 0)
			{
				temperature0 = 0;
				specificHeat0 = 0;
			}

			return ((specificHeat1 - specificHeat0) / (temperature1 - temperature0) * temperature -
				((specificHeat1 - specificHeat0) / (temperature1 - temperature0) * temperature0 - specificHeat0));
		}
		else
		{
			return _specificHeats[(_specificHeatsSizeRows * 2) - 1];
		}
	}

	return 0;
}

double GetAlpha(double heatConductivity, double specificHeat)
{
	return (heatConductivity / (specificHeat * _density));
}

double GetHeatTransferCoefficient(double temperature, bool isParticle, int particleIndex)
{
	if (!isParticle)
	{
		if (_dimensionNumber > 0)
		{
			double heatTransferCoefficient0 = 0;
			double heatTransferCoefficient1 = 0;
			double temperature0 = 0;
			double temperature1 = 0;
			int i = 0;

			while ((i < (_dimensionNumber * 2)) &&
				(_htcValues[i + 0] <= temperature))
			{
				heatTransferCoefficient0 = _htcValues[i + 1];
				temperature0 = _htcValues[i + 0];
				i += 2;
			}

			if (i < (_dimensionNumber * 2))
			{
				heatTransferCoefficient1 = _htcValues[i + 1];
				temperature1 = _htcValues[i + 0];

				if (i == 0)
				{
					temperature0 = 0;
					heatTransferCoefficient0 = 0;
				}

				return ((heatTransferCoefficient1 - heatTransferCoefficient0) / (temperature1 - temperature0) * temperature -
					((heatTransferCoefficient1 - heatTransferCoefficient0) / (temperature1 - temperature0) * temperature0 - heatTransferCoefficient0));
			}
			else
			{
				return _htcValues[(_dimensionNumber * 2) - 1];
			}
		}

		return 0;
	}
	else
	{
		if (_dimensionNumber > 0)
		{
			double heatTransferCoefficient0 = 0;
			double heatTransferCoefficient1 = 0;
			double temperature0 = 0;
			double temperature1 = 0;
			int i = 0;

			while ((i < (_dimensionNumber * 2)) &&
				(_htcValues[i + 0] <= temperature))
			{
				heatTransferCoefficient0 = _particlePosition[(particleIndex * _dimensionNumber) + (i / 2)];
				temperature0 = _htcValues[i + 0];
				i += 2;
			}

			if (i < (_dimensionNumber * 2))
			{
				heatTransferCoefficient1 = _particlePosition[(particleIndex * _dimensionNumber) + (i / 2)];
				temperature1 = _htcValues[i + 0];

				if (i == 0)
				{
					temperature0 = 0;
					heatTransferCoefficient0 = 0;
				}

				return ((heatTransferCoefficient1 - heatTransferCoefficient0) / (temperature1 - temperature0) * temperature -
					((heatTransferCoefficient1 - heatTransferCoefficient0) / (temperature1 - temperature0) * temperature0 - heatTransferCoefficient0));
			}
			else
			{
				return _particlePosition[(particleIndex * _dimensionNumber) + (_dimensionNumber - 1)];
			}
		}

		return 0;
	}
}

void WriteReferenceCooldownLogToFile()
{
	if (_referenceCooldownCurve != NULL && _referenceCooldownCurveSizeRows > 0)
	{
		std::ofstream fileStream("ReferenceCooldownLog.txt");

		if (fileStream.is_open())
		{
			fileStream << "time ";
			fileStream << "temperature ";
			fileStream << "\r\n";

			for (int i = 0; i < (_referenceCooldownCurveSizeRows * 2); i += 2)
			{
				std::stringstream ss1;
				ss1 << std::fixed << std::setprecision(3) << _referenceCooldownCurve[i + 0];
				std::string timeData = ss1.str();
				std::replace(timeData.begin(), timeData.end(), ',', '.');
				std::stringstream ss2;
				ss2 << std::fixed << std::setprecision(3) << _referenceCooldownCurve[i + 1];
				std::string temperatureData = ss2.str();
				replace(temperatureData.begin(), temperatureData.end(), ',', '.');
				fileStream << timeData << " ";
				fileStream << temperatureData << " ";
				fileStream << "\r\n";
			}

			fileStream.close();
		}
	}
}

void WritePsoGlobalBestLogToFile()
{
	if (_globalBestFitness != NULL && _globalBestPosition != NULL &&
		_globalBestSize > 0)
	{
		std::ofstream fileStream("ParticleSwarmOptimizationGlobalBestLog.txt");

		if (fileStream.is_open())
		{
			fileStream << "fitness ";

			for (int i = 0; i < _dimensionNumber; i += 1)
			{
				std::stringstream ss;
				ss << (i + 1);
				fileStream << "htc" << ss.str() << " ";
			}

			fileStream << "\r\n";

			for (int i = 0; i < _globalBestSize; i += 1)
			{
				std::stringstream ss1;
				ss1 << std::fixed << std::setprecision(3) << _globalBestFitness[i];
				std::string value = ss1.str();
				std::replace(value.begin(), value.end(), ',', '.');
				fileStream << value << " ";

				for (int j = 0; j < _dimensionNumber; j += 1)
				{
					std::stringstream ss2;
					ss2 << std::fixed << std::setprecision(3) << _globalBestPosition[(i * _dimensionNumber) + j];
					value = ss2.str();
					replace(value.begin(), value.end(), ',', '.');
					fileStream << value << " ";
				}

				fileStream << "\r\n";
			}

			fileStream.close();
		}
	}
}

void WriteExitResultAndTimeLogToFile(long elaspedMilliseconds)
{
	std::ofstream fileStream("ExitResultAndTimeLog.txt");

	if (fileStream.is_open())
	{
		std::stringstream ss;
		ss << elaspedMilliseconds;
		std::stringstream ss2;
		ss2 << _epoch;
		fileStream << "Exit reason: " << _exitReason << ", elapsed time: " << ss.str() << "ms, iterations: " << ss2.str();
		fileStream.close();
	}
}

void InitializeParticle(int particleIndex)
{
	GetInitialPosition(particleIndex);
	GetInitialVelocity(particleIndex);
	_particleFitness[particleIndex] = ObjectiveFunction(particleIndex);
	_particleBestFitness[particleIndex] = _particleFitness[particleIndex];

	for (int j = 0; j < _dimensionNumber; j += 1)
	{
		_particleBestPosition[(particleIndex * _dimensionNumber) + j] = _particlePosition[(particleIndex * _dimensionNumber) + j];
	}
}

double ObjectiveFunction(int particleIndex)
{
	int calculatedCooldownCurveSizeRows = _referenceCooldownCurveSizeRows;
	double* calculatedCooldownCurve = new double[calculatedCooldownCurveSizeRows * 2];
	double* currentTemperature = SetCurrentTemperature();
	double* previousTemperature = SetPreviousTemperature();
	double* g = SetG();

	for (int i = 0; i < (calculatedCooldownCurveSizeRows * 2); i += 2)
	{
		currentTemperature = CalculateCooldownCurve1D(true, particleIndex, currentTemperature, previousTemperature, g);
		calculatedCooldownCurve[i + 0] = ((i / 2) * _timeDifference);
		calculatedCooldownCurve[i + 1] = currentTemperature[_monitoredIndex];

		for (int j = 0; j < _horizontalSplitting; j += 1)
		{
			previousTemperature[j] = currentTemperature[j];
		}
	}

	delete[] currentTemperature;
	delete[] previousTemperature;
	delete[] g;

	if (calculatedCooldownCurve != NULL && calculatedCooldownCurveSizeRows > 0)
	{
		if (_referenceCooldownCurveSizeRows == calculatedCooldownCurveSizeRows)
		{
			double sum = 0;

			for (int i = 0; i < (_referenceCooldownCurveSizeRows * 2); i += 2)
			{
				std::stringstream ss1;
				ss1 << std::fixed << std::setprecision(1) << _referenceCooldownCurve[i + 1];
				std::string value1 = ss1.str();
				std::stringstream ss2;
				ss2 << std::fixed << std::setprecision(1) << calculatedCooldownCurve[i + 1];
				std::string value2 = ss2.str();
				sum += pow((atof(value1.c_str()) - atof(value2.c_str())), 2);
			}

			delete[] calculatedCooldownCurve;

			return sum;
		}
		else
		{
			delete[] calculatedCooldownCurve;

			return DBL_MAX;
		}
	}
	else
	{
		delete[] calculatedCooldownCurve;

		return DBL_MAX;
	}
}

void Optimize()
{
	_epoch = 0;
	int staticEpochs = 0;

	while (_epoch < _maxEpochs && staticEpochs < _maxStaticEpochs)
	{
		double beta = (0.9 - 0.55) * (_maxEpochs - _epoch) / _maxEpochs + 0.55;
		bool isErrorImproved = false;
		double mBest = GetSwarmAverageBestPosition();
		int threadCount = _particleSwarmSize;
		std::vector<std::thread> taskArray;

		for (int i = 0; i < _particleSwarmSize; i += 1)
		{
			taskArray.push_back(std::thread(OptimizePosition, i, mBest, beta));
		}

		for (int i = 0; i < _particleSwarmSize; i += 1)
		{
			taskArray[i].join();
		}

		double improovedFitness = 0;
		double* improovedPosition = new double[_dimensionNumber];

		for (int i = 0; i < _particleSwarmSize; i += 1)
		{
			if (_particleFitness[i] < _globalBestFitness[_globalBestSize - 1])
			{
				improovedFitness = _particleBestFitness[i];

				for (int j = 0; j < _dimensionNumber; j += 1)
				{
					improovedPosition[j] = _particleBestPosition[(i * _dimensionNumber) + j];
				}

				isErrorImproved = true;
				staticEpochs = 0;
			}
		}

		if (!isErrorImproved)
		{
			staticEpochs += 1;
		}
		else
		{
			double* tempFitness = new double[_globalBestSize];
			double* tempPosition = new double[_globalBestSize * _dimensionNumber];

			for (int i = 0; i < _globalBestSize; i += 1)
			{
				tempFitness[i] = _globalBestFitness[i];

				for (int j = 0; j < _dimensionNumber; j += 1)
				{
					tempPosition[(i * _dimensionNumber) + j] = _globalBestPosition[(i * _dimensionNumber) + j];
				}
			}

			delete[] _globalBestFitness;
			delete[] _globalBestPosition;
			_globalBestSize += 1;
			_globalBestFitness = new double[_globalBestSize];
			_globalBestPosition = new double[_globalBestSize * _dimensionNumber];

			for (int i = 0; i < _globalBestSize; i += 1)
			{
				if (i < _globalBestSize - 1)
				{
					_globalBestFitness[i] = tempFitness[i];

					for (int j = 0; j < _dimensionNumber; j += 1)
					{
						_globalBestPosition[(i * _dimensionNumber) + j] = tempPosition[(i * _dimensionNumber) + j];
					}
				}
				else
				{
					_globalBestFitness[i] = improovedFitness;

					for (int j = 0; j < _dimensionNumber; j += 1)
					{
						_globalBestPosition[(i * _dimensionNumber) + j] = improovedPosition[j];
					}
				}
			}

			delete[] tempFitness;
			delete[] tempPosition;

			if (fabs(_globalBestFitness[_globalBestSize - 1]) <= _particleEpsilon)
			{
				std::stringstream ss;
				ss << _particleEpsilon;
				_exitReason = "The particle swarm optimization reached an acceptable fitness value. The value was given at: " + ss.str();

				return;
			}
		}

		delete[] improovedPosition;
		_epoch += 1;

		if (_globalBestSize >= 5)
		{
			if ((_globalBestFitness[_globalBestSize - 1] / _globalBestFitness[_globalBestSize - 5]) < 0.003)
			{
				_exitReason = "The PSO global fitness value changed less then 0.003% over the last 5 global value refresh.";

				return;
			}
		}

		if (((double)_globalBestSize / (double)_epoch) < 0.03)
		{
			std::string value = "";

			switch (_particleOptimalisationType)
			{
				case 1:
					value = "Clerc";
					break;
				case 2:
					value = "Quantum";
					break;
			}

			_exitReason = "The convergence of the " + value + " PSO is too slow, it may not find the optimum";

			return;
		}

		if (staticEpochs >= _maxStaticEpochs)
		{
			std::stringstream ss;
			ss << _maxStaticEpochs;
			_exitReason = "Static iteration limit reached, limit was: " + ss.str();

			return;
		}

		if (_epoch >= _maxEpochs)
		{
			std::stringstream ss;
			ss << _maxEpochs;
			_exitReason = "Iteration limit reached, limit was: " + ss.str();

			return;
		}
	}
}

void UpdateRing()
{
	int* particleIndex = Shuffle(GetIntegerRange(0, _particleSwarmSize), _particleSwarmSize);

	for (int i = 0; i < _particleSwarmSize; i += 1)
	{
		if (_particleInformers == NULL)
		{
			_particleInformers = new int[_particleSwarmSize * _particleInformerNumber];
		}

		int numberOfinformers = (_particleInformerNumber / 2);
		int currentInformer = 0;

		for (int n = 1; n <= numberOfinformers; n += 1)
		{
			int p = i - n;

			if (p < 0)
			{
				p = _particleSwarmSize + p;
			}

			_particleInformers[(i * _particleInformerNumber) + currentInformer] = particleIndex[p];
			currentInformer += 1;
		}

		numberOfinformers += (_particleInformerNumber % 2);

		for (int n = 1; n <= numberOfinformers; n += 1)
		{
			int p = i + n;

			if (p >= _particleSwarmSize)
			{
				p = p - _particleSwarmSize;
			}

			_particleInformers[(i * _particleInformerNumber) + currentInformer] = particleIndex[p];
			currentInformer += 1;
		}
	}

	delete[] particleIndex;
}

void OptimizePosition(int particleIndex, double mbest, double beta)
{
	switch (_particleOptimalisationType)
	{
		case 1:
			UpdateClercPosition(particleIndex);
			break;
		case 2:
			UpdateQuantumPosition(particleIndex, mbest, beta);
			break;
	}

	UpdateErrorValues(particleIndex);
}

void UpdateClercPosition(int particleIndex)
{
	UpdateVelocity(particleIndex, GetBestLocalPositions(particleIndex));
	UpdatePosition(particleIndex);
}

void UpdateErrorValues(int particleIndex)
{
	double* checkPosition = new double[_dimensionNumber];

	for (int i = 0; i < _dimensionNumber; i += 1)
	{
		checkPosition[i] = _particlePosition[(particleIndex * _dimensionNumber) + i];
	}

	if (CheckIsInRange(checkPosition))
	{
		_particleFitness[particleIndex] = ObjectiveFunction(particleIndex);

		if (_particleFitness[particleIndex] < _particleBestFitness[particleIndex])
		{
			_particleBestFitness[particleIndex] = _particleFitness[particleIndex];

			for (int i = 0; i < _dimensionNumber; i += 1)
			{
				_particleBestPosition[(particleIndex * _dimensionNumber) + i] = _particlePosition[(particleIndex * _dimensionNumber) + i];
			}
		}
	}

	delete[] checkPosition;
}

void UpdatePosition(int particleIndex)
{
	for (int i = 0; i < _dimensionNumber; i += 1)
	{
		_particlePosition[(particleIndex * _dimensionNumber) + i] = (_particlePosition[(particleIndex * _dimensionNumber) + i] + _particleVelocity[(particleIndex * _dimensionNumber) + i]);
	}
}

void UpdateQuantumPosition(int particleIndex, double mbest, double beta)
{
	for (int i = 0; i < _dimensionNumber; i += 1)
	{
		double fi = ((double)rand() / (RAND_MAX));
		double p = fi * _particleBestPosition[(particleIndex * _dimensionNumber) + i] + (1 - fi) * _globalBestPosition[(_globalBestSize - 1 - _dimensionNumber) + i];

		if (fi > 0.5)
		{
			_particlePosition[(particleIndex * _dimensionNumber) + i] = p - beta * fabs(mbest - _particlePosition[(particleIndex * _dimensionNumber) + i]) * (-log10(fi));
		}
		else
		{
			_particlePosition[(particleIndex * _dimensionNumber) + i] = p + beta * fabs(mbest - _particlePosition[(particleIndex * _dimensionNumber) + i]) * (-log10(fi));
		}

		if (_particlePosition[(particleIndex * _dimensionNumber) + i] < _rangeMin)
		{
			_particlePosition[(particleIndex * _dimensionNumber) + i] = 2 * _rangeMin - _particlePosition[(particleIndex * _dimensionNumber) + i];
		}

		if (_particlePosition[(particleIndex * _dimensionNumber) + i] > _rangeMax)
		{
			_particlePosition[(particleIndex * _dimensionNumber) + i] = 2 * _rangeMax - _particlePosition[(particleIndex * _dimensionNumber) + i];
		}
	}
}

void UpdateVelocity(int particleIndex, double* bestLocalPosition)
{
	for (int i = 0; i < _dimensionNumber; i += 1)
	{
		// particle.velocity[i] = (w * paritcle.velocity[i]) + (constant1 * rand() * (particle.bestPosition[i] - particle.position[i])) + (constant2 * rand() * (bestPosition[i] - particle.position[i]))
		_particleVelocity[(particleIndex * _dimensionNumber) + i] = (_weight * _particleVelocity[(particleIndex * _dimensionNumber) + i]) +
			(_particleConstant1 * ((double)rand() / (RAND_MAX)) * (_particleBestPosition[(particleIndex * _dimensionNumber) + i] -
				_particlePosition[(particleIndex * _dimensionNumber) + i])) + (_particleConstant2 * ((double)rand() / (RAND_MAX)) *
				(bestLocalPosition[i] - _particlePosition[(particleIndex * _dimensionNumber) + i]));;
	}
}

bool CheckIsInRange(double* position)
{
	for (int i = 0; i < _dimensionNumber; i += 1)
	{
		if (position[i] < _rangeMin || position[i] > _rangeMax)
		{
			return false;
		}
	}

	return true;
}

double* GetBestLocalPositions(int particleIndex)
{
	double* returnValue = new double[_dimensionNumber];
	int bestIndex = particleIndex;

	for (int i = 0; i < _particleInformerNumber; i += 1)
	{
		if (_particleBestFitness[_particleInformers[(particleIndex * _particleInformerNumber) + i]] < _particleBestFitness[bestIndex])
		{
			bestIndex = _particleInformers[(particleIndex * _particleInformerNumber) + i];
		}
	}

	for (int i = 0; i < _dimensionNumber; i += 1)
	{
		returnValue[i] = _particleBestPosition[(bestIndex * _dimensionNumber) + i];
	}

	return returnValue;
}

void GetInitialPosition(int particalIndex)
{
	for (int i = 0; i < _dimensionNumber; i += 1)
	{
		// The given values are null.
		if (_particleInitialPosition == NULL)
		{
			_particlePosition[(particalIndex * _dimensionNumber) + i] = ((_rangeMax - _rangeMin) * ((double)rand() / (RAND_MAX)) + _rangeMin);
		}
		// There is only one given value.
		else if (_particleInitialPositionSize == 1)
		{
			_particlePosition[(particalIndex * _dimensionNumber) + i] = _particleInitialPosition[0];
		}
		// There are as many given value as the position dimensions.
		else if (_particleInitialPositionSize == _dimensionNumber)
		{
			_particlePosition[(particalIndex * _dimensionNumber) + i] = _particleInitialPosition[i];
		}
		// The current position can be set from the given values.
		else if (i < _particleInitialPositionSize)
		{
			_particlePosition[(particalIndex * _dimensionNumber) + i] = _particleInitialPosition[i];
		}
		// The current position can't be set from the given values.
		else
		{
			_particlePosition[(particalIndex * _dimensionNumber) + i] = ((_rangeMax - _rangeMin) * ((double)rand() / (RAND_MAX)) + _rangeMin);
		}
	}
}

void GetInitialVelocity(int particalIndex)
{
	for (int i = 0; i < _dimensionNumber; i += 1)
	{
		// The given values are null.
		if (_particleInitialVelocity == NULL)
		{
			_particleVelocity[(particalIndex * _dimensionNumber) + i] = ((_rangeMax - _rangeMin) * ((double)rand() / (RAND_MAX)) + _rangeMin);
		}
		// There is only one given value.
		else if (_particleInitialVelocitySize == 1)
		{
			_particleVelocity[(particalIndex * _dimensionNumber) + i] = _particleInitialVelocity[0];
		}
		// There are as many given value as the velocity dimensions.
		else if (_particleInitialVelocitySize == _dimensionNumber)
		{
			_particleVelocity[(particalIndex * _dimensionNumber) + i] = _particleInitialVelocity[i];
		}
		// The current velocity can be set from the given values.
		else if (i < _particleInitialVelocitySize)
		{
			_particleVelocity[(particalIndex * _dimensionNumber) + i] = _particleInitialVelocity[i];
		}
		// The current velocity can't be set from the given values.
		else
		{
			_particleVelocity[(particalIndex * _dimensionNumber) + i] = ((_rangeMax - _rangeMin) * ((double)rand() / (RAND_MAX)) + _rangeMin);
		}
	}
}

int* GetIntegerRange(int startIndex, int count)
{
	int* returnValue = new int[count];

	for (int i = 0; i < count; i += 1)
	{
		returnValue[i] = (startIndex + i);
	}

	return returnValue;
}

int* Shuffle(int* index, int indexSize)
{
	for (int i = 0; i < indexSize; i += 1)
	{
		int n = ((indexSize - i) * ((double)rand() / RAND_MAX) + i);

		if (n >= 100)
		{
			n = 99;
		}

		int temp = index[n];
		index[n] = index[i];
		index[i] = temp;
	}

	return index;
}

double GetSwarmAverageBestPosition()
{
	double sum = 0;

	for (int i = 0; i < _particleSwarmSize; i += 1)
	{
		for (int j = 0; j < _dimensionNumber; j += 1)
		{
			sum += _particleBestPosition[(i * _dimensionNumber) + j];
		}
	}

	return ((sum * 1.0) / _particleSwarmSize);
}

std::vector<std::string> explode(std::string originalString, char delimeter)
{
	std::vector<std::string> elems;
	split(originalString, delimeter, std::back_inserter(elems));

	return elems;
}

template<typename Out>
void split(std::string content, char delimeter, Out result)
{
	std::stringstream ss(content);
	std::string item = "";

	while (std::getline(ss, item, delimeter))
	{
		*(result++) = item;
	}
}

#pragma endregion

class Particle
{
public:
	double* position;
	double* velocity;
	int index;
	double fitness;
	double* bestPosition;
	double bestFitness;
};
