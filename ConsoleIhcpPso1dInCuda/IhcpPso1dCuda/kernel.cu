#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
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
#include <random>
#include "MersenneTwister.h"

#pragma region Fields

double* _initialTemperature;
int* _horizontalSplitting;
int* _verticalSplitting;
int* _monitoredIndex;
double* _timeDifference;
int* _dimensionNumber;
double* _htcValues;
double* _specificHeats;
int* _specificHeatsSizeRows;
double* _heatConductivities;
int* _heatConductivitiesSizeRows;
int* _rangeMin;
int* _rangeMax;
double* _simulationTime;
double* _height;
double* _radius;
double* _density;
double* _dX;
double* _dY;
double* _referenceCooldownCurve;
int* _referenceCooldownCurveSizeRows;
double* _tK;
double* _particleInitialVelocity;
int* _particleInitialVelocitySize;
double* _particleInitialPosition;
int* _particleInitialPositionSize;
int* _particleInformerNumber;
double* _particleConstant1;
double* _particleConstant2;
double* _particleEpsilon;
double* _particlePosition;
double* _particleVelocity;
double* _particleBestPosition;
double* _particleFitness;
double* _particleBestFitness;
double* _globalBestFitness;
double* _globalBestPosition;
int* _globalBestSize;
int* _particleOptimalisationType;
int* _particleSwarmSize;
int* _maxEpochs;
int* _maxStaticEpochs;
double* _weight;
int* _particleInformers;
std::string _exitReason;
int* _epoch;
cudaError_t _cudaStatus;
std::string _cudaError;

MersenneTwister* _random;

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
void WritePreIterationResult();
void WriteIterationResult(int iteration);
void WriteCudaErrorToFile();
void Optimize();
void UpdateRing();
int* GetIntegerRange(int startIndex, int count);
int* Shuffle(int* index, int indexSize);
double GetSwarmAverageBestPosition();
std::vector<std::string> explode(std::string originalString, char delimeter);
template<typename Out> void split(std::string s, char delim, Out result);
void InitializeParticle();
void OptimizePosition(double* beta, double* mBest);

__global__ void InitializeParticleKernel(int* dimensionNumber, double* particlePosition,
	double* particleVelocity, double* particleInitialPosition, int* particleInitialPositionSize, 
	double* particleInitialVelocity, int* particleInitialVelocitySize, int* rangeMin, int* rangeMax,
	double* randomNumbers1, double* randomNumbers2)
{
	int index = (blockIdx.x * blockDim.x) + threadIdx.x;

	// Set initial position.
	for (int i = 0; i < dimensionNumber[0]; i += 1)
	{
		// The given values are null.
		if (particleInitialPosition == NULL)
		{
			particlePosition[(index * dimensionNumber[0]) + i] = ((rangeMax[0] - rangeMin[0]) * randomNumbers1[(index * dimensionNumber[0]) + i] + rangeMin[0]);
		}
		// There is only one given value.
		else if (particleInitialPositionSize[0] == 1)
		{
			particlePosition[(index * dimensionNumber[0]) + i] = particleInitialPosition[0];
		}
		// There are as many given value as the position dimensions.
		else if (particleInitialPositionSize[0] == dimensionNumber[0])
		{
			particlePosition[(index * dimensionNumber[0]) + i] = particleInitialPosition[i];
		}
		// The current position can be set from the given values.
		else if (i < particleInitialPositionSize[0])
		{
			particlePosition[(index * dimensionNumber[0]) + i] = particleInitialPosition[i];
		}
		// The current position can't be set from the given values.
		else
		{
			particlePosition[(index * dimensionNumber[0]) + i] = ((rangeMax[0] - rangeMin[0]) * randomNumbers1[(index * dimensionNumber[0]) + i] + rangeMin[0]);
		}
	}

	// Set initial velocity.
	for (int i = 0; i < dimensionNumber[0]; i += 1)
	{
		// The given values are null.
		if (particleInitialVelocity == NULL)
		{
			particleVelocity[(index * dimensionNumber[0]) + i] = ((rangeMax[0] - rangeMin[0]) * randomNumbers2[(index * dimensionNumber[0]) + i] + rangeMin[0]);
		}
		// There is only one given value.
		else if (particleInitialVelocitySize[0] == 1)
		{
			particleVelocity[(index * dimensionNumber[0]) + i] = particleInitialVelocity[0];
		}
		// There are as many given value as the velocity dimensions.
		else if (particleInitialVelocitySize[0] == dimensionNumber[0])
		{
			particleVelocity[(index * dimensionNumber[0]) + i] = particleInitialVelocity[i];
		}
		// The current velocity can be set from the given values.
		else if (i < particleInitialVelocitySize[0])
		{
			particleVelocity[(index * dimensionNumber[0]) + i] = particleInitialVelocity[i];
		}
		// The current velocity can't be set from the given values.
		else
		{
			particleVelocity[(index * dimensionNumber[0]) + i] = ((rangeMax[0] - rangeMin[0]) * randomNumbers2[(index * dimensionNumber[0]) + i] + rangeMin[0]);
		}
	}
}

__global__
void OptimizePositionKernel(int* particleOptimalisationType, double* mbest, double* beta,
	int* dimensionNumber, int* particleInformerNumber, int* particleInformers, double* particleBestFitness,
	double* particleBestPosition, double* weight, double* particleVelocity, double* particlePosition,
	int* rangeMin, int* rangeMax, int* globalBestSize, double* globalBestPosition, double* particleConstant1,
	double* particleConstant2, double* randomNumbers1, double* randomNumbers2, double* randomNumbers3,
	double* bestLocalPosition)
{
	// Get the iteration number.
	int index = (blockIdx.x * blockDim.x) + threadIdx.x;

	// Update the particle's position based on the type of the PSO.
	switch (particleOptimalisationType[0])
	{
		// Cleck.
		case 1:
			int bestIndex = index;

			// Get the minimum fitness from the particle or it's informers.
			for (int i = 0; i < particleInformerNumber[0]; i += 1)
			{
				if (particleBestFitness[particleInformers[(index * particleInformerNumber[0]) + i]] < particleBestFitness[bestIndex])
				{
					bestIndex = particleInformers[(index * particleInformerNumber[0]) + i];
				}
			}

			// Set the potinion of the found minimum fitness particle.
			for (int i = 0; i < dimensionNumber[0]; i += 1)
			{
				bestLocalPosition[i] = particleBestPosition[(bestIndex * dimensionNumber[0]) + i];
			}

			// Update the particle's velocity.
			for (int i = 0; i < dimensionNumber[0]; i += 1)
			{
				particleVelocity[(index * dimensionNumber[0]) + i] = (weight[0] * particleVelocity[(index * dimensionNumber[0]) + i]) +
					(particleConstant1[0] * randomNumbers1[(index * dimensionNumber[0]) + i] * (particleBestPosition[(index * dimensionNumber[0]) + i] -
						particlePosition[(index * dimensionNumber[0]) + i])) + (particleConstant2[0] * randomNumbers2[(index * dimensionNumber[0]) + i] *
						(bestLocalPosition[i] - particlePosition[(index * dimensionNumber[0]) + i]));
			}

			// Update the particle's position.
			for (int i = 0; i < dimensionNumber[0]; i += 1)
			{
				particlePosition[(index * dimensionNumber[0]) + i] = (particlePosition[(index * dimensionNumber[0]) + i] + particleVelocity[(index * dimensionNumber[0]) + i]);
			}
			break;
		// Quantum.
		case 2:
			// Update the particle's position.
			for (int i = 0; i < dimensionNumber[0]; i += 1)
			{
				double fi = randomNumbers3[(index * dimensionNumber[0]) + i];
				double p = fi * particleBestPosition[(index * dimensionNumber[0]) + i] + (1 - fi) * globalBestPosition[(globalBestSize[0] - 1 - dimensionNumber[0]) + i];

				if (fi > 0.5)
				{
					particlePosition[(index * dimensionNumber[0]) + i] = p - beta[0] * fabs(mbest[0] - particlePosition[(index * dimensionNumber[0]) + i]) * (-log10(fi));
				}
				else
				{
					particlePosition[(index * dimensionNumber[0]) + i] = p + beta[0] * fabs(mbest[0] - particlePosition[(index * dimensionNumber[0]) + i]) * (-log10(fi));
				}

				if (particlePosition[(index * dimensionNumber[0]) + i] < rangeMin[0])
				{
					particlePosition[(index * dimensionNumber[0]) + i] = 2 * rangeMin[0] - particlePosition[(index * dimensionNumber[0]) + i];
				}

				if (particlePosition[(index * dimensionNumber[0]) + i] > rangeMax[0])
				{
					particlePosition[(index * dimensionNumber[0]) + i] = 2 * rangeMax[0] - particlePosition[(index * dimensionNumber[0]) + i];
				}
			}
			break;
	}
}

__global__
void ObjectiveFunctionKernel(int* dimensionNumber, double* particleBestFitness,
	double* particleBestPosition, double* particlePosition,
	int* rangeMin, int* rangeMax, double* particleFitness, double* initialTemperature, int* horizontalSplitting,
	int* referenceCooldownCurveSizeRows, double* referenceCooldownCurve, int* monitoredIndex,
	double* timeDifference, double* dX, double* tK, int* heatConductivitiesSizeRows, double* heatConductivities,
	int* specificHeatsSizeRows, double* specificHeats, double* density, double* htcValues,
	double* currentTemperature, double* previousTemperature, double* g)
{
	// Get the iteration number.
	int index = (blockIdx.x * blockDim.x) + threadIdx.x;
	bool isInRange = true;

	// Check if the position of the particle is the specified range.
	for (int i = 0; i < dimensionNumber[0]; i += 1)
	{
		if (particlePosition[(index * dimensionNumber[0]) + i] < rangeMin[0] ||
			particlePosition[(index * dimensionNumber[0]) + i] > rangeMax[0])
		{
			isInRange = false;
		}
	}

	// The position of the particle is in the specified range.
	if (isInRange)
	{
		particleFitness[index] = 0;
		
		// Simulate cooldown and create curve from the position of the particle as heat transfer coefficients.
		for (int i = 0; i < (referenceCooldownCurveSizeRows[0] * 2); i += 2)
		{
			double heatConductivity = 0;
			double specificHeat = 0;
			double heatTransferCoefficient = 0;
			double temperature = currentTemperature[(index * horizontalSplitting[0]) + (horizontalSplitting[0] - 1)];

			if (heatConductivitiesSizeRows[0] > 0)
			{
				double heatConductivity0 = 0;
				double heatConductivity1 = 0;
				double temperature0 = 0;
				double temperature1 = 0;
				int iterator = 0;

				while ((iterator < (heatConductivitiesSizeRows[0] * 2)) &&
					(heatConductivities[iterator + 0] <= temperature))
				{
					heatConductivity0 = heatConductivities[iterator + 1];
					temperature0 = heatConductivities[iterator + 0];
					iterator += 2;
				}

				if (iterator < (heatConductivitiesSizeRows[0] * 2))
				{
					heatConductivity1 = heatConductivities[iterator + 1];
					temperature1 = heatConductivities[iterator + 0];

					if (iterator == 0)
					{
						temperature0 = 0;
						heatConductivity0 = 0;
					}

					heatConductivity = ((heatConductivity1 - heatConductivity0) / (temperature1 - temperature0) * temperature -
						((heatConductivity1 - heatConductivity0) / (temperature1 - temperature0) * temperature0 - heatConductivity0));
				}
				else
				{
					heatConductivity = heatConductivities[(heatConductivitiesSizeRows[0] * 2) - 1];
				}
			}

			if (specificHeatsSizeRows[0] > 0)
			{
				double specificHeat0 = 0;
				double specificHeat1 = 0;
				double temperature0 = 0;
				double temperature1 = 0;
				int iterator = 0;

				while ((iterator < (specificHeatsSizeRows[0] * 2)) &&
					(specificHeats[iterator + 0] <= temperature))
				{
					specificHeat0 = specificHeats[iterator + 1];
					temperature0 = specificHeats[iterator + 0];
					iterator += 2;
				}

				if (iterator < (specificHeatsSizeRows[0] * 2))
				{
					specificHeat1 = specificHeats[iterator + 1];
					temperature1 = specificHeats[iterator + 0];

					if (iterator == 0)
					{
						temperature0 = 0;
						specificHeat0 = 0;
					}

					specificHeat = ((specificHeat1 - specificHeat0) / (temperature1 - temperature0) * temperature -
						((specificHeat1 - specificHeat0) / (temperature1 - temperature0) * temperature0 - specificHeat0));
				}
				else
				{
					specificHeat = specificHeats[(specificHeatsSizeRows[0] * 2) - 1];
				}
			}

			double alpha = (heatConductivity / (specificHeat * density[0]));
			currentTemperature[(index * horizontalSplitting[0]) + 0] = previousTemperature[(index * horizontalSplitting[0]) + 0] + timeDifference[0] * alpha * (1 / (dX[0] * dX[0]) * 2 *
				(previousTemperature[(index * horizontalSplitting[0]) + 1] - previousTemperature[(index * horizontalSplitting[0]) + 0]) + g[(index * horizontalSplitting[0]) + 0] / heatConductivity);

			if (dimensionNumber[0] > 0)
			{
				double heatTransferCoefficient0 = 0;
				double heatTransferCoefficient1 = 0;
				double temperature0 = 0;
				double temperature1 = 0;
				int iterator = 0;

				while ((iterator < (dimensionNumber[0] * 2)) &&
					(htcValues[iterator + 0] <= temperature))
				{
					heatTransferCoefficient0 = particlePosition[(index * dimensionNumber[0]) + (iterator / 2)];
					temperature0 = htcValues[iterator + 0];
					iterator += 2;
				}

				if (iterator < (dimensionNumber[0] * 2))
				{
					heatTransferCoefficient1 = particlePosition[(index * dimensionNumber[0]) + (iterator / 2)];
					temperature1 = htcValues[iterator + 0];

					if (iterator == 0)
					{
						temperature0 = 0;
						heatTransferCoefficient0 = 0;
					}

					heatTransferCoefficient = ((heatTransferCoefficient1 - heatTransferCoefficient0) / (temperature1 - temperature0) * temperature -
						((heatTransferCoefficient1 - heatTransferCoefficient0) / (temperature1 - temperature0) * temperature0 - heatTransferCoefficient0));
				}
				else
				{
					heatTransferCoefficient = particlePosition[(index * dimensionNumber[0]) + (dimensionNumber[0] - 1)];
				}
			}

			currentTemperature[(index * horizontalSplitting[0]) + (horizontalSplitting[0] - 1)] = previousTemperature[(index * horizontalSplitting[0]) + (horizontalSplitting[0] - 1)] +
				timeDifference[0] * alpha * (1 / (dX[0] * dX[0]) * 2 * (previousTemperature[(index * horizontalSplitting[0]) + (horizontalSplitting[0] - 2)] -
					previousTemperature[(index * horizontalSplitting[0]) + (horizontalSplitting[0] - 1)] - dX[0] / heatConductivity * (heatTransferCoefficient *
					(previousTemperature[(index * horizontalSplitting[0]) + (horizontalSplitting[0] - 1)] - tK[0]))) + 1 / (horizontalSplitting[0] * dX[0]) *
						(-1 / heatConductivity) * (heatTransferCoefficient * (previousTemperature[(index * horizontalSplitting[0]) + (horizontalSplitting[0] - 1)] -
							tK[0])) + g[(index * horizontalSplitting[0]) + (horizontalSplitting[0] - 1)] / heatConductivity);

			for (int j = 1; j < (horizontalSplitting[0] - 1); j += 1)
			{
				currentTemperature[(index * horizontalSplitting[0]) + j] = previousTemperature[(index * horizontalSplitting[0]) + j] + timeDifference[0] * alpha * (1 / (dX[0] * dX[0]) *
					(previousTemperature[(index * horizontalSplitting[0]) + (j - 1)] + previousTemperature[(index * horizontalSplitting[0]) + (j + 1)] - 2 *
						previousTemperature[(index * horizontalSplitting[0]) + j]) + 1 / (j * dX[0]) * 1 / (2 * dX[0]) * (previousTemperature[(index * horizontalSplitting[0]) + (j + 1)] -
							previousTemperature[(index * horizontalSplitting[0]) + (j - 1)]) + g[(index * horizontalSplitting[0]) + j] / heatConductivity);
			}

			particleFitness[index] += (currentTemperature[(index * horizontalSplitting[0]) + monitoredIndex[0]] - referenceCooldownCurve[i + 1]) * (currentTemperature[(index * horizontalSplitting[0]) + monitoredIndex[0]] - referenceCooldownCurve[i + 1]);

			// Set the current temperature values as the previous temperature values for the next iteration.
			for (int j = 0; j < horizontalSplitting[0]; j += 1)
			{
				previousTemperature[(index * horizontalSplitting[0]) + j] = currentTemperature[(index * horizontalSplitting[0]) + j];
			}
		}

		if (particleFitness[index] < particleBestFitness[index])
		{
			particleBestFitness[index] = particleFitness[index];

			for (int j = 0; j < dimensionNumber[0]; j += 1)
			{
				particleBestPosition[(index * dimensionNumber[0]) + j] = particlePosition[(index * dimensionNumber[0]) + j];
			}
		}
	}
}

int main()
{
	cudaSetDevice(0);
	ReadDataFromFile("ConfigurationIn.txt");
	CalculateReferenceCooldownCurve();
	WriteReferenceCooldownLogToFile();
	// Set particle swarm initial values.
	InitializeParticle();

	if (_cudaStatus != cudaSuccess)
	{
		WriteCudaErrorToFile();
		cudaDeviceReset();

		return -1;
	}

	//WritePreIterationResult();
	// Set informers.
	UpdateRing();

	// Set initial best values.
	for (int i = 0; i < _particleSwarmSize[0]; i += 1)
	{
		if (_particleBestFitness[i] < _globalBestFitness[_globalBestSize[0] - 1])
		{
			_globalBestFitness[_globalBestSize[0] - 1] = _particleBestFitness[i];

			for (int j = 0; j < _dimensionNumber[0]; j += 1)
			{
				_globalBestPosition[(_globalBestSize[0] - 1) + j] = _particlePosition[(i * _dimensionNumber[0]) + j];
			}
		}
	}

	clock_t start = clock();
	Optimize();
	clock_t finish = clock();

	if (_cudaStatus != cudaSuccess)
	{
		WriteCudaErrorToFile();
		cudaDeviceReset();

		return -1;
	}

	WritePsoGlobalBestLogToFile();
	WriteExitResultAndTimeLogToFile((long)(finish - start));
	cudaDeviceReset();

	return 0;
}

void InitializeParticle()
{
	#pragma region InitializeParticle

	//cudaDeviceReset();
	//cudaSetDevice(0);
	double* _randomNumbers1 = new double[_dimensionNumber[0] * _particleSwarmSize[0]];
	double* _randomNumbers2 = new double[_dimensionNumber[0] * _particleSwarmSize[0]];

	for (int i = 0; i < (_dimensionNumber[0] * _particleSwarmSize[0]); i += 1)
	{
		_randomNumbers1[i] = _random->rnd2();//((double)rand() / RAND_MAX);
		_randomNumbers2[i] = _random->rnd2();//((double)rand() / RAND_MAX);
	}

	// Create device variables.
	int* dimensionNumber1;
	double* particlePosition1;
	double* particleVelocity1;
	double* particleInitialPosition;
	int* particleInitialPositionSize;
	double* particleInitialVelocity;
	int* particleInitialVelocitySize;
	int* rangeMin;
	int* rangeMax;
	double* randomNumbers1;
	double* randomNumbers2;

	// Allocate device memory for variables.
	cudaMalloc((void**)&dimensionNumber1, sizeof(int));
	cudaMalloc((void**)&particlePosition1, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double));
	cudaMalloc((void**)&particleVelocity1, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double));
	cudaMalloc((void**)&particleInitialPosition, _particleInitialPositionSize[0] * sizeof(double));
	cudaMalloc((void**)&particleInitialPositionSize, sizeof(int));
	cudaMalloc((void**)&particleInitialVelocity, _particleInitialVelocitySize[0] * sizeof(double));
	cudaMalloc((void**)&particleInitialVelocitySize, sizeof(int));
	cudaMalloc((void**)&rangeMin, sizeof(int));
	cudaMalloc((void**)&rangeMax, sizeof(int));
	cudaMalloc((void**)&randomNumbers1, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double));
	cudaMalloc((void**)&randomNumbers2, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double));

	// Copy host variable values to device variables.
	cudaMemcpy(dimensionNumber1, _dimensionNumber, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(particlePosition1, _particlePosition, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(particleVelocity1, _particleVelocity, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(particleInitialPosition, _particleInitialPosition, _particleInitialPositionSize[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(particleInitialPositionSize, _particleInitialPositionSize, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(particleInitialVelocity, _particleInitialVelocity, _particleInitialVelocitySize[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(particleInitialVelocitySize, _particleInitialVelocitySize, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(rangeMin, _rangeMin, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(rangeMax, _rangeMax, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(randomNumbers1, _randomNumbers1, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(randomNumbers2, _randomNumbers2, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyHostToDevice);

	// Free memory.
	delete[] _randomNumbers1;
	delete[] _randomNumbers2;

	// Call kernel function.
	InitializeParticleKernel<<<1, _particleSwarmSize[0]>>>(dimensionNumber1, 
		particlePosition1, particleVelocity1, particleInitialPosition, 
		particleInitialPositionSize, particleInitialVelocity, 
		particleInitialVelocitySize, rangeMin, rangeMax, randomNumbers1, 
		randomNumbers2);

	_cudaStatus = cudaGetLastError();
	if (_cudaStatus != cudaSuccess) {
		std::stringstream ss;
		ss << cudaGetErrorString(_cudaStatus);
		_cudaError += "InitializeParticle launch failed: " + ss.str() + "\n";
		fprintf(stderr, "InitializeParticle launch failed: %s\n", cudaGetErrorString(_cudaStatus));
	}
	// Wait for all threads to finish.
	_cudaStatus = cudaDeviceSynchronize();
	if (_cudaStatus != cudaSuccess) {
		std::stringstream ss;
		ss << _cudaStatus;
		_cudaError += "cudaDeviceSynchronize returned error code " + ss.str() + " after launching InitializeParticleKernel!\n";
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching InitializeParticleKernel!\n", _cudaStatus);
	}

	// Copy device variable values to host variables.
	cudaMemcpy(_particlePosition, particlePosition1, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(_particleVelocity, particleVelocity1, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyDeviceToHost);

	// Free device variables.
	cudaFree(dimensionNumber1);
	cudaFree(particlePosition1);
	cudaFree(particleVelocity1);
	cudaFree(particleInitialPosition);
	cudaFree(particleInitialPositionSize);
	cudaFree(particleInitialVelocity);
	cudaFree(particleInitialVelocitySize);
	cudaFree(rangeMin);
	cudaFree(rangeMax);
	cudaFree(randomNumbers1);
	cudaFree(randomNumbers2);
	//cudaDeviceReset();

	#pragma endregion
	
	#pragma region ObjectiveFunction

	//cudaDeviceReset();
	//cudaSetDevice(0);
	double* _currentTemperature = new double[_horizontalSplitting[0] * _particleSwarmSize[0]];
	double* _previousTemperature = new double[_horizontalSplitting[0] * _particleSwarmSize[0]];
	double* _g = new double[_horizontalSplitting[0] * _particleSwarmSize[0]];

	for (int i = 0; i < (_horizontalSplitting[0] * _particleSwarmSize[0]); i += 1)
	{
		_currentTemperature[i] = _initialTemperature[0];
		_previousTemperature[i] = _initialTemperature[0];
		_g[i] = 0;
	}

	// Create device variables.
	int* dimensionNumber2;
	double* particleBestFitness2;
	double* particleBestPosition2;
	double* particlePosition2;
	int* rangeMin2;
	int* rangeMax2;
	double* particleFitness;
	double* initialTemperature;
	int* horizontalSplitting;
	int* referenceCooldownCurveSizeRows;
	double* referenceCooldownCurve;
	int* monitoredIndex;
	double* timeDifference;
	double* dX;
	double* tK;
	int* heatConductivitiesSizeRows;
	double* heatConductivities;
	int* specificHeatsSizeRows;
	double* specificHeats;
	double* density;
	double* htcValues;
	double* currentTemperature;
	double* previousTemperature;
	double* g;

	// Allocate device memory for variables.
	cudaMalloc((void**)&dimensionNumber2, sizeof(int));
	cudaMalloc((void**)&particleBestFitness2, _particleSwarmSize[0] * sizeof(double));
	cudaMalloc((void**)&particleBestPosition2, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double));
	cudaMalloc((void**)&particlePosition2, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double));
	cudaMalloc((void**)&rangeMin2, sizeof(int));
	cudaMalloc((void**)&rangeMax2, sizeof(int));
	cudaMalloc((void**)&particleFitness, _particleSwarmSize[0] * sizeof(double));
	cudaMalloc((void**)&initialTemperature, sizeof(double));
	cudaMalloc((void**)&horizontalSplitting, sizeof(int));
	cudaMalloc((void**)&referenceCooldownCurveSizeRows, sizeof(int));
	cudaMalloc((void**)&referenceCooldownCurve, _referenceCooldownCurveSizeRows[0] * 2 * sizeof(double));
	cudaMalloc((void**)&monitoredIndex, sizeof(int));
	cudaMalloc((void**)&timeDifference, sizeof(double));
	cudaMalloc((void**)&dX, sizeof(double));
	cudaMalloc((void**)&tK, sizeof(double));
	cudaMalloc((void**)&heatConductivitiesSizeRows, sizeof(int));
	cudaMalloc((void**)&heatConductivities, _heatConductivitiesSizeRows[0] * 2 * sizeof(double));
	cudaMalloc((void**)&specificHeatsSizeRows, sizeof(int));
	cudaMalloc((void**)&specificHeats, _specificHeatsSizeRows[0] * 2 * sizeof(double));
	cudaMalloc((void**)&density, sizeof(double));
	cudaMalloc((void**)&htcValues, _dimensionNumber[0] * 2 * sizeof(double));
	cudaMalloc((void**)&currentTemperature, _horizontalSplitting[0] * _particleSwarmSize[0] * sizeof(double));
	cudaMalloc((void**)&previousTemperature, _horizontalSplitting[0] * _particleSwarmSize[0] * sizeof(double));
	cudaMalloc((void**)&g, _horizontalSplitting[0] * _particleSwarmSize[0] * sizeof(double));

	// Copy host variable values to device variables.
	cudaMemcpy(dimensionNumber2, _dimensionNumber, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(particleBestFitness2, _particleBestFitness, _particleSwarmSize[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(particleBestPosition2, _particleBestPosition, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(particlePosition2, _particlePosition, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(rangeMin2, _rangeMin, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(rangeMax2, _rangeMax, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(particleFitness, _particleFitness, _particleSwarmSize[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(initialTemperature, _initialTemperature, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(horizontalSplitting, _horizontalSplitting, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(referenceCooldownCurveSizeRows, _referenceCooldownCurveSizeRows, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(referenceCooldownCurve, _referenceCooldownCurve, _referenceCooldownCurveSizeRows[0] * 2 * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(monitoredIndex, _monitoredIndex, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(timeDifference, _timeDifference, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dX, _dX, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(tK, _tK, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(heatConductivitiesSizeRows, _heatConductivitiesSizeRows, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(heatConductivities, _heatConductivities, _heatConductivitiesSizeRows[0] * 2 * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(specificHeatsSizeRows, _specificHeatsSizeRows, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(specificHeats, _specificHeats, _specificHeatsSizeRows[0] * 2 * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(density, _density, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(htcValues, _htcValues, _dimensionNumber[0] * 2 * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(currentTemperature, _currentTemperature, _horizontalSplitting[0] * _particleSwarmSize[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(previousTemperature, _previousTemperature, _horizontalSplitting[0] * _particleSwarmSize[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(g, _g, _horizontalSplitting[0] * _particleSwarmSize[0] * sizeof(double), cudaMemcpyHostToDevice);

	// Free memory.
	delete[] _currentTemperature;
	delete[] _previousTemperature;
	delete[] _g;

	// Call kernel function.
	ObjectiveFunctionKernel << <1, _particleSwarmSize[0] >> > (dimensionNumber2, 
		particleBestFitness2, particleBestPosition2, particlePosition2,
		rangeMin2, rangeMax2, particleFitness, initialTemperature, horizontalSplitting,
		referenceCooldownCurveSizeRows, referenceCooldownCurve, monitoredIndex,
		timeDifference, dX, tK, heatConductivitiesSizeRows, heatConductivities,
		specificHeatsSizeRows, specificHeats, density, htcValues,
		currentTemperature, previousTemperature, g);

	_cudaStatus = cudaGetLastError();
	if (_cudaStatus != cudaSuccess) {
		std::stringstream ss;
		ss << cudaGetErrorString(_cudaStatus);
		_cudaError += "ObjectiveFunctionKernel launch failed: " + ss.str() + "\n";
		fprintf(stderr, "ObjectiveFunctionKernel launch failed: %s\n", cudaGetErrorString(_cudaStatus));
	}
	// Wait for all threads to finish.
	_cudaStatus = cudaDeviceSynchronize();
	if (_cudaStatus != cudaSuccess) {
		std::stringstream ss;
		ss << _cudaStatus;
		_cudaError += "cudaDeviceSynchronize returned error code " + ss.str() + " after launching ObjectiveFunctionKernel!\n";
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching ObjectiveFunctionKernel!\n", _cudaStatus);
	}

	// Copy device variable values to host variables.
	cudaMemcpy(_particleBestFitness, particleBestFitness2, _particleSwarmSize[0] * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(_particleBestPosition, particleBestPosition2, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(_particlePosition, particlePosition2, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(_particleFitness, particleFitness, _particleSwarmSize[0] * sizeof(double), cudaMemcpyDeviceToHost);

	// Free device variables.
	cudaFree(dimensionNumber2);
	cudaFree(particleBestFitness2);
	cudaFree(particleBestPosition2);
	cudaFree(particlePosition2);
	cudaFree(rangeMin2);
	cudaFree(rangeMax2);
	cudaFree(particleFitness);
	cudaFree(initialTemperature);
	cudaFree(horizontalSplitting);
	cudaFree(referenceCooldownCurveSizeRows);
	cudaFree(referenceCooldownCurve);
	cudaFree(monitoredIndex);
	cudaFree(timeDifference);
	cudaFree(dX);
	cudaFree(tK);
	cudaFree(heatConductivitiesSizeRows);
	cudaFree(heatConductivities);
	cudaFree(specificHeatsSizeRows);
	cudaFree(specificHeats);
	cudaFree(density);
	cudaFree(htcValues);
	cudaFree(currentTemperature);
	cudaFree(previousTemperature);
	cudaFree(g);
	//cudaDeviceReset();

	for (int i = 0; i < _particleSwarmSize[0]; i += 1)
	{
		_particleBestFitness[i] = _particleFitness[i];

		for (int j = 0; j < _dimensionNumber[0]; j += 1)
		{
			_particleBestPosition[(i * _dimensionNumber[0]) + j] = _particlePosition[(i * _dimensionNumber[0]) + j];
		}
	}

	#pragma endregion
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
					_initialTemperature = new double[1];
					_initialTemperature[0] = atof(currentLine.c_str());
				}
				else if (currentLine.find("horizontal splitting") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_horizontalSplitting = new int[1];
					_horizontalSplitting[0] = atoi(currentLine.c_str());
				}
				else if (currentLine.find("vertical splitting") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_verticalSplitting = new int[1];
					_verticalSplitting[0] = atoi(currentLine.c_str());
				}
				else if (currentLine.find("monitored index") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_monitoredIndex = new int[1];
					_monitoredIndex[0] = atoi(currentLine.c_str());
				}
				else if (currentLine.find("range min") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_rangeMin = new int[1];
					_rangeMin[0] = atoi(currentLine.c_str());
				}
				else if (currentLine.find("range max") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_rangeMax = new int[1];
					_rangeMax[0] = atoi(currentLine.c_str());
				}
				else if (currentLine.find("simulation time") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_simulationTime = new double[1];
					_simulationTime[0] = atof(currentLine.c_str());
				}
				else if (currentLine.find("height") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_height = new double[1];
					_height[0] = atof(currentLine.c_str());
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
					_timeDifference = new double[1];
					_timeDifference[0] = atof(currentLine.c_str());
				}
				else if (currentLine.find("pso epsilon") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_particleEpsilon = new double[1];
					_particleEpsilon[0] = atof(currentLine.c_str());
				}
				else if (currentLine.find("radius") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_radius = new double[1];
					_radius[0] = atof(currentLine.c_str());
				}
				else if (currentLine.find("density") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_density = new double[1];
					_density[0] = atof(currentLine.c_str());
				}
				else if (currentLine.find("optimization type") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					_particleOptimalisationType = new int[1];

					if (currentLine == "clerc")
					{
						_particleOptimalisationType[0] = 1;
					}
					else if (currentLine == "quantum")
					{
						_particleOptimalisationType[0] = 2;
					}
					else
					{
						_particleOptimalisationType[0] = 1;
					}
				}
				else if (currentLine.find("swarm size") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_particleSwarmSize = new int[1];
					_particleSwarmSize[0] = atoi(currentLine.c_str());
				}
				else if (currentLine.find("max epochs") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_maxEpochs = new int[1];
					_maxEpochs[0] = atoi(currentLine.c_str());
				}
				else if (currentLine.find("max static epochs") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_maxStaticEpochs = new int[1];
					_maxStaticEpochs[0] = atoi(currentLine.c_str());
				}
				else if (currentLine.find("weight") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_weight = new double[1];
					_weight[0] = atof(currentLine.c_str());
				}
				else if (currentLine.find("constant1") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_particleConstant1 = new double[1];
					_particleConstant1[0] = atof(currentLine.c_str());
				}
				else if (currentLine.find("constant2") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_particleConstant2 = new double[1];
					_particleConstant2[0] = atof(currentLine.c_str());
				}
				else if (currentLine.find("informer number") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_particleInformerNumber = new int[1];
					_particleInformerNumber[0] = atoi(currentLine.c_str());
				}
				else if (currentLine.find("initial position") != std::string::npos)
				{
					getline(fileStream, currentLine);
					remove(currentLine.begin(), currentLine.end(), ' ');
					replace(currentLine.begin(), currentLine.end(), ',', '.');
					_particleInitialPositionSize = new int[1];
					_particleInitialPositionSize[0] = atoi(currentLine.c_str());

					if (_particleInitialPositionSize[0] > 0)
					{
						_particleInitialPosition = new double[_particleInitialPositionSize[0]];

						for (int i = 0; i < _particleInitialPositionSize[0]; i += 1)
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
					_particleInitialVelocitySize = new int[1];
					_particleInitialVelocitySize[0] = atoi(currentLine.c_str());

					if (_particleInitialVelocitySize[0] > 0)
					{
						_particleInitialVelocity = new double[_particleInitialVelocitySize[0]];

						for (int i = 0; i < _particleInitialVelocitySize[0]; i += 1)
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
					_heatConductivitiesSizeRows = new int[1];
					_heatConductivitiesSizeRows[0] = atoi(currentLine.c_str());

					if (_heatConductivitiesSizeRows[0] > 0)
					{
						_heatConductivities = new double[_heatConductivitiesSizeRows[0] * 2];

						for (int i = 0; i < (_heatConductivitiesSizeRows[0] * 2); i += 2)
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
					_specificHeatsSizeRows = new int[1];
					_specificHeatsSizeRows[0] = atoi(currentLine.c_str());

					if (_specificHeatsSizeRows[0] > 0)
					{
						_specificHeats = new double[_specificHeatsSizeRows[0] * 2];

						for (int i = 0; i < (_specificHeatsSizeRows[0] * 2); i += 2)
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
					_dimensionNumber = new int[1];
					_dimensionNumber[0] = atoi(currentLine.c_str());

					if (_dimensionNumber[0] > 0)
					{
						_htcValues = new double[_dimensionNumber[0] * 2];

						for (int i = 0; i < (_dimensionNumber[0] * 2); i += 2)
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

	_dX = new double[1];
	_dX[0] = (_radius[0] / 1000) / _horizontalSplitting[0];
	_dY = new double[1];
	_dY[0] = (_height[0] / 1000) / _verticalSplitting[0];
	_tK = new double[1];
	_tK[0] = 0;
	_cudaStatus = cudaSuccess;
	_exitReason = "";
	// Get particle swarm initial values.
	_particlePosition = new double[_particleSwarmSize[0] * _dimensionNumber[0]];
	_particleVelocity = new double[_particleSwarmSize[0] * _dimensionNumber[0]];
	_particleBestPosition = new double[_particleSwarmSize[0] * _dimensionNumber[0]];
	_particleFitness = new double[_particleSwarmSize[0]];
	_particleBestFitness = new double[_particleSwarmSize[0]];
	_globalBestSize = new int[1];
	_globalBestSize[0] = 1;
	_globalBestFitness = new double[_globalBestSize[0]];
	_globalBestFitness[_globalBestSize[0] - 1] = DBL_MAX;
	_globalBestPosition = new double[_globalBestSize[0] * _dimensionNumber[0]];
	_epoch = new int[1];
	_epoch[0] = 0;
	_random = new MersenneTwister();
}

void CalculateReferenceCooldownCurve()
{
	_referenceCooldownCurveSizeRows = new int[1];
	_referenceCooldownCurveSizeRows[0] = (int)(_simulationTime[0] / _timeDifference[0]) + 1;
	_referenceCooldownCurve = new double[_referenceCooldownCurveSizeRows[0] * 2];
	double* currentTemperature = new double[_horizontalSplitting[0]];
	double* previousTemperature = new double[_horizontalSplitting[0]];
	double* g = new double[_horizontalSplitting[0]];

	for (int i = 0; i < _horizontalSplitting[0]; i += 1)
	{
		currentTemperature[i] = _initialTemperature[0];
		previousTemperature[i] = _initialTemperature[0];
		g[i] = 0;
	}
	
	for (int i = 0; i < (_referenceCooldownCurveSizeRows[0] * 2); i += 2)
	{
		currentTemperature = CalculateCooldownCurve1D(false, -1, currentTemperature, previousTemperature, g);
		_referenceCooldownCurve[i + 0] = ((i / 2) * _timeDifference[0]);
		_referenceCooldownCurve[i + 1] = currentTemperature[_monitoredIndex[0]];

		for (int j = 0; j < _horizontalSplitting[0]; j += 1)
		{
			previousTemperature[j] = currentTemperature[j];
		}
	}

	delete[] currentTemperature;
	delete[] previousTemperature;
	delete[] g;
}

double* CalculateCooldownCurve1D(bool isParticle, int particleIndex, 
	double* currentTemperature, double* previousTemperature, double* g)
{
	double temperature = currentTemperature[_horizontalSplitting[0] - 1];
	double heatConductivity = GetHeatConductivity(temperature);
	double alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
	currentTemperature[0] = previousTemperature[0] + _timeDifference[0] * alpha * (1 / (_dX[0] * _dX[0]) * 2 *
		(previousTemperature[1] - previousTemperature[0]) + g[0] / heatConductivity);
	double heatTransferCoefficient = GetHeatTransferCoefficient(temperature, isParticle, particleIndex);
	currentTemperature[_horizontalSplitting[0] - 1] = previousTemperature[_horizontalSplitting[0] - 1] +
		_timeDifference[0] * alpha * (1 / (_dX[0] * _dX[0]) * 2 * (previousTemperature[_horizontalSplitting[0] - 2] -
			previousTemperature[_horizontalSplitting[0] - 1] - _dX[0] / heatConductivity * (heatTransferCoefficient *
			(previousTemperature[_horizontalSplitting[0] - 1] - _tK[0]))) + 1 / (_horizontalSplitting[0] * _dX[0]) *
				(-1 / heatConductivity) * (heatTransferCoefficient * (previousTemperature[_horizontalSplitting[0] - 1] -
					_tK[0])) + g[_horizontalSplitting[0] - 1] / heatConductivity);

	for (int i = 1; i < _horizontalSplitting[0] - 1; i += 1)
	{
		currentTemperature[i] = previousTemperature[i] + _timeDifference[0] * alpha * (1 / (_dX[0] * _dX[0]) *
			(previousTemperature[i - 1] + previousTemperature[i + 1] - 2 * previousTemperature[i]) +
			1 / (i * _dX[0]) * 1 / (2 * _dX[0]) * (previousTemperature[i + 1] - previousTemperature[i - 1]) +
			g[i] / heatConductivity);
	}

	return currentTemperature;
}

double GetHeatConductivity(double temperature)
{
	if (_heatConductivitiesSizeRows[0] > 0)
	{
		double heatConductivity0 = 0;
		double heatConductivity1 = 0;
		double temperature0 = 0;
		double temperature1 = 0;
		int i = 0;

		while ((i < (_heatConductivitiesSizeRows[0] * 2)) &&
			(_heatConductivities[i + 0] <= temperature))
		{
			heatConductivity0 = _heatConductivities[i + 1];
			temperature0 = _heatConductivities[i + 0];
			i += 2;
		}

		if (i < (_heatConductivitiesSizeRows[0] * 2))
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
			return _heatConductivities[(_heatConductivitiesSizeRows[0] * 2) - 1];
		}
	}

	return 0;
}

double GetSpecificHeat(double temperature)
{
	if (_specificHeatsSizeRows[0] > 0)
	{
		double specificHeat0 = 0;
		double specificHeat1 = 0;
		double temperature0 = 0;
		double temperature1 = 0;
		int i = 0;

		while ((i < (_specificHeatsSizeRows[0] * 2)) &&
			(_specificHeats[i + 0] <= temperature))
		{
			specificHeat0 = _specificHeats[i + 1];
			temperature0 = _specificHeats[i + 0];
			i += 2;
		}

		if (i < (_specificHeatsSizeRows[0] * 2))
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
			return _specificHeats[(_specificHeatsSizeRows[0] * 2) - 1];
		}
	}

	return 0;
}

double GetAlpha(double heatConductivity, double specificHeat)
{
	return (heatConductivity / (specificHeat * _density[0]));
}

double GetHeatTransferCoefficient(double temperature, bool isParticle, int particleIndex)
{
	if (!isParticle)
	{
		if (_dimensionNumber[0] > 0)
		{
			double heatTransferCoefficient0 = 0;
			double heatTransferCoefficient1 = 0;
			double temperature0 = 0;
			double temperature1 = 0;
			int i = 0;

			while ((i < (_dimensionNumber[0] * 2)) &&
				(_htcValues[i + 0] <= temperature))
			{
				heatTransferCoefficient0 = _htcValues[i + 1];
				temperature0 = _htcValues[i + 0];
				i += 2;
			}

			if (i < (_dimensionNumber[0] * 2))
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
				return _htcValues[(_dimensionNumber[0] * 2) - 1];
			}
		}

		return 0;
	}
	else
	{
		if (_dimensionNumber[0] > 0)
		{
			double heatTransferCoefficient0 = 0;
			double heatTransferCoefficient1 = 0;
			double temperature0 = 0;
			double temperature1 = 0;
			int i = 0;

			while ((i < (_dimensionNumber[0] * 2)) &&
				(_htcValues[i + 0] <= temperature))
			{
				heatTransferCoefficient0 = _particlePosition[(particleIndex * _dimensionNumber[0]) + (i / 2)];
				temperature0 = _htcValues[i + 0];
				i += 2;
			}

			if (i < (_dimensionNumber[0] * 2))
			{
				heatTransferCoefficient1 = _particlePosition[(particleIndex * _dimensionNumber[0]) + (i / 2)];
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
				return _particlePosition[(particleIndex * _dimensionNumber[0]) + (_dimensionNumber[0] - 1)];
			}
		}

		return 0;
	}
}

void WriteReferenceCooldownLogToFile()
{
	if (_referenceCooldownCurve != NULL && _referenceCooldownCurveSizeRows[0] > 0)
	{
		std::ofstream fileStream("ReferenceCooldownLog.txt");

		if (fileStream.is_open())
		{
			fileStream << "time ";
			fileStream << "temperature ";
			fileStream << "\r\n";

			for (int i = 0; i < (_referenceCooldownCurveSizeRows[0] * 2); i += 2)
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
		_globalBestSize[0] > 0)
	{
		std::ofstream fileStream("ParticleSwarmOptimizationGlobalBestLog.txt");

		if (fileStream.is_open())
		{
			fileStream << "fitness ";

			for (int i = 0; i < _dimensionNumber[0]; i += 1)
			{
				std::stringstream ss;
				ss << (i + 1);
				fileStream << "htc" << ss.str() << " ";
			}

			fileStream << "\r\n";

			for (int i = 0; i < _globalBestSize[0]; i += 1)
			{
				std::stringstream ss1;
				ss1 << std::fixed << std::setprecision(3) << _globalBestFitness[i];
				std::string value = ss1.str();
				std::replace(value.begin(), value.end(), ',', '.');
				fileStream << value << " ";

				for (int j = 0; j < _dimensionNumber[0]; j += 1)
				{
					std::stringstream ss2;
					ss2 << std::fixed << std::setprecision(3) << _globalBestPosition[(i * _dimensionNumber[0]) + j];
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
		ss2 << _epoch[0];
		fileStream << "Exit reason: " << _exitReason << ", elapsed time: " << ss.str() << "ms, iterations: " << ss2.str();
		fileStream.close();
	}
}

void WriteCudaErrorToFile()
{
	std::ofstream fileStream("CudaError.txt");

	if (fileStream.is_open())
	{
		fileStream << _cudaError;
		fileStream.close();
	}
}

void WritePreIterationResult()
{
	std::ofstream fileStream("PreIteration.txt");

	if (fileStream.is_open())
	{
		fileStream << "fitness ";

		for (int i = 0; i < _dimensionNumber[0]; i += 1)
		{
			std::stringstream ss;
			ss << (i + 1);
			fileStream << "position" << ss.str() << " ";
		}

		for (int i = 0; i < _dimensionNumber[0]; i += 1)
		{
			std::stringstream ss;
			ss << (i + 1);
			fileStream << "velocity" << ss.str() << " ";
		}

		fileStream << "bestFitness ";

		for (int i = 0; i < _dimensionNumber[0]; i += 1)
		{
			std::stringstream ss;
			ss << (i + 1);
			fileStream << "bestPosition" << ss.str() << " ";
		}

		fileStream << "\r\n";

		for (int i = 0; i < _particleSwarmSize[0]; i += 1)
		{
			std::stringstream ss1;
			ss1 << std::fixed << std::setprecision(3) << _particleFitness[i];
			fileStream << ss1.str() << " ";

			for (int j = 0; j < _dimensionNumber[0]; j += 1)
			{
				std::stringstream ss;
				ss << std::fixed << std::setprecision(3) << _particlePosition[(i * _dimensionNumber[0]) + j];
				fileStream << ss.str() << " ";
			}

			for (int j = 0; j < _dimensionNumber[0]; j += 1)
			{
				std::stringstream ss;
				ss << std::fixed << std::setprecision(3) << _particleVelocity[(i * _dimensionNumber[0]) + j];
				fileStream << ss.str() << " ";
			}

			std::stringstream ss2;
			ss2 << std::fixed << std::setprecision(3) << _particleBestFitness[i];
			fileStream << ss2.str() << " ";

			for (int j = 0; j < _dimensionNumber[0]; j += 1)
			{
				std::stringstream ss;
				ss << std::fixed << std::setprecision(3) << _particleBestPosition[(i * _dimensionNumber[0]) + j];
				fileStream << ss.str() << " ";
			}

			fileStream << "\r\n";
		}

		fileStream.close();
	}
}

void WriteIterationResult(int iteration)
{
	std::stringstream ss0;
	ss0 << iteration;
	std::ofstream fileStream("Iteration" + ss0.str() + ".txt");

	if (fileStream.is_open())
	{
		fileStream << "fitness ";

		for (int i = 0; i < _dimensionNumber[0]; i += 1)
		{
			std::stringstream ss;
			ss << (i + 1);
			fileStream << "position" << ss.str() << " ";
		}

		for (int i = 0; i < _dimensionNumber[0]; i += 1)
		{
			std::stringstream ss;
			ss << (i + 1);
			fileStream << "velocity" << ss.str() << " ";
		}

		fileStream << "bestFitness ";

		for (int i = 0; i < _dimensionNumber[0]; i += 1)
		{
			std::stringstream ss;
			ss << (i + 1);
			fileStream << "bestPosition" << ss.str() << " ";
		}

		fileStream << "\r\n";

		for (int i = 0; i < _particleSwarmSize[0]; i += 1)
		{
			std::stringstream ss1;
			ss1 << std::fixed << std::setprecision(3) << _particleFitness[i];
			fileStream << ss1.str() << " ";

			for (int j = 0; j < _dimensionNumber[0]; j += 1)
			{
				std::stringstream ss;
				ss << std::fixed << std::setprecision(3) << _particlePosition[(i * _dimensionNumber[0]) + j];
				fileStream << ss.str() << " ";
			}

			for (int j = 0; j < _dimensionNumber[0]; j += 1)
			{
				std::stringstream ss;
				ss << std::fixed << std::setprecision(3) << _particleVelocity[(i * _dimensionNumber[0]) + j];
				fileStream << ss.str() << " ";
			}

			std::stringstream ss2;
			ss2 << std::fixed << std::setprecision(3) << _particleBestFitness[i];
			fileStream << ss2.str() << " ";

			for (int j = 0; j < _dimensionNumber[0]; j += 1)
			{
				std::stringstream ss;
				ss << std::fixed << std::setprecision(3) << _particleBestPosition[(i * _dimensionNumber[0]) + j];
				fileStream << ss.str() << " ";
			}

			fileStream << "\r\n";
		}

		fileStream.close();
	}
}

void OptimizePosition(double* beta, double* mBest)
{
	#pragma region OptimizePosition

	//cudaDeviceReset();
	//cudaSetDevice(0);
	double* _randomNumbers1 = new double[_particleSwarmSize[0] * _dimensionNumber[0]];
	double* _randomNumbers2 = new double[_particleSwarmSize[0] * _dimensionNumber[0]];
	double* _randomNumbers3 = new double[_particleSwarmSize[0] * _dimensionNumber[0]];
	double* _bestLocalPosition = new double[_dimensionNumber[0]];

	for (int i = 0; i < (_particleSwarmSize[0] * _dimensionNumber[0]); i += 1)
	{
		_randomNumbers1[i] = _random->rnd2();//((double)rand() / RAND_MAX);
		_randomNumbers2[i] = _random->rnd2();//((double)rand() / RAND_MAX);
		_randomNumbers3[i] = _random->rnd2();//((double)rand() / RAND_MAX);
	}

	// Create device variables.
	int* particleOptimalisationType;
	double* dMBest;
	double* dBeta;
	int* dimensionNumber1;
	int* particleInformerNumber;
	int* particleInformers;
	double* particleBestFitness1;
	double* particleBestPosition1;
	double* weight;
	double* particleVelocity;
	double* particlePosition1;
	int* rangeMin1;
	int* rangeMax1;
	int* globalBestSize;
	double* globalBestPosition;
	double* particleConstant1;
	double* particleConstant2;
	double* randomNumbers1;
	double* randomNumbers2;
	double* randomNumbers3;
	double* bestLocalPosition;

	// Allocate device memory for variables.
	cudaMalloc((void**)&particleOptimalisationType, sizeof(int));
	cudaMalloc((void**)&dMBest, sizeof(double));
	cudaMalloc((void**)&dBeta, sizeof(double));
	cudaMalloc((void**)&dimensionNumber1, sizeof(int));
	cudaMalloc((void**)&particleInformerNumber, sizeof(int));
	cudaMalloc((void**)&particleInformers, _particleInformerNumber[0] * _particleSwarmSize[0] * sizeof(int));
	cudaMalloc((void**)&particleBestFitness1, _particleSwarmSize[0] * sizeof(double));
	cudaMalloc((void**)&particleBestPosition1, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double));
	cudaMalloc((void**)&weight, sizeof(double));
	cudaMalloc((void**)&particleVelocity, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double));
	cudaMalloc((void**)&particlePosition1, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double));
	cudaMalloc((void**)&rangeMin1, sizeof(int));
	cudaMalloc((void**)&rangeMax1, sizeof(int));
	cudaMalloc((void**)&globalBestSize, sizeof(int));
	cudaMalloc((void**)&globalBestPosition, _globalBestSize[0] * _dimensionNumber[0] * sizeof(double));
	cudaMalloc((void**)&particleConstant1, sizeof(double));
	cudaMalloc((void**)&particleConstant2, sizeof(double));
	cudaMalloc((void**)&randomNumbers1, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double));
	cudaMalloc((void**)&randomNumbers2, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double));
	cudaMalloc((void**)&randomNumbers3, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double));
	cudaMalloc((void**)&bestLocalPosition, _dimensionNumber[0] * sizeof(double));

	// Copy host variable values to device variables.
	cudaMemcpy(particleOptimalisationType, _particleOptimalisationType, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dMBest, mBest, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dBeta, beta, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dimensionNumber1, _dimensionNumber, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(particleInformerNumber, _particleInformerNumber, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(particleInformers, _particleInformers, _particleInformerNumber[0] * _particleSwarmSize[0] * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(particleBestFitness1, _particleBestFitness, _particleSwarmSize[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(particleBestPosition1, _particleBestPosition, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(weight, _weight, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(particleVelocity, _particleVelocity, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(particlePosition1, _particlePosition, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(rangeMin1, _rangeMin, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(rangeMax1, _rangeMax, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(globalBestSize, _globalBestSize, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(globalBestPosition, _globalBestPosition, _globalBestSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(particleConstant1, _particleConstant1, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(particleConstant2, _particleConstant2, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(randomNumbers1, _randomNumbers1, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(randomNumbers2, _randomNumbers2, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(randomNumbers3, _randomNumbers3, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(bestLocalPosition, _bestLocalPosition, _dimensionNumber[0] * sizeof(double), cudaMemcpyHostToDevice);

	// Free memory.
	delete[] _randomNumbers1;
	delete[] _randomNumbers2;
	delete[] _randomNumbers3;
	delete[] _bestLocalPosition;

	// Call kernel function.
	OptimizePositionKernel<<<1, _particleSwarmSize[0]>>>(particleOptimalisationType, 
		dMBest, dBeta, dimensionNumber1, particleInformerNumber, particleInformers, 
		particleBestFitness1, particleBestPosition1, weight, particleVelocity, 
		particlePosition1, rangeMin1, rangeMax1, globalBestSize, globalBestPosition, 
		particleConstant1, particleConstant2, randomNumbers1, randomNumbers2, 
		randomNumbers3, bestLocalPosition);

	_cudaStatus = cudaGetLastError();
	if (_cudaStatus != cudaSuccess) {
		std::stringstream ss;
		ss << cudaGetErrorString(_cudaStatus);
		_cudaError += "OptimizePosition launch failed: " + ss.str() + "\n";
		fprintf(stderr, "OptimizePosition launch failed: %s\n", cudaGetErrorString(_cudaStatus));
	}
	// Wait for all threads to finish.
	cudaDeviceSynchronize();
	if (_cudaStatus != cudaSuccess) {
		std::stringstream ss;
		ss << _cudaStatus;
		_cudaError += "cudaDeviceSynchronize returned error code " + ss.str() + " after launching OptimizePositionKernel!\n";
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching OptimizePositionKernel!\n", _cudaStatus);
	}

	// Copy device variable values to host variables.
	cudaMemcpy(_particleBestFitness, particleBestFitness1, _particleSwarmSize[0] * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(_particleBestPosition, particleBestPosition1, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(_particleVelocity, particleVelocity, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(_particlePosition, particlePosition1, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyDeviceToHost);

	// Free device variables.
	cudaFree(particleOptimalisationType);
	cudaFree(dMBest);
	cudaFree(dBeta);
	cudaFree(dimensionNumber1);
	cudaFree(particleInformerNumber);
	cudaFree(particleInformers);
	cudaFree(particleBestFitness1);
	cudaFree(particleBestPosition1);
	cudaFree(weight);
	cudaFree(particleVelocity);
	cudaFree(particlePosition1);
	cudaFree(rangeMin1);
	cudaFree(rangeMax1);
	cudaFree(globalBestSize);
	cudaFree(globalBestPosition);
	cudaFree(particleConstant1);
	cudaFree(particleConstant2);
	cudaFree(randomNumbers1);
	cudaFree(randomNumbers2);
	cudaFree(randomNumbers3);
	cudaFree(bestLocalPosition);
	//cudaDeviceReset();

	#pragma endregion

	#pragma region ObjectiveFunction

	//cudaDeviceReset();
	//cudaSetDevice(0);
	double* _currentTemperature = new double[_horizontalSplitting[0] * _particleSwarmSize[0]];
	double* _previousTemperature = new double[_horizontalSplitting[0] * _particleSwarmSize[0]];
	double* _g = new double[_horizontalSplitting[0] * _particleSwarmSize[0]];

	for (int i = 0; i < (_horizontalSplitting[0] * _particleSwarmSize[0]); i += 1)
	{
		_currentTemperature[i] = _initialTemperature[0];
		_previousTemperature[i] = _initialTemperature[0];
		_g[i] = 0;
	}

	// Create device variables.
	int* dimensionNumber2;
	double* particleBestFitness2;
	double* particleBestPosition2;
	double* particlePosition2;
	int* rangeMin2;
	int* rangeMax2;
	double* particleFitness;
	double* initialTemperature;
	int* horizontalSplitting;
	int* referenceCooldownCurveSizeRows;
	double* referenceCooldownCurve;
	int* monitoredIndex;
	double* timeDifference;
	double* dX;
	double* tK;
	int* heatConductivitiesSizeRows;
	double* heatConductivities;
	int* specificHeatsSizeRows;
	double* specificHeats;
	double* density;
	double* htcValues;
	double* currentTemperature;
	double* previousTemperature;
	double* g;

	// Allocate device memory for variables.
	cudaMalloc((void**)&dimensionNumber2, sizeof(int));
	cudaMalloc((void**)&particleBestFitness2, _particleSwarmSize[0] * sizeof(double));
	cudaMalloc((void**)&particleBestPosition2, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double));
	cudaMalloc((void**)&particlePosition2, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double));
	cudaMalloc((void**)&rangeMin2, sizeof(int));
	cudaMalloc((void**)&rangeMax2, sizeof(int));
	cudaMalloc((void**)&particleFitness, _particleSwarmSize[0] * sizeof(double));
	cudaMalloc((void**)&initialTemperature, sizeof(double));
	cudaMalloc((void**)&horizontalSplitting, sizeof(int));
	cudaMalloc((void**)&referenceCooldownCurveSizeRows, sizeof(int));
	cudaMalloc((void**)&referenceCooldownCurve, _referenceCooldownCurveSizeRows[0] * 2 * sizeof(double));
	cudaMalloc((void**)&monitoredIndex, sizeof(int));
	cudaMalloc((void**)&timeDifference, sizeof(double));
	cudaMalloc((void**)&dX, sizeof(double));
	cudaMalloc((void**)&tK, sizeof(double));
	cudaMalloc((void**)&heatConductivitiesSizeRows, sizeof(int));
	cudaMalloc((void**)&heatConductivities, _heatConductivitiesSizeRows[0] * 2 * sizeof(double));
	cudaMalloc((void**)&specificHeatsSizeRows, sizeof(int));
	cudaMalloc((void**)&specificHeats, _specificHeatsSizeRows[0] * 2 * sizeof(double));
	cudaMalloc((void**)&density, sizeof(double));
	cudaMalloc((void**)&htcValues, _dimensionNumber[0] * 2 * sizeof(double));
	cudaMalloc((void**)&currentTemperature, _horizontalSplitting[0] * _particleSwarmSize[0] * sizeof(double));
	cudaMalloc((void**)&previousTemperature, _horizontalSplitting[0] * _particleSwarmSize[0] * sizeof(double));
	cudaMalloc((void**)&g, _horizontalSplitting[0] * _particleSwarmSize[0] * sizeof(double));

	// Copy host variable values to device variables.
	cudaMemcpy(dimensionNumber2, _dimensionNumber, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(particleBestFitness2, _particleBestFitness, _particleSwarmSize[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(particleBestPosition2, _particleBestPosition, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(particlePosition2, _particlePosition, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(rangeMin2, _rangeMin, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(rangeMax2, _rangeMax, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(particleFitness, _particleFitness, _particleSwarmSize[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(initialTemperature, _initialTemperature, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(horizontalSplitting, _horizontalSplitting, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(referenceCooldownCurveSizeRows, _referenceCooldownCurveSizeRows, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(referenceCooldownCurve, _referenceCooldownCurve, _referenceCooldownCurveSizeRows[0] * 2 * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(monitoredIndex, _monitoredIndex, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(timeDifference, _timeDifference, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dX, _dX, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(tK, _tK, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(heatConductivitiesSizeRows, _heatConductivitiesSizeRows, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(heatConductivities, _heatConductivities, _heatConductivitiesSizeRows[0] * 2 * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(specificHeatsSizeRows, _specificHeatsSizeRows, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(specificHeats, _specificHeats, _specificHeatsSizeRows[0] * 2 * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(density, _density, sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(htcValues, _htcValues, _dimensionNumber[0] * 2 * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(currentTemperature, _currentTemperature, _horizontalSplitting[0] * _particleSwarmSize[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(previousTemperature, _previousTemperature, _horizontalSplitting[0] * _particleSwarmSize[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(g, _g, _horizontalSplitting[0] * _particleSwarmSize[0] * sizeof(double), cudaMemcpyHostToDevice);

	// Free memory.
	delete[] _currentTemperature;
	delete[] _previousTemperature;
	delete[] _g;

	// Call kernel function.
	ObjectiveFunctionKernel << <1, _particleSwarmSize[0] >> > (dimensionNumber2, 
		particleBestFitness2, particleBestPosition2, particlePosition2,
		rangeMin2, rangeMax2, particleFitness, initialTemperature, horizontalSplitting,
		referenceCooldownCurveSizeRows, referenceCooldownCurve, monitoredIndex,
		timeDifference, dX, tK, heatConductivitiesSizeRows, heatConductivities,
		specificHeatsSizeRows, specificHeats, density, htcValues,
		currentTemperature, previousTemperature, g);

	_cudaStatus = cudaGetLastError();
	if (_cudaStatus != cudaSuccess) {
		std::stringstream ss;
		ss << cudaGetErrorString(_cudaStatus);
		_cudaError += "ObjectiveFunctionKernel launch failed: " + ss.str() + "\n";
		fprintf(stderr, "ObjectiveFunctionKernel launch failed: %s\n", cudaGetErrorString(_cudaStatus));
	}
	// Wait for all threads to finish.
	_cudaStatus = cudaDeviceSynchronize();
	if (_cudaStatus != cudaSuccess) {
		std::stringstream ss;
		ss << _cudaStatus;
		_cudaError += "cudaDeviceSynchronize returned error code " + ss.str() + " after launching ObjectiveFunctionKernel!\n";
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching ObjectiveFunctionKernel!\n", _cudaStatus);
	}

	// Copy device variable values to host variables.
	cudaMemcpy(_particleBestFitness, particleBestFitness2, _particleSwarmSize[0] * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(_particleBestPosition, particleBestPosition2, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(_particlePosition, particlePosition2, _particleSwarmSize[0] * _dimensionNumber[0] * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(_particleFitness, particleFitness, _particleSwarmSize[0] * sizeof(double), cudaMemcpyDeviceToHost);

	// Free device variables.
	cudaFree(dimensionNumber2);
	cudaFree(particleBestFitness2);
	cudaFree(particleBestPosition2);
	cudaFree(particlePosition2);
	cudaFree(rangeMin2);
	cudaFree(rangeMax2);
	cudaFree(particleFitness);
	cudaFree(initialTemperature);
	cudaFree(horizontalSplitting);
	cudaFree(referenceCooldownCurveSizeRows);
	cudaFree(referenceCooldownCurve);
	cudaFree(monitoredIndex);
	cudaFree(timeDifference);
	cudaFree(dX);
	cudaFree(tK);
	cudaFree(heatConductivitiesSizeRows);
	cudaFree(heatConductivities);
	cudaFree(specificHeatsSizeRows);
	cudaFree(specificHeats);
	cudaFree(density);
	cudaFree(htcValues);
	cudaFree(currentTemperature);
	cudaFree(previousTemperature);
	cudaFree(g);
	//cudaDeviceReset();

	#pragma endregion
}

void Optimize()
{
	int staticEpochs = 0;

	while (_epoch[0] < _maxEpochs[0] && staticEpochs < _maxStaticEpochs[0])
	{
		double* beta = new double[1];
		beta[0] = (0.9 - 0.55) * (_maxEpochs[0] - _epoch[0]) / _maxEpochs[0] + 0.55;
		bool isErrorImproved = false;
		double* mBest = new double[1];
		mBest[0] = GetSwarmAverageBestPosition();
		OptimizePosition(beta, mBest);

		if (_cudaStatus != cudaSuccess)
		{
			return;
		}

		//WriteIterationResult(_epoch[0]);
		double improovedFitness = 0;
		double* improovedPosition = new double[_dimensionNumber[0]];

		for (int i = 0; i < _particleSwarmSize[0]; i += 1)
		{
			if (_particleFitness[i] < _globalBestFitness[_globalBestSize[0] - 1])
			{
				improovedFitness = _particleBestFitness[i];

				for (int j = 0; j < _dimensionNumber[0]; j += 1)
				{
					improovedPosition[j] = _particleBestPosition[(i * _dimensionNumber[0]) + j];
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
			double* tempFitness = new double[_globalBestSize[0]];
			double* tempPosition = new double[_globalBestSize[0] * _dimensionNumber[0]];

			for (int i = 0; i < _globalBestSize[0]; i += 1)
			{
				tempFitness[i] = _globalBestFitness[i];

				for (int j = 0; j < _dimensionNumber[0]; j += 1)
				{
					tempPosition[(i * _dimensionNumber[0]) + j] = _globalBestPosition[(i * _dimensionNumber[0]) + j];
				}
			}

			delete[] _globalBestFitness;
			delete[] _globalBestPosition;
			_globalBestSize[0] += 1;
			_globalBestFitness = new double[_globalBestSize[0]];
			_globalBestPosition = new double[_globalBestSize[0] * _dimensionNumber[0]];

			for (int i = 0; i < _globalBestSize[0]; i += 1)
			{
				if (i < _globalBestSize[0] - 1)
				{
					_globalBestFitness[i] = tempFitness[i];

					for (int j = 0; j < _dimensionNumber[0]; j += 1)
					{
						_globalBestPosition[(i * _dimensionNumber[0]) + j] = tempPosition[(i * _dimensionNumber[0]) + j];
					}
				}
				else
				{
					_globalBestFitness[i] = improovedFitness;
					std::string positions = "New best found at iteration: ";
					std::stringstream ss1;
					ss1 << _epoch[0];
					positions += ss1.str() + ", with value of: ";
					std::stringstream ss2;
					ss2 << improovedFitness;
					positions += ss2.str() + ", at ";

					for (int j = 0; j < _dimensionNumber[0]; j += 1)
					{
						std::stringstream ss3;
						ss3 << (j + 1);
						std::stringstream ss4;
						ss4 << improovedPosition[j];
						positions += ("position" + ss3.str() + ": " + ss4.str() + " ");
						_globalBestPosition[(i * _dimensionNumber[0]) + j] = improovedPosition[j];
					}

					positions += "\r\n";
					printf(positions.c_str());
				}
			}

			delete[] tempFitness;
			delete[] tempPosition;

			if (fabs(_globalBestFitness[_globalBestSize[0] - 1]) <= _particleEpsilon[0])
			{
				std::stringstream ss;
				ss << _particleEpsilon[0];
				_exitReason = "The particle swarm optimization reached an acceptable fitness value. The value was given at: " + ss.str();

				return;
			}
		}

		delete[] improovedPosition;
		_epoch[0] += 1;

		if (_globalBestSize[0] >= 5)
		{
			if ((_globalBestFitness[_globalBestSize[0] - 1] / _globalBestFitness[_globalBestSize[0] - 5]) < 0.003)
			{
				_exitReason = "The PSO global fitness value changed less then 0.003% over the last 5 global value refresh.";

				return;
			}
		}

		/*if (((double)_globalBestSize[0] / (double)_epoch[0]) < 0.03)
		{
			std::string value = "";

			switch (_particleOptimalisationType[0])
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
		}*/

		if (staticEpochs >= _maxStaticEpochs[0])
		{
			std::stringstream ss;
			ss << _maxStaticEpochs[0];
			_exitReason = "Static iteration limit reached, limit was: " + ss.str();

			return;
		}

		if (_epoch[0] >= _maxEpochs[0])
		{
			std::stringstream ss;
			ss << _maxEpochs[0];
			_exitReason = "Iteration limit reached, limit was: " + ss.str();

			return;
		}
	}
}

void UpdateRing()
{
	int* particleIndex = Shuffle(GetIntegerRange(0, _particleSwarmSize[0]), _particleSwarmSize[0]);

	for (int i = 0; i < _particleSwarmSize[0]; i += 1)
	{
		if (_particleInformers == NULL)
		{
			_particleInformers = new int[_particleSwarmSize[0] * _particleInformerNumber[0]];
		}

		int numberOfinformers = (_particleInformerNumber[0] / 2);
		int currentInformer = 0;

		for (int n = 1; n <= numberOfinformers; n += 1)
		{
			int p = (i - n);

			if (p < 0)
			{
				p = (_particleSwarmSize[0] + p);
			}

			_particleInformers[(i * _particleInformerNumber[0]) + currentInformer] = particleIndex[p];
			currentInformer += 1;
		}

		numberOfinformers += (_particleInformerNumber[0] % 2);

		for (int n = 1; n <= numberOfinformers; n += 1)
		{
			int p = (i + n);

			if (p >= _particleSwarmSize[0])
			{
				p = (p - _particleSwarmSize[0]);
			}

			_particleInformers[(i * _particleInformerNumber[0]) + currentInformer] = particleIndex[p];
			currentInformer += 1;
		}
	}

	delete[] particleIndex;
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
		int n = ((indexSize - i) * _random->rnd2()/*((double)rand() / RAND_MAX)*/ + i);

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

	for (int i = 0; i < _particleSwarmSize[0]; i += 1)
	{
		for (int j = 0; j < _dimensionNumber[0]; j += 1)
		{
			sum += _particleBestPosition[(i * _dimensionNumber[0]) + j];
		}
	}

	return ((sum * 1.0) / _particleSwarmSize[0]);
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

#pragma endregion