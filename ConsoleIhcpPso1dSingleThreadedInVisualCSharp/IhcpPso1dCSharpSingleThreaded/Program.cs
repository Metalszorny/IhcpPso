using System;
using System.Diagnostics;
using System.IO;

namespace IhcpPso1dCSharpSingleThreaded
{
    /// <summary>
    /// Calculates the reference cool-down curve from the hypothetical input data and searches for a similar cool-down curve within the given error value using a one dimensional cool-down model.
    /// </summary>
    class Program
    {
        #region Fields

        #region HTC

        // The initial temperature of the object that's cooling down.
        private static double _initialTemperature;

        // The number of how many parts the body is split horizontally.
        private static int _horizontalSplitting;

        // The number of how many parts the body is split vertically.
        private static int _verticalSplitting;

        // The tK value of the cool-down model.
        private static double _tK;

        // The index of the horizontal part that is monitored.
        private static int _monitoredIndex;

        // The differencial value of time that is incemented in the cooldown.
        private static double _timeDifference;

        // The number of dimensions of the heat transfer coefficients.
        private static int _dimensionNumber;

        // The heat transfer coefficient and temperature values for the reference cooldown curve.
        private static double[] _htcValues;

        // The specific heat and temperature values for the cooldown curves.
        private static double[] _specificHeats;

        // The heat conductivity and temperature values for the cooldown curves.
        private static double[] _heatConductivities;

        // The minimum value of the range for the heat conductivity values.
        private static int _rangeMin;

        // The maximum value of the range for the heat conductivity values.
        private static int _rangeMax;

        // The limit value of time needed for the simulation of the object cooldown.
        private static double _simulationTime;

        // The height of the solid metal cylindrical body.
        private static double _height;

        // The radius of the solid metal cylindrical body.
        private static double _radius;

        // The density of the solid metal cylindrical body.
        private static double _density;

        // The identical distance the radius is split by the horizontal splitting number.
        private static double _dX;

        // The identical distance the height is split by the vertical splitting number.
        private static double _dY;

        // The row count of the heat conductivity and temperature values for the cooldown curves.
        private static int _heatConductivitiesSizeRows;

        // The row count of the specific heat and temperature values for the cooldown curves.
        private static int _specificHeatsSizeRows;

        // The reference cool-down curve.
        private static double[] _referenceCooldownCurve;

        // The row count of the reference cool-down curve.
        private static int _referenceCooldownCurveSizeRows;

        // The calculated cool-down curve.
        private static double[] _calculatedCooldownCurve;

        // The row count of the calculated cool-down curve.
        private static int _calculatedCooldownCurveSizeRows;

        #endregion HTC

        #region PSO

        // The initial velocity values for the particals, these were read for the file.
        private static double[] _particleInitialVelocity;

        // The initial position values for the particals, these were read for the file.
        private static double[] _particleInitialPosition;

        // The number of informers (or neighbors) a particle have.
        private static int _particleInformerNumber;

        // The constant1 number of the Clerk optimization.
        private static double _particleConstant1;

        // The constant2 number of the Clerk optimization.
        private static double _particleConstant2;

        // The fitness limit that the optimization aims to reach.
        private static double _particleEpsilon;

        // The positions of the particles.
        private static double[] _particlePosition;

        // The velocity of the particles.
        private static double[] _particleVelocity;

        // The best positions of the particles.
        private static double[] _particleBestPosition;

        // The fitness values of the particles.
        private static double[] _particleFitness;

        // The best fitness values of the particles.
        private static double[] _particleBestFitness;

        // The size of the global best values.
        private static int _globalBestSize;

        // The global best fitness value of the optimazition.
        private static double[] _globalBestFitness;

        // The global best positions of the optimazition.
        private static double[] _globalBestPosition;

        // The type of the optimization.
        private static int _particleOptimalisationType;

        // The size of the particle swarm.
        private static int _particleSwarmSize;

        // The maximum number of the optimization iterations.
        private static int _maxEpochs;

        // The maximum number of the optimization static iterations.
        private static int _maxStaticEpochs;

        // The weight of the particle.
        private static double _weight;

        // The informers (or neighbors) of the particles.
        private static int[] _particleInformers;

        // The size of the particle's position.
        private static int _particleInitialPositionSize;

        // The size of the particle's velocity.
        private static int _particleInitialVelocitySize;

        // The interation counter of the optimazion.
        private static int _epoch;

        // The stopwatch to mesure the elapsed time during the optimization.
        private static Stopwatch _stopwatch;

        // The termination reason.
        private static string _exitReason;

        // The power exponent for the difference equation.
        private static int _powerExponent;

        #endregion PSO

        // Pseudo random number generator.
        private static Random _random;

        #endregion Fields

        #region Methods

        /// <summary>
        /// The main entry point of the application. This method controls the process of the application flow.
        /// </summary>
        /// <param name="args">The console line input arguments.</param>
        static void Main(string[] args)
        {
            ReadDataFromFile("Resources\\ConfigurationIn.txt");
            CalculateReferenceCooldownCurve();
            WriteReferenceCooldownLogToFile();

            // Set particle swarm initial values.
            for (int i = 0; i < _particleSwarmSize; i++)
            {
                InitializeParticle(i);
            }

            // Set informers.
            UpdateRing();

            // Set initial best values.
            for (int i = 0; i < _particleSwarmSize; i++)
            {
                if (_particleBestFitness[i] < _globalBestFitness[_globalBestSize - 1])
                {
                    _globalBestFitness[_globalBestSize - 1] = _particleBestFitness[i];

                    for (int j = 0; j < _dimensionNumber; j++)
                    {
                        _globalBestPosition[(_globalBestSize - 1) + j] = _particlePosition[(i * _dimensionNumber) + j];
                    }
                }
            }

            _stopwatch.Start();
            Optimize();
            _stopwatch.Stop();
            WritePsoGlobalBestLogToFile();
            WriteExitResultAndTimeLogToFile(_stopwatch.ElapsedMilliseconds);
            CalculateCalculatedCooldownCurve();
            WriteCalculatedCooldownLogToFile();
        }

        /// <summary>
        /// Reads the input values from the given file and initializes the fields used for the calculations.
        /// </summary>
        /// <param name="dataFilePath">The path and name to the input data file.</param>
        private static void ReadDataFromFile(string dataFilePath)
        {
            using (Stream stream = File.Open(dataFilePath, FileMode.Open, FileAccess.Read))
            {
                using (StreamReader streamReader = new StreamReader(stream))
                {
                    string currentLine = streamReader.ReadLine().Trim();

                    while (!string.IsNullOrEmpty(currentLine) &&
                        !string.IsNullOrWhiteSpace(currentLine))
                    {
                        if (currentLine.Contains("#"))
                        {
                            if (currentLine.ToLower().Contains("initial temperature"))
                            {
                                _initialTemperature = Convert.ToDouble(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("horizontal splitting"))
                            {
                                _horizontalSplitting = Convert.ToInt32(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("vertical splitting"))
                            {
                                _verticalSplitting = Convert.ToInt32(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("monitored index"))
                            {
                                _monitoredIndex = Convert.ToInt32(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("power exponent"))
                            {
                                _powerExponent = Convert.ToInt32(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("range min"))
                            {
                                _rangeMin = Convert.ToInt32(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("range max"))
                            {
                                _rangeMax = Convert.ToInt32(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("simulation time"))
                            {
                                _simulationTime = Convert.ToDouble(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("height"))
                            {
                                _height = Convert.ToDouble(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("htcin"))
                            {
                                dataFilePath = streamReader.ReadLine().Trim();
                            }
                            else if (currentLine.ToLower().Contains("time difference"))
                            {
                                _timeDifference = Convert.ToDouble(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("pso epsilon"))
                            {
                                _particleEpsilon = Convert.ToDouble(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("radius"))
                            {
                                _radius = Convert.ToDouble(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("density"))
                            {
                                _density = Convert.ToDouble(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("optimization type"))
                            {
                                string stringDataLine = streamReader.ReadLine().Trim();

                                if (stringDataLine == "clerc")
                                {
                                    _particleOptimalisationType = 1;
                                }
                                else if (stringDataLine == "quantum")
                                {
                                    _particleOptimalisationType = 2;
                                }
                                else
                                {
                                    _particleOptimalisationType = 1;
                                }
                            }
                            else if (currentLine.ToLower().Contains("swarm size"))
                            {
                                _particleSwarmSize = Convert.ToInt32(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("max epochs"))
                            {
                                _maxEpochs = Convert.ToInt32(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("max static epochs"))
                            {
                                _maxStaticEpochs = Convert.ToInt32(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("weight"))
                            {
                                _weight = Convert.ToDouble(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("constant1"))
                            {
                                _particleConstant1 = Convert.ToDouble(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("constant2"))
                            {
                                _particleConstant2 = Convert.ToDouble(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("informer number"))
                            {
                                _particleInformerNumber = Convert.ToInt32(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("initial position"))
                            {
                                _particleInitialPositionSize = Convert.ToInt32(streamReader.ReadLine().Trim().Replace(".", ","));

                                if (_particleInitialPositionSize > 0)
                                {
                                    _particleInitialPosition = new double[_particleInitialPositionSize];

                                    for (int i = 0; i < _particleInitialPositionSize; i++)
                                    {
                                        _particleInitialPosition[i] = Convert.ToDouble(streamReader.ReadLine().Trim().Replace(".", ","));
                                    }
                                }
                            }
                            else if (currentLine.ToLower().Contains("initial velocity"))
                            {
                                _particleInitialVelocitySize = Convert.ToInt32(streamReader.ReadLine().Trim().Replace(".", ","));

                                if (_particleInitialVelocitySize > 0)
                                {
                                    _particleInitialVelocity = new double[_particleInitialVelocitySize];

                                    for (int i = 0; i < _particleInitialVelocitySize; i++)
                                    {
                                        _particleInitialVelocity[i] = Convert.ToDouble(streamReader.ReadLine().Trim().Replace(".", ","));
                                    }
                                }
                            }
                            else if (currentLine.ToLower().Contains("heat conductivity"))
                            {
                                _heatConductivitiesSizeRows = Convert.ToInt32(streamReader.ReadLine().Trim().Replace(".", ","));

                                if (_heatConductivitiesSizeRows > 0)
                                {
                                    _heatConductivities = new double[(_heatConductivitiesSizeRows * 2)];

                                    for (int i = 0; i < (_heatConductivitiesSizeRows * 2); i += 2)
                                    {
                                        string value = streamReader.ReadLine().Trim().Replace(".", ",").Replace("\t", " ");
                                        string[] result = value.Split(' ');
                                        _heatConductivities[i + 0] = Convert.ToDouble(result[0]);
                                        _heatConductivities[i + 1] = Convert.ToDouble(result[1]);
                                    }
                                }
                            }
                            else if (currentLine.ToLower().Contains("specific heat"))
                            {
                                _specificHeatsSizeRows = Convert.ToInt32(streamReader.ReadLine().Trim().Replace(".", ","));

                                if (_specificHeatsSizeRows > 0)
                                {
                                    _specificHeats = new double[(_specificHeatsSizeRows * 2)];

                                    for (int i = 0; i < (_specificHeatsSizeRows * 2); i += 2)
                                    {
                                        string value = streamReader.ReadLine().Trim().Replace(".", ",").Replace("\t", " ");
                                        string[] result = value.Split(' ');
                                        _specificHeats[i + 0] = Convert.ToDouble(result[0]);
                                        _specificHeats[i + 1] = Convert.ToDouble(result[1]);
                                    }
                                }
                            }
                        }

                        currentLine = streamReader.ReadLine();
                    }
                }
            }

            using (Stream stream = File.Open(dataFilePath, FileMode.Open, FileAccess.Read))
            {
                using (StreamReader streamReader = new StreamReader(stream))
                {
                    string currentLine = streamReader.ReadLine().Trim();

                    while (!string.IsNullOrEmpty(currentLine) &&
                        !string.IsNullOrWhiteSpace(currentLine))
                    {
                        if (currentLine.Contains("#"))
                        {
                            if (currentLine.ToLower().Contains("htc values"))
                            {
                                _dimensionNumber = Convert.ToInt32(streamReader.ReadLine());

                                if (_dimensionNumber > 0)
                                {
                                    _htcValues = new double[(_dimensionNumber * 2)];

                                    for (int i = 0; i < (_dimensionNumber * 2); i += 2)
                                    {
                                        string value = streamReader.ReadLine().Trim().Replace(".", ",");
                                        string[] result = value.Split(' ');
                                        _htcValues[i + 0] = Convert.ToDouble(result[0]);
                                        _htcValues[i + 1] = Convert.ToDouble(result[1]);
                                    }
                                }
                            }
                        }

                        currentLine = streamReader.ReadLine();
                    }
                }
            }

            _dX = ((_radius / 1000) / _horizontalSplitting);
            _dY = ((_height / 1000) / _verticalSplitting);
            _tK = 0;
            _stopwatch = new Stopwatch();
            _random = new Random();
            _exitReason = "";
            // Get particle swarm initial values.
            _particlePosition = new double[(_particleSwarmSize * _dimensionNumber)];
            _particleVelocity = new double[(_particleSwarmSize * _dimensionNumber)];
            _particleBestPosition = new double[(_particleSwarmSize * _dimensionNumber)];
            _particleFitness = new double[_particleSwarmSize];
            _particleBestFitness = new double[_particleSwarmSize];
            _globalBestSize = 1;
            _globalBestFitness = new double[_globalBestSize];
            _globalBestFitness[(_globalBestSize - 1)] = double.MaxValue;
            _globalBestPosition = new double[(_globalBestSize * _dimensionNumber)];
        }

        /// <summary>
        /// Sets the default values for the current temperature field.
        /// </summary>
        /// <returns>The current temperature field filled with the default values.</returns>
        private static double[] SetInitialCurrentTemperature()
        {
            double[] returnValue = new double[_horizontalSplitting];

            for (int i = 0; i < _horizontalSplitting; i++)
            {
                returnValue[i] = _initialTemperature;
            }

            return returnValue;
        }

        /// <summary>
        /// Sets the default values for the previous temperature field.
        /// </summary>
        /// <returns>The previous temperature field filled with the default values.</returns>
        private static double[] SetInitialPreviousTemperature()
        {
            double[] returnValue = new double[_horizontalSplitting];

            for (int i = 0; i < _horizontalSplitting; i++)
            {
                returnValue[i] = _initialTemperature;
            }

            return returnValue;
        }

        /// <summary>
        /// Sets the default values for the g field.
        /// </summary>
        /// <returns>The g field filled with the default values.</returns>
        private static double[] SetInitialG()
        {
            double[] returnValue = new double[_horizontalSplitting];

            for (int i = 0; i < _horizontalSplitting; i++)
            {
                returnValue[i] = 0;
            }

            return returnValue;
        }

        /// <summary>
        /// Calculates the reference (or master) cool-down curve.
        /// </summary>
        private static void CalculateReferenceCooldownCurve()
        {
            _referenceCooldownCurveSizeRows = ((int)(_simulationTime / _timeDifference) + 1);
            _referenceCooldownCurve = new double[(_referenceCooldownCurveSizeRows * 2)];
            double[] currentTemperature = SetInitialCurrentTemperature();
            double[] previousTemperature = SetInitialPreviousTemperature();
            double[] g = SetInitialG();

            for (int i = 0; i < (_referenceCooldownCurveSizeRows * 2); i += 2)
            {
                currentTemperature = CalculateCooldownCurve1D(false, -1, false, currentTemperature, previousTemperature, g);
                _referenceCooldownCurve[i + 0] = ((i / 2) * _timeDifference);
                _referenceCooldownCurve[i + 1] = currentTemperature[_monitoredIndex];

                for (int j = 0; j < _horizontalSplitting; j++)
                {
                    previousTemperature[j] = currentTemperature[j];
                }
            }

            currentTemperature = null;
            previousTemperature = null;
            g = null;
        }

        /// <summary>
        /// Calculates the new values of the cool-down's next iteration using the one dimension cool-down model.
        /// </summary>
        /// <param name="isParticle">Indicates if the calculation is on a particle.</param>
        /// <param name="particleIndex">The particle's index in the swarm.</param>
        /// <param name="currentTemperature">The current temperature values of the calculation.</param>
        /// <param name="previousTemperature">The previous temperature values of the calculation.</param>
        /// <param name="g">The g values of the calculation.</param>
        /// <returns>The next iteration's cool-down values.</returns>
        private static double[] CalculateCooldownCurve1D(bool isParticle, int particleIndex, bool useGlobal, double[] currentTemperature, double[] previousTemperature, double[] g)
        {
            double heatConductivity = GetHeatConductivity(currentTemperature[_horizontalSplitting - 1]);
            double alpha = GetAlpha(heatConductivity, GetSpecificHeat(currentTemperature[_horizontalSplitting - 1]));
            // tempTi[0] = tempTi-1[0] + dt * aplha * (1 / (dX * dX) * 2 * (tempTi-1[1] - tempTi-1[0]) + g[0] / HC)
            currentTemperature[0] = previousTemperature[0] + _timeDifference * alpha * (1 / (_dX * _dX) * 2 *
                (previousTemperature[1] - previousTemperature[0]) + g[0] / heatConductivity);
            double heatTransferCoefficient = GetHeatTransferCoefficient(currentTemperature[_horizontalSplitting - 1], isParticle, particleIndex, useGlobal);
            // tempTi[n] = tempTi-1[n] + dt * alpha * (1 / (dX * dX) * 2 * (tempTi-1[n - 1] - tempTi-1[n] - dX / HC * (HTC * (tempTi-1[n] - tK))) + 1 / (horizontalSplitting * dX) * (-1 / HC) * (HTC * (tempTi-1[n] - tK)) + g[n] / HC)
            currentTemperature[_horizontalSplitting - 1] = previousTemperature[_horizontalSplitting - 1] +
                _timeDifference * alpha * (1 / (_dX * _dX) * 2 * (previousTemperature[_horizontalSplitting - 2] -
                previousTemperature[_horizontalSplitting - 1] - _dX / heatConductivity * (heatTransferCoefficient *
                (previousTemperature[_horizontalSplitting - 1] - _tK))) + 1 / (_horizontalSplitting * _dX) *
                (-1 / heatConductivity) * (heatTransferCoefficient * (previousTemperature[_horizontalSplitting - 1] - 
				_tK)) + g[_horizontalSplitting - 1] / heatConductivity);

            for (int i = 1; i < _horizontalSplitting - 1; i++)
            {
                // tempTi[i] = tempTi-1[i] + dt * alpha * (1 / (dX * dX) * (tempTi-1[i - 1] + tempTi-1[i + 1] - 2 * tempTi-1[i]) + 1 / (i * dX) * 1 / (2 * dX) * (tempTi-1[i + 1] - tempTi-1[i - 1]) + g[i] / HC)
                currentTemperature[i] = previousTemperature[i] + _timeDifference * alpha * (1 / (_dX * _dX) *
                    (previousTemperature[i - 1] + previousTemperature[i + 1] - 2 * previousTemperature[i]) + 
                    1 / (i * _dX) * 1 / (2 * _dX) * (previousTemperature[i + 1] - previousTemperature[i - 1]) + 
                    g[i] / heatConductivity);
            }

            return currentTemperature;
        }

        /// <summary>
        /// Gets the heat conductivity value for the given temperature.
        /// </summary>
        /// <param name="temperature">The temperature to match the heat conductivity to.</param>
        /// <returns>The heat conductivity value for the given temperature.</returns>
        private static double GetHeatConductivity(double temperature)
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

                    // Linear-interpolate the return value.
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

        /// <summary>
        /// Gets the specific heat value for the given temperature.
        /// </summary>
        /// <param name="temperature">The temperature to match the specific heat to.</param>
        /// <returns>The specific heat value for the given temperature.</returns>
        private static double GetSpecificHeat(double temperature)
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

                    // Linear-interpolate the return value.
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

        /// <summary>
        /// Gets the alpha value for the given heat conductivity and specific heat.
        /// </summary>
        /// <param name="heatConductivity">The heat conductivity to match the alpha to.</param>
        /// <param name="specificHeat">The specific heat to match the alpha to.</param>
        /// <returns>The alpha value for the given heat conductivity and specific heat.</returns>
        private static double GetAlpha(double heatConductivity, double specificHeat)
        {
            return (heatConductivity / (specificHeat * _density));
        }

        /// <summary>
        /// Gets the heat transfer coefficient value for the given temperature.
        /// </summary>
        /// <param name="temperature">The temperature to match the heat transfer coefficient to.</param>
        /// <param name="isParticle">Indicates if the calculation is on a particle.</param>
        /// <param name="particleIndex">The particle's index in the swarm.</param>
        /// <returns></returns>
        private static double GetHeatTransferCoefficient(double temperature, bool isParticle, int particleIndex, bool useGlobal)
        {
            if (!useGlobal)
            {
                // The calculation is based on the temperature parameter.
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

                            // Linear-interpolate the return value.
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
                // The calculation is based on the particle's position parameter.
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

                            // Linear-interpolate the return value.
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
            // Use the global best position.
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
                        heatTransferCoefficient0 = _globalBestPosition[(_globalBestSize * _dimensionNumber) + (i / 2)];
                        temperature0 = _htcValues[i + 0];
                        i += 2;
                    }

                    if (i < (_dimensionNumber * 2))
                    {
                        heatTransferCoefficient1 = _globalBestPosition[(_globalBestSize * _dimensionNumber) + (i / 2)];
                        temperature1 = _htcValues[i + 0];

                        if (i == 0)
                        {
                            temperature0 = 0;
                            heatTransferCoefficient0 = 0;
                        }

                        // Linear-interpolate the return value.
                        return ((heatTransferCoefficient1 - heatTransferCoefficient0) / (temperature1 - temperature0) * temperature -
                            ((heatTransferCoefficient1 - heatTransferCoefficient0) / (temperature1 - temperature0) * temperature0 - heatTransferCoefficient0));
                    }
                    else
                    {
                        return _globalBestPosition[(_globalBestSize * _dimensionNumber) + (_dimensionNumber - 1)];
                    }
                }

                return 0;
            }
        }

        /// <summary>
        /// Writes the reference cool-down curve values to a file.
        /// </summary>
        private static void WriteReferenceCooldownLogToFile()
        {
            if (_referenceCooldownCurve != null && _referenceCooldownCurveSizeRows > 0)
            {
                using (Stream stream = File.Open("ReferenceCooldownLog.txt", FileMode.Create, FileAccess.Write))
                {
                    using (StreamWriter streamWriter = new StreamWriter(stream))
                    {
                        streamWriter.Write("time ");
                        streamWriter.Write("temperature ");
                        streamWriter.Write("\r\n");

                        for (int i = 0; i < (_referenceCooldownCurveSizeRows * 2); i += 2)
                        {
                            streamWriter.Write(_referenceCooldownCurve[i + 0].ToString("000.000").Replace(",", ".") + " ");
                            streamWriter.Write(_referenceCooldownCurve[i + 1].ToString("000.000").Replace(",", ".") + " ");
                            streamWriter.Write("\r\n");
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Writes the particle swarm optimization's global best fitness and position values to a file.
        /// </summary>
        private static void WritePsoGlobalBestLogToFile()
        {
            if (_globalBestFitness != null && _globalBestPosition != null &&
                _globalBestSize > 0)
            {
                using (Stream stream = File.Open("ParticleSwarmOptimizationGlobalBestLog.txt", FileMode.Create, FileAccess.Write))
                {
                    using (StreamWriter streamWriter = new StreamWriter(stream))
                    {
                        streamWriter.Write("fitness ");

                        for (int i = 0; i < _dimensionNumber; i++)
                        {
                            streamWriter.Write("htc" + (i + 1) + " ");
                        }

                        streamWriter.Write("\r\n");

                        for (int i = 0; i < _globalBestSize; i++)
                        {
                            streamWriter.Write(_globalBestFitness[i].ToString("#.000").Replace(",", ".") + " ");

                            for (int j = 0; j < _dimensionNumber; j++)
                            {
                                streamWriter.Write(_globalBestPosition[(i * _dimensionNumber) + j].ToString("00000.000").Replace(",", ".") + " ");
                            }

                            streamWriter.Write("\r\n");
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Writes the termination reason, the time elapsed of the optimization and the interaation count to a file.
        /// </summary>
        /// <param name="elaspedMilliseconds"></param>
        private static void WriteExitResultAndTimeLogToFile(long elaspedMilliseconds)
        {
            using (Stream stream = File.Open("ExitResultAndTimeLog.txt", FileMode.Create, FileAccess.Write))
            {
                using (StreamWriter streamWriter = new StreamWriter(stream))
                {
                    streamWriter.WriteLine("Exit reason: " + _exitReason + ", elapsed time: " + elaspedMilliseconds + "ms, iterations: " + _epoch);
                }
            }
        }

        /// <summary>
        /// Initializes the particle given by it's index in the swarm.
        /// </summary>
        /// <param name="particleIndex">The particle's index in the swarm.</param>
        private static void InitializeParticle(int particleIndex)
        {
            GetInitialPosition(particleIndex);
            GetInitialVelocity(particleIndex);
            _particleFitness[particleIndex] = ObjectiveFunction(particleIndex);
            _particleBestFitness[particleIndex] = _particleFitness[particleIndex];

            for (int j = 0; j < _dimensionNumber; j++)
            {
                _particleBestPosition[(particleIndex * _dimensionNumber) + j] = _particlePosition[(particleIndex * _dimensionNumber) + j];
            }
        }

        /// <summary>
        /// Calculates the calculated cool-down curve from the particle's position given by the particle's index and summs the difference with the reference cool-down curve.
        /// </summary>
        /// <param name="particleIndex">The particle's index in the swarm.</param>
        /// <returns>The sum of the difference between the reference cool-down curve and the calculated cool-down curve.</returns>
        private static double ObjectiveFunction(int particleIndex)
        {
            int calculatedCooldownCurveSizeRows = _referenceCooldownCurveSizeRows;
            double[] calculatedCooldownCurve = new double[calculatedCooldownCurveSizeRows * 2];
            double[] currentTemperature = SetInitialCurrentTemperature();
            double[] previousTemperature = SetInitialPreviousTemperature();
            double[] g = SetInitialG();

            for (int i = 0; i < (calculatedCooldownCurveSizeRows * 2); i += 2)
            {
                currentTemperature = CalculateCooldownCurve1D(true, particleIndex, false, currentTemperature, previousTemperature, g);
                calculatedCooldownCurve[i + 0] = ((i / 2) * _timeDifference);
                calculatedCooldownCurve[i + 1] = currentTemperature[_monitoredIndex];

                for (int j = 0; j < _horizontalSplitting; j++)
                {
                    previousTemperature[j] = currentTemperature[j];
                }
            }

            currentTemperature = null;
            previousTemperature = null;
            g = null;

            if (calculatedCooldownCurve != null && calculatedCooldownCurveSizeRows > 0)
            {
                if (_referenceCooldownCurveSizeRows == calculatedCooldownCurveSizeRows)
                {
                    double sum = 0;

                    for (int i = 0; i < _referenceCooldownCurveSizeRows * 2; i += 2)
                    {
                        sum += Math.Pow(_referenceCooldownCurve[i + 1] - calculatedCooldownCurve[i + 1], _powerExponent);
                    }

                    calculatedCooldownCurve = null;

                    return sum;
                }
                else
                {
                    calculatedCooldownCurve = null;

                    return double.MaxValue;
                }
            }
            else
            {
                calculatedCooldownCurve = null;

                return double.MaxValue;
            }
        }

        /// <summary>
        /// Optimizes the particle's postions and checks for a termination reason.
        /// </summary>
        private static void Optimize()
        {
            _epoch = 0;
            int staticEpochs = 0;

            while (_epoch < _maxEpochs && staticEpochs < _maxStaticEpochs)
            {
                double beta = (0.9 - 0.55) * (_maxEpochs - _epoch) / _maxEpochs + 0.55;
                bool isErrorImproved = false;
                double mBest = GetSwarmAverageBestPosition();

                for (int i = 0; i < _particleSwarmSize; i++)
                {
                    OptimizePosition(i, mBest, beta);
                }

                double improovedFitness = 0;
                double[] improovedPosition = new double[_dimensionNumber];

                for (int i = 0; i < _particleSwarmSize; i++)
                {
                    if (_particleFitness[i] < _globalBestFitness[_globalBestSize - 1])
                    {
                        improovedFitness = _particleBestFitness[i];

                        for (int j = 0; j < _dimensionNumber; j++)
                        {
                            improovedPosition[j] = _particleBestPosition[(i * _dimensionNumber) + j];
                        }

                        isErrorImproved = true;
                        staticEpochs = 0;
                    }
                }

                if (!isErrorImproved)
                {
                    staticEpochs++;
                }
                else
                {
                    double[] tempFitness = new double[_globalBestSize];
                    double[] tempPosition = new double[_globalBestSize * _dimensionNumber];

                    for (int i = 0; i < _globalBestSize; i++)
                    {
                        tempFitness[i] = _globalBestFitness[i];

                        for (int j = 0; j < _dimensionNumber; j++)
                        {
                            tempPosition[(i * _dimensionNumber) + j] = _globalBestPosition[(i * _dimensionNumber) + j];
                        }
                    }

                    _globalBestFitness = null;
                    _globalBestPosition = null;
                    _globalBestSize++;
                    _globalBestFitness = new double[_globalBestSize];
                    _globalBestPosition = new double[_globalBestSize * _dimensionNumber];

                    for (int i = 0; i < _globalBestSize; i++)
                    {
                        if (i < _globalBestSize - 1)
                        {
                            _globalBestFitness[i] = tempFitness[i];

                            for (int j = 0; j < _dimensionNumber; j++)
                            {
                                _globalBestPosition[(i * _dimensionNumber) + j] = tempPosition[(i * _dimensionNumber) + j];
                            }
                        }
                        else
                        {
                            _globalBestFitness[i] = improovedFitness;

                            for (int j = 0; j < _dimensionNumber; j++)
                            {
                                _globalBestPosition[(i * _dimensionNumber) + j] = improovedPosition[j];
                            }
                        }
                    }

                    string positions = "";

                    for (int i = 0; i < _dimensionNumber; i++)
                    {
                        positions += "position" + (i + 1) +  ": " + improovedPosition[i] + ", ";
                    }

                    Console.WriteLine("New best found: " + improovedFitness + " with " + positions + "at iteration: " + _epoch + ", count: " + _globalBestSize);
                    tempFitness = null;
                    tempPosition = null;
                    
                    if (Math.Abs(_globalBestFitness[_globalBestSize - 1]) <= _particleEpsilon)
                    {
                        _exitReason = "The particle swarm optimization reached an acceptable fitness value. The value was given at: " + _particleEpsilon;

                        return;
                    }
                }

                improovedPosition = null;
                _epoch++;

                if (_globalBestSize >= 5)
                {
                    if ((_globalBestFitness[_globalBestSize - 1] / _globalBestFitness[_globalBestSize - 5]) < 0.003)
                    {
                        _exitReason = "The PSO global fitness value changed less then 0.003% over the last 5 global value refresh.";

                        return;
                    }
                }

                if (((_globalBestSize * 1.0) / (_epoch * 1.0)) < 0.03)
                {
                    _exitReason = "The convergence of the " + (_particleOptimalisationType == 1 ? "Clerc" : "Quantum") + " PSO is too slow, it may not find the optimum";

                    return;
                }

                if (staticEpochs >= _maxStaticEpochs)
                {
                    _exitReason = "Static iteration limit reached, limit was: " + _maxStaticEpochs;

                    return;
                }

                if (_epoch >= _maxEpochs)
                {
                    _exitReason = "Iteration limit reached, limit was: " + _maxEpochs;

                    return;
                }
            }
        }

        /// <summary>
        /// Sets the neighbors for the particles.
        /// </summary>
        private static void UpdateRing()
        {
            int[] particleIndex = Shuffle(GetIntegerRange(0, _particleSwarmSize), _particleSwarmSize);

            for (int i = 0; i < _particleSwarmSize; i++)
            {
                if (_particleInformers == null)
                {
                    _particleInformers = new int[_particleSwarmSize * _particleInformerNumber];
                }

                int numberOfinformers = (_particleInformerNumber / 2);
                int currentInformer = 0;

                for (int n = 1; n <= numberOfinformers; n++)
                {
                    int p = (i - n);

                    if (p < 0)
                    {
                        p = (_particleSwarmSize + p);
                    }
                    
                    _particleInformers[(i * _particleInformerNumber) + currentInformer] = particleIndex[p];
                    currentInformer++;
                }

                numberOfinformers += (_particleInformerNumber % 2);

                for (int n = 1; n <= numberOfinformers; n++)
                {
                    int p = i + n;

                    if (p >= _particleSwarmSize)
                    {
                        p = (p - _particleSwarmSize);
                    }
                    
                    _particleInformers[(i * _particleInformerNumber) + currentInformer] = particleIndex[p];
                    currentInformer++;
                }
            }

            particleIndex = null;
        }

        /// <summary>
        /// Optimizes the position of the particle given by the particle's index.
        /// </summary>
        /// <param name="particleIndex">The particle's index in the swarm.</param>
        /// <param name="mbest">The modus of the global best positions.</param>
        /// <param name="beta">The beta value.</param>
        private static void OptimizePosition(int particleIndex, double mbest, double beta)
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

        /// <summary>
        /// Updates the particle's position and velocity if the optimization type is Clerc.
        /// </summary>
        /// <param name="particleIndex">The particle's index in the swarm.</param>
        private static void UpdateClercPosition(int particleIndex)
        {
            UpdateVelocity(particleIndex, GetBestLocalPositions(particleIndex));
            UpdatePosition(particleIndex);
        }

        /// <summary>
        /// Updates the position's fitness value.
        /// </summary>
        /// <param name="particleIndex">The particle's index in the swarm.</param>
        private static void UpdateErrorValues(int particleIndex)
        {
            double[] checkPosition = new double[_dimensionNumber];

            for (int i = 0; i < _dimensionNumber; i++)
            {
                checkPosition[i] = _particlePosition[(particleIndex * _dimensionNumber) + i];
            }

            if (CheckIsInRange(checkPosition))
            {
                _particleFitness[particleIndex] = ObjectiveFunction(particleIndex);
                
                if (_particleFitness[particleIndex] < _particleBestFitness[particleIndex])
                {
                    _particleBestFitness[particleIndex] = _particleFitness[particleIndex];

                    for (int i = 0; i < _dimensionNumber; i++)
                    {
                        _particleBestPosition[(particleIndex * _dimensionNumber) + i] = _particlePosition[(particleIndex * _dimensionNumber) +  i];
                    }
                }
            }

            checkPosition = null;
        }

        /// <summary>
        /// Updates the particle's position.
        /// </summary>
        /// <param name="particleIndex">The particle's index in the swarm.</param>
        private static void UpdatePosition(int particleIndex)
        {
            for (int i = 0; i < _dimensionNumber; i++)
            {
                _particlePosition[(particleIndex * _dimensionNumber) + i] = 
                    (_particlePosition[(particleIndex * _dimensionNumber ) + i] + 
                    _particleVelocity[(particleIndex * _dimensionNumber) + i]);
            }
        }

        /// <summary>
        /// Updates the particle's position if the optimization type is Quantum.
        /// </summary>
        /// <param name="particleIndex">The particle's index in the swarm.</param>
        /// <param name="mbest">The modus of the global best positions.</param>
        /// <param name="beta">The beta value.</param>
        private static void UpdateQuantumPosition(int particleIndex, double mbest, double beta)
        {
            for (int i = 0; i < _dimensionNumber; i++)
            {
                double fi = _random.NextDouble();
                double p = fi * _particleBestPosition[(particleIndex * _dimensionNumber) + i] + 
                    (1 - fi) * _globalBestPosition[(_globalBestSize - 1 - _dimensionNumber) + i];

                if (fi > 0.5)
                {
                    _particlePosition[(particleIndex * _dimensionNumber) + i] = p - beta * 
                        Math.Abs(mbest - _particlePosition[(particleIndex * _dimensionNumber) + i]) * (-Math.Log10(fi));
                }
                else
                {
                    _particlePosition[(particleIndex * _dimensionNumber) + i] = p + beta * 
                        Math.Abs(mbest - _particlePosition[(particleIndex * _dimensionNumber) + i]) * (-Math.Log10(fi));
                }

                if (_particlePosition[(particleIndex * _dimensionNumber) + i] < _rangeMin)
                {
                    _particlePosition[(particleIndex * _dimensionNumber) + i] = 2 * 
                        _rangeMin - _particlePosition[(particleIndex * _dimensionNumber) + i];
                }

                if (_particlePosition[(particleIndex * _dimensionNumber) + i] > _rangeMax)
                {
                    _particlePosition[(particleIndex * _dimensionNumber) + i] = 2 * 
                        _rangeMax - _particlePosition[(particleIndex * _dimensionNumber) + i];
                }
            }
        }

        /// <summary>
        /// Updates the particle's velocity.
        /// </summary>
        /// <param name="particleIndex">The particle's index in the swarm.</param>
        /// <param name="bestLocalPosition">The positions of the particle with the best fitness value in the neigborhood.</param>
        private static void UpdateVelocity(int particleIndex, double[] bestLocalPosition)
        {
            for (int i = 0; i < _dimensionNumber; i++)
            {
                // particle.velocity[i] = (w * paritcle.velocity[i]) + (constant1 * rand() * (particle.bestPosition[i] - particle.position[i])) + (constant2 * rand() * bestPosition[i] - particle.position[i])
                _particleVelocity[(particleIndex * _dimensionNumber) + i] = 
                    (_weight * _particleVelocity[(particleIndex * _dimensionNumber) + i]) +
                    (_particleConstant1 * _random.NextDouble() * 
                    (_particleBestPosition[(particleIndex * _dimensionNumber) + i] -
                    _particlePosition[(particleIndex * _dimensionNumber) + i])) + 
                    (_particleConstant2 * _random.NextDouble() *
                    (bestLocalPosition[i] - _particlePosition[(particleIndex * _dimensionNumber) + i]));
            }
        }

        /// <summary>
        /// Checks if all of the given position values are in the specified range.
        /// </summary>
        /// <param name="position">The position values that needs to be checked.</param>
        /// <returns>True if the position values are in range, false if at least one of the position value is not in the specified range.</returns>
        private static bool CheckIsInRange(double[] position)
        {
			for (int i = 0; i < _dimensionNumber; i++)
            {
                if (position[i] < _rangeMin || position[i] > _rangeMax)
                {
                    return false;
                }
            }

            return true;
        }

        /// <summary>
        /// Gets the positions of a particle with the best fitness value from the neighborhood.
        /// </summary>
        /// <param name="particleIndex">The particle's index in the swarm.</param>
        /// <returns>The positions of a particle with the best fitness value from the neighborhood.</returns>
        private static double[] GetBestLocalPositions(int particleIndex)
        {
            double[] returnValue = new double[_dimensionNumber];
            int bestIndex = particleIndex;
            
            for (int i = 0; i < _particleInformerNumber; i++)
            {
                if (_particleBestFitness[_particleInformers[(particleIndex * _particleInformerNumber) + i]] < _particleBestFitness[bestIndex])
                {
                    bestIndex = _particleInformers[(particleIndex * _particleInformerNumber) + i];
                }
            }

            for (int i = 0; i < _dimensionNumber; i++)
            {
                returnValue[i] = _particleBestPosition[(bestIndex * _dimensionNumber) + i];
            }

            return returnValue;
        }

        /// <summary>
        /// Sets the default value of the particle's position given by the particle index.
        /// </summary>
        /// <param name="particleIndex">The particle's index in the swarm.</param>
        private static void GetInitialPosition(int particleIndex)
        {
            for (int i = 0; i < _dimensionNumber; i++)
            {
                // The given values are null.
                if (_particleInitialPosition == null)
                {
                    _particlePosition[(particleIndex * _dimensionNumber) + i] = ((_rangeMax - _rangeMin) * _random.NextDouble() + _rangeMin);
                }
                // There is only one given value.
                else if (_particleInitialPositionSize == 1)
                {
                    _particlePosition[(particleIndex * _dimensionNumber) + i] = _particleInitialPosition[0];
                }
                // There are as many given value as the position dimensions.
                else if (_particleInitialPositionSize == _dimensionNumber)
                {
                    _particlePosition[(particleIndex * _dimensionNumber) + i] = _particleInitialPosition[i];
                }
                // The current position can be set from the given values.
                else if (i < _particleInitialPositionSize)
                {
                    _particlePosition[(particleIndex * _dimensionNumber) + i] = _particleInitialPosition[i];
                }
                // The current position can't be set from the given values.
                else
                {
                    _particlePosition[(particleIndex * _dimensionNumber) + i] = ((_rangeMax - _rangeMin) * _random.NextDouble() + _rangeMin);
                }
            }
        }

        /// <summary>
        /// Sets the default value of the particle's velocity given by the particle index.
        /// </summary>
        /// <param name="particleIndex">The particle's index in the swarm.</param>
        private static void GetInitialVelocity(int particleIndex)
        {
            for (int i = 0; i < _dimensionNumber; i++)
            {
                // The given values are null.
                if (_particleInitialVelocity == null)
                {
                    _particleVelocity[(particleIndex * _dimensionNumber) + i] = ((_rangeMax - _rangeMin) * _random.NextDouble() + _rangeMin);
                }
                //  There is only one given value.
                else if (_particleInitialVelocitySize == 1)
                {
                    _particleVelocity[(particleIndex * _dimensionNumber) + i] = _particleInitialVelocity[0];
                }
                // There are as many given value as the velocity dimensions.
                else if (_particleInitialVelocitySize == _dimensionNumber)
                {
                    _particleVelocity[(particleIndex * _dimensionNumber) + i] = _particleInitialVelocity[i];
                }
                // The current velocity can be set from the given values.
                else if (i < _particleInitialVelocitySize)
                {
                    _particleVelocity[(particleIndex * _dimensionNumber) + i] = _particleInitialVelocity[i];
                }
                // The current velocity can't be set from the given values.
                else
                {
                    _particleVelocity[(particleIndex * _dimensionNumber) + i] = ((_rangeMax - _rangeMin) * _random.NextDouble() + _rangeMin);
                }
            }
        }

        /// <summary>
        /// Creates an array of 32 bit integer numbers in the given range.
        /// </summary>
        /// <param name="startIndex">The starting number from where to start the range.</param>
        /// <param name="count">The number of numbers to include in the range.</param>
        /// <returns>An array of 32 bit integer numbers in the given range.</returns>
        private static int[] GetIntegerRange(int startIndex, int count)
        {
            int[] returnValue = new int[count];

            for (int i = 0; i < count; i++)
            {
                returnValue[i] = (startIndex + i);
            }

            return returnValue;
        }

        /// <summary>
        /// Shuffles the given array of 32 bit integers to a random order.
        /// </summary>
        /// <param name="index">The array of 32 bit integer numbers to shuffle.</param>
        /// <param name="indexSize">The length of the array.</param>
        /// <returns>The shuffled array of 32 bit integer numbers.</returns>
        private static int[] Shuffle(int[] index, int indexSize)
        {
            for (int i = 0; i < indexSize; i++)
            {
                int n = _random.Next(i, indexSize);
                int temp = index[n];
                index[n] = index[i];
                index[i] = temp;
            }
            
            return index;
        }

        /// <summary>
        /// Gets the modus of the best positions of the particles.
        /// </summary>
        /// <returns>The modus of the best positions of the particles.</returns>
        private static double GetSwarmAverageBestPosition()
        {
            double sum = 0;

            for (int i = 0; i < _particleSwarmSize; i++)
            {
                for (int j = 0; j < _dimensionNumber; j++)
                {
                    sum += _particleBestPosition[(i * _dimensionNumber) + j];
                }
            }

            return ((sum * 1.0) / _particleSwarmSize);
        }

        /// <summary>
        /// Calculates the calculated cool-down curve.
        /// </summary>
        private static void CalculateCalculatedCooldownCurve()
        {
            _calculatedCooldownCurveSizeRows = _referenceCooldownCurveSizeRows;
            _calculatedCooldownCurve = new double[_calculatedCooldownCurveSizeRows * 2];
            double[] currentTemperature = SetInitialCurrentTemperature();
            double[] previousTemperature = SetInitialPreviousTemperature();
            double[] g = SetInitialG();

            for (int i = 0; i < (_calculatedCooldownCurveSizeRows * 2); i += 2)
            {
                currentTemperature = CalculateCooldownCurve1D(false, -1, true, currentTemperature, previousTemperature, g);
                _calculatedCooldownCurve[i + 0] = ((i / 2) * _timeDifference);
                _calculatedCooldownCurve[i + 1] = currentTemperature[_monitoredIndex];

                for (int j = 0; j < _horizontalSplitting; j++)
                {
                    previousTemperature[j] = currentTemperature[j];
                }
            }

            currentTemperature = null;
            previousTemperature = null;
            g = null;
        }

        /// <summary>
        /// Writes the calculated cool-down curve values to a file.
        /// </summary>
        private static void WriteCalculatedCooldownLogToFile()
        {
            if (_calculatedCooldownCurve != null && _calculatedCooldownCurveSizeRows > 0)
            {
                using (Stream stream = File.Open("CalculatedCooldownLog.txt", FileMode.Create, FileAccess.Write))
                {
                    using (StreamWriter streamWriter = new StreamWriter(stream))
                    {
                        streamWriter.Write("time ");
                        streamWriter.Write("temperature ");
                        streamWriter.Write("\r\n");

                        for (int i = 0; i < (_calculatedCooldownCurveSizeRows * 2); i += 2)
                        {
                            streamWriter.Write(_calculatedCooldownCurve[i + 0].ToString("000.000").Replace(",", ".") + " ");
                            streamWriter.Write(_calculatedCooldownCurve[i + 1].ToString("000.000").Replace(",", ".") + " ");
                            streamWriter.Write("\r\n");
                        }
                    }
                }
            }
        }

        #endregion Methods
    }
}
