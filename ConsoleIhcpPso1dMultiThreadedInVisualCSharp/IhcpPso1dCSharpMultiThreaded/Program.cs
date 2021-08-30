using System;
using System.Diagnostics;
using System.IO;
using System.Threading.Tasks;

namespace IhcpPso1dCSharpMultiThreaded
{
    /// <summary>
    /// 
    /// </summary>
    class Program
    {
        #region Fields

        #region HTC

        // 
        private static double _initialTemperature;

        // 
        private static int _horizontalSplitting;

        // 
        private static int _verticalSplitting;

        // 
        private static int _monitoredIndex;

        // 
        private static double _timeDifference;

        // 
        private static int _dimensionNumber;

        // 
        private static double[] _htcValues;

        // 
        private static double[] _specificHeats;

        // 
        private static double[] _heatConductivities;

        // 
        private static int _rangeMin;

        // 
        private static int _rangeMax;

        // 
        private static double _simulationTime;

        // 
        private static double _height;

        // 
        private static double _radius;

        // 
        private static double _density;

        // 
        private static double _dX;

        // 
        private static double _dY;

        // 
        private static double[] _referenceCooldownCurve;

        // 
        private static double _tK;

        // 
        private static int _heatConductivitiesSizeRows;

        // 
        private static int _specificHeatsSizeRows;

        // 
        private static int _referenceCooldownCurveSizeRows;

        #endregion HTC

        #region PSO

        // The initial velocity values for the particals, these were read for the file.
        private static double[] _particleInitialVelocity;

        // The initial position values for the particals, these were read for the file.
        private static double[] _particleInitialPosition;

        // 
        private static int _particleInformerNumber;

        // 
        private static double _particleConstant1;

        // 
        private static double _particleConstant2;

        // 
        private static double _particleEpsilon;

        // 
        private static double[] _particlePosition;

        // 
        private static double[] _particleVelocity;

        // 
        private static double[] _particleBestPosition;

        // 
        private static double[] _particleFitness;

        // 
        private static double[] _particleBestFitness;

        // 
        private static int _globalBestSize;

        // 
        private static double[] _globalBestFitness;

        // 
        private static double[] _globalBestPosition;

        // 
        private static int _particleOptimalisationType;

        // 
        private static int _particleSwarmSize;

        // 
        private static int _maxEpochs;

        // 
        private static int _maxStaticEpochs;

        // 
        private static double _weight;

        // 
        private static int[] _particleInformers;

        // 
        private static int _particleInitialPositionSize;

        // 
        private static int _particleInitialVelocitySize;

        // 
        private static int _epoch;

        #endregion PSO

        // 
        private static Random _random;

        // 
        private static Stopwatch _stopwatch;

        // 
        private static string _exitReason;

        #endregion Fields

        #region Methods

        /// <summary>
        /// 
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args)
        {
            ReadDataFromFile("Resources\\ConfigurationIn.txt");
            CalculateReferenceCooldownCurve();
            WriteReferenceCooldownLogToFile();
            Task<Particle>[] taskArray = new Task<Particle>[_particleSwarmSize];

            // Set particle swarm initial values.
            for (int i = 0; i < _particleSwarmSize; i++)
            {
                Particle currentParticle = new Particle();
                currentParticle.bestFitness = _particleBestFitness[i];
                currentParticle.index = i;
                currentParticle.fitness = _particleFitness[i];
                currentParticle.bestPosition = new double[_dimensionNumber];
                currentParticle.position = new double[_dimensionNumber];
                currentParticle.velocity = new double[_dimensionNumber];

                for (int j = 0; j < _dimensionNumber; j++)
                {
                    currentParticle.bestPosition[j] = _particleBestPosition[(i * _dimensionNumber) + j];
                    currentParticle.position[j] = _particlePosition[(i * _dimensionNumber) + j];
                    currentParticle.velocity[j] = _particleVelocity[(i * _dimensionNumber) + j];
                }

                Task<Particle> newTask = new Task<Particle>(() =>
                {
                    return InitializeParticle(currentParticle);
                });
                newTask.Start();
                taskArray[i] = newTask;
            }

            Task.WaitAll(taskArray);

            for (int i = 0; i < _particleSwarmSize; i++)
            {
                Particle currentParticle = taskArray[i].Result;
                _particleFitness[currentParticle.index] = currentParticle.fitness;
                _particleBestFitness[currentParticle.index] = currentParticle.bestFitness;

                for (int j = 0; j < _dimensionNumber; j++)
                {
                    _particlePosition[(currentParticle.index * _dimensionNumber) + j] = currentParticle.position[j];
                    _particleVelocity[(currentParticle.index * _dimensionNumber) + j] = currentParticle.velocity[j];
                    _particleBestPosition[(currentParticle.index * _dimensionNumber) + j] = currentParticle.bestPosition[j];
                }
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
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="dataFilePath"></param>
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
                                    _heatConductivities = new double[_heatConductivitiesSizeRows * 2];

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
                                    _specificHeats = new double[_specificHeatsSizeRows * 2];

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
                                    _htcValues = new double[_dimensionNumber * 2];

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

            _dX = (_radius / 1000) / _horizontalSplitting;
            _dY = (_height / 1000) / _verticalSplitting;
            _tK = 0;
			_stopwatch = new Stopwatch();
            _random = new Random();
            _exitReason = "";
			// Get particle swarm initial values.
            _particlePosition = new double[_particleSwarmSize * _dimensionNumber];
            _particleVelocity = new double[_particleSwarmSize * _dimensionNumber];
            _particleBestPosition = new double[_particleSwarmSize * _dimensionNumber];
            _particleFitness = new double[_particleSwarmSize];
            _particleBestFitness = new double[_particleSwarmSize];
            _globalBestSize = 1;
            _globalBestFitness = new double[_globalBestSize];
            _globalBestFitness[_globalBestSize - 1] = double.MaxValue;
            _globalBestPosition = new double[_globalBestSize * _dimensionNumber];
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
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
        /// 
        /// </summary>
        /// <returns></returns>
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
        /// 
        /// </summary>
        /// <returns></returns>
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
        /// 
        /// </summary>
        private static void CalculateReferenceCooldownCurve()
        {
            _referenceCooldownCurveSizeRows = (int)(_simulationTime / _timeDifference) + 1;
            _referenceCooldownCurve = new double[_referenceCooldownCurveSizeRows * 2];
            double[] currentTemperature = SetInitialCurrentTemperature();
            double[] previousTemperature = SetInitialPreviousTemperature();
            double[] g = SetInitialG();

            for (int i = 0; i < (_referenceCooldownCurveSizeRows * 2); i += 2)
            {
                currentTemperature = CalculateCooldownCurve1D(false, -1, currentTemperature, previousTemperature, g);
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
        /// 
        /// </summary>
        /// <param name="isParticle"></param>
        /// <param name="particleIndex"></param>
        /// <param name="currentTemperature"></param>
        /// <param name="previousTemperature"></param>
        /// <param name="g"></param>
        /// <returns></returns>
        private static double[] CalculateCooldownCurve1D(bool isParticle, int particleIndex, double[] currentTemperature, double[] previousTemperature, double[] g)
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
        /// 
        /// </summary>
        /// <param name="currentParticle"></param>
        /// <param name="currentTemperature"></param>
        /// <param name="previousTemperature"></param>
        /// <param name="g"></param>
        /// <returns></returns>
        private static double[] CalculateCooldownCurve1D(Particle currentParticle, double[] currentTemperature, double[] previousTemperature, double[] g)
        {
            double heatConductivity = GetHeatConductivity(currentTemperature[_horizontalSplitting - 1]);
            double alpha = GetAlpha(heatConductivity, GetSpecificHeat(currentTemperature[_horizontalSplitting - 1]));
            // tempTi[0] = tempTi-1[0] + dt * aplha * (1 / (dX * dX) * 2 * (tempTi-1[1] - tempTi-1[0]) + g[0] / HC)
            currentTemperature[0] = previousTemperature[0] + _timeDifference * alpha * (1 / (_dX * _dX) * 2 *
                (previousTemperature[1] - previousTemperature[0]) + g[0] / heatConductivity);
            double heatTransferCoefficient = GetHeatTransferCoefficient(currentTemperature[_horizontalSplitting - 1], currentParticle);
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
        /// 
        /// </summary>
        /// <param name="temperature"></param>
        /// <returns></returns>
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
        /// 
        /// </summary>
        /// <param name="temperature"></param>
        /// <returns></returns>
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
        /// 
        /// </summary>
        /// <param name="heatConductivity"></param>
        /// <param name="specificHeat"></param>
        /// <returns></returns>
        private static double GetAlpha(double heatConductivity, double specificHeat)
        {
            return (heatConductivity / (specificHeat * _density));
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="temperature"></param>
        /// <param name="isParticle"></param>
        /// <param name="particleIndex"></param>
        /// <returns></returns>
        private static double GetHeatTransferCoefficient(double temperature, bool isParticle, int particleIndex)
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

        /// <summary>
        /// 
        /// </summary>
        /// <param name="temperature"></param>
        /// <param name="currentParticle"></param>
        /// <returns></returns>
        private static double GetHeatTransferCoefficient(double temperature, Particle currentParticle)
        {
            if (currentParticle == null)
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
                        heatTransferCoefficient0 = currentParticle.position[(i / 2)];
                        temperature0 = _htcValues[i + 0];
                        i += 2;
                    }

                    if (i < (_dimensionNumber * 2))
                    {
                        heatTransferCoefficient1 = currentParticle.position[(i / 2)];
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
                        return currentParticle.position[(_dimensionNumber - 1)];
                    }
                }

                return 0;
            }
        }

        /// <summary>
        /// 
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
        /// 
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
        /// 
        /// </summary>
        /// <param name="exitResult"></param>
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
        /// 
        /// </summary>
        /// <param name="currentParticle"></param>
        /// <returns></returns>
        private static Particle InitializeParticle(Particle currentParticle)
        {
            currentParticle = GetInitialVelocity(GetInitialPosition(currentParticle));
            currentParticle.fitness = ObjectiveFunction(currentParticle);
            currentParticle.bestFitness = currentParticle.fitness;

            for (int j = 0; j < _dimensionNumber; j++)
            {
                currentParticle.bestPosition[j] = currentParticle.position[j];
            }

            return currentParticle;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="currentParticle"></param>
        /// <returns></returns>
        public static double ObjectiveFunction(Particle currentParticle)
        {
            int calculatedCooldownCurveSizeRows = _referenceCooldownCurveSizeRows;
            double[] calculatedCooldownCurve = new double[calculatedCooldownCurveSizeRows * 2];
            double[] currentTemperature = SetInitialCurrentTemperature();
            double[] previousTemperature = SetInitialPreviousTemperature();
            double[] g = SetInitialG();

            for (int i = 0; i < (calculatedCooldownCurveSizeRows * 2); i += 2)
            {
                currentTemperature = CalculateCooldownCurve1D(currentParticle, currentTemperature, previousTemperature, g);
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

                    for (int i = 0; i < (_referenceCooldownCurveSizeRows * 2); i += 2)
                    {
                        sum += Math.Pow(Convert.ToDouble(_referenceCooldownCurve[i + 1].ToString("000.0")) -
                            Convert.ToDouble(calculatedCooldownCurve[i + 1].ToString("000.0")), 2);
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
        /// 
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
                Task<Particle>[] taskArray = new Task<Particle>[_particleSwarmSize];

                for (int i = 0; i < _particleSwarmSize; i++)
                {
                    Particle currentParticle = new Particle();
                    currentParticle.bestFitness = _particleBestFitness[i];
                    currentParticle.index = i;
                    currentParticle.fitness = _particleFitness[i];
                    currentParticle.bestPosition = new double[_dimensionNumber];
                    currentParticle.position = new double[_dimensionNumber];
                    currentParticle.velocity = new double[_dimensionNumber];
                    
                    for (int j = 0; j < _dimensionNumber; j++)
                    {
                        currentParticle.bestPosition[j] = _particleBestPosition[(i * _dimensionNumber) + j];
                        currentParticle.position[j] = _particlePosition[(i * _dimensionNumber) + j];
                        currentParticle.velocity[j] = _particleVelocity[(i * _dimensionNumber) + j];
                    }

                    Task<Particle> newTask = new Task<Particle>(() =>
                    {
                        return OptimizePosition(currentParticle, mBest, beta);
                    });
                    newTask.Start();
                    taskArray[i] = newTask;
                }

                Task.WaitAll(taskArray);

                for (int i = 0; i < _particleSwarmSize; i++)
                {
                    Particle currentParticle = taskArray[i].Result;
                    _particleFitness[currentParticle.index] = currentParticle.fitness;
                    _particleBestFitness[currentParticle.index] = currentParticle.bestFitness;

                    for (int j = 0; j < _dimensionNumber; j++)
                    {
                        _particlePosition[(currentParticle.index * _dimensionNumber) + j] = currentParticle.position[j];
                        _particleVelocity[(currentParticle.index * _dimensionNumber) + j] = currentParticle.velocity[j];
                        _particleBestPosition[(currentParticle.index * _dimensionNumber) + j] = currentParticle.bestPosition[j];
                    }
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
                        positions += "position" + (i + 1) + ": " + improovedPosition[i] + ", ";
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

                if (((double)_globalBestSize / (double)_epoch) < 0.03)
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
        /// 
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
                    int p = (i + n);

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
        /// 
        /// </summary>
        /// <param name="currentParticle"></param>
        /// <param name="mbest"></param>
        /// <param name="beta"></param>
        /// <returns></returns>
        private static Particle OptimizePosition(Particle currentParticle, double mbest, double beta)
        {
            switch (_particleOptimalisationType)
            {
                case 1:
                    currentParticle = UpdateClercPosition(currentParticle);
                    break;
                case 2:
                    currentParticle = UpdateQuantumPosition(currentParticle, mbest, beta);
                    break;
            }

            return UpdateErrorValues(currentParticle);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="currentParticle"></param>
        /// <returns></returns>
        private static Particle UpdateClercPosition(Particle currentParticle)
        {
            return UpdatePosition(UpdateVelocity(currentParticle, GetBestLocalPositions(currentParticle.index)));
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="currentParticle"></param>
        /// <returns></returns>
        private static Particle UpdateErrorValues(Particle currentParticle)
        {
            if (CheckIsInRange(currentParticle.position))
            {
                currentParticle.fitness = ObjectiveFunction(currentParticle);

                if (currentParticle.fitness < currentParticle.bestFitness)
                {
                    currentParticle.bestFitness = currentParticle.fitness;

                    for (int i = 0; i < _dimensionNumber; i++)
                    {
                        currentParticle.bestPosition[i] = currentParticle.position[i];
                    }
                }
            }

            return currentParticle;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="currentParticle"></param>
        /// <returns></returns>
        private static Particle UpdatePosition(Particle currentParticle)
        {
            for (int i = 0; i < _dimensionNumber; i++)
            {
                currentParticle.position[i] = (currentParticle.position[i] + currentParticle.velocity[i]);
            }

            return currentParticle;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="particleIndex"></param>
        /// <param name="mbest"></param>
        /// <param name="beta"></param>
        /// <returns></returns>
        private static Particle UpdateQuantumPosition(Particle currentParticle, double mbest, double beta)
        {
            for (int i = 0; i < _dimensionNumber; i++)
            {
                double fi = _random.NextDouble();
                double p = fi * currentParticle.bestPosition[i] + (1 - fi) * 
                    _globalBestPosition[(_globalBestSize - 1 - _dimensionNumber) + i];

                if (fi > 0.5)
                {
                    currentParticle.position[i] = p - beta * 
                        Math.Abs(mbest - currentParticle.position[i]) * (-Math.Log10(fi));
                }
                else
                {
                    currentParticle.position[i] = p + beta * 
                        Math.Abs(mbest - currentParticle.position[i]) * (-Math.Log10(fi));
                }

                if (currentParticle.position[i] < _rangeMin)
                {
                    currentParticle.position[i] = 2 * _rangeMin - currentParticle.position[i];
                }

                if (currentParticle.position[i] > _rangeMax)
                {
                    currentParticle.position[i] = 2 * _rangeMax - currentParticle.position[i];
                }
            }

            return currentParticle;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="particleIndex"></param>
        /// <param name="bestLocalPosition"></param>
        private static Particle UpdateVelocity(Particle currentParticle, double[] bestLocalPosition)
        {
            for (int i = 0; i < _dimensionNumber; i++)
            {
                // particle.velocity[i] = (w * paritcle.velocity[i]) + (constant1 * rand() * (particle.bestPosition[i] - particle.position[i])) + (constant2 * rand() * bestPosition[i] - particle.position[i])
                currentParticle.velocity[i] = (_weight * currentParticle.velocity[i]) +
                    (_particleConstant1 * _random.NextDouble() * (currentParticle.bestPosition[i] -
                    currentParticle.position[i])) + (_particleConstant2 * _random.NextDouble() *
                    (bestLocalPosition[i] - currentParticle.position[i]));
            }

            return currentParticle;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="position"></param>
        /// <returns></returns>
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
        /// 
        /// </summary>
        /// <param name="particleIndex"></param>
        /// <returns></returns>
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
        /// 
        /// </summary>
        /// <returns></returns>
        public static Particle GetInitialPosition(Particle currentParticle)
        {
            for (int i = 0; i < _dimensionNumber; i++)
            {
                // The given values are null.
                if (_particleInitialPosition == null)
                {
                    currentParticle.position[i] = ((_rangeMax - _rangeMin) * _random.NextDouble() + _rangeMin);
                }
                // There is only one given value.
                else if (_particleInitialPositionSize == 1)
                {
                    currentParticle.position[i] = _particleInitialPosition[0];
                }
                // There are as many given value as the position dimensions.
                else if (_particleInitialPositionSize == _dimensionNumber)
                {
                    currentParticle.position[i] = _particleInitialPosition[i];
                }
                // The current position can be set from the given values.
                else if (i < _particleInitialPositionSize)
                {
                    currentParticle.position[i] = _particleInitialPosition[i];
                }
                // The current position can't be set from the given values.
                else
                {
                    currentParticle.position[i] = ((_rangeMax - _rangeMin) * _random.NextDouble() + _rangeMin);
                }
            }

            return currentParticle;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public static Particle GetInitialVelocity(Particle currentParticle)
        {
            for (int i = 0; i < _dimensionNumber; i++)
            {
                // The given values are null.
                if (_particleInitialVelocity == null)
                {
                    currentParticle.velocity[i] = ((_rangeMax - _rangeMin) * _random.NextDouble() + _rangeMin);
                }
                // There is only one given value.
                else if (_particleInitialVelocitySize == 1)
                {
                    currentParticle.velocity[i] = _particleInitialVelocity[0];
                }
                // There are as many given value as the velocity dimensions.
                else if (_particleInitialVelocitySize == _dimensionNumber)
                {
                    currentParticle.velocity[i] = _particleInitialVelocity[i];
                }
                // The current velocity can be set from the given values.
                else if (i < _particleInitialVelocitySize)
                {
                    currentParticle.velocity[i] = _particleInitialVelocity[i];
                }
                // The current velocity can't be set from the given values.
                else
                {
                    currentParticle.velocity[i] = ((_rangeMax - _rangeMin) * _random.NextDouble() + _rangeMin);
                }
            }

            return currentParticle;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="startIndex"></param>
        /// <param name="count"></param>
        /// <returns></returns>
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
        /// 
        /// </summary>
        /// <param name="index"></param>
        /// <param name="indexSize"></param>
        /// <returns></returns>
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
        /// 
        /// </summary>
        /// <returns></returns>
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

        #endregion Methods
    }

    /// <summary>
    /// 
    /// </summary>
    class Particle
    {
        // 
        public double[] position;

        // 
        public double[] velocity;

        // 
        public int index;

        // 
        public double fitness;

        // 
        public double[] bestPosition;

        // 
        public double bestFitness;
    }
}
