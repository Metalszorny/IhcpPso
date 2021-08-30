using System;
using System.Diagnostics;
using System.IO;

namespace IhcpPso2dCSharpSingleThreaded
{
    /// <summary>
    /// 
    /// </summary>
    class Program
    {
        #region Fields

        #region HTC

        // The initial temperature of the object that's cooling down.
        private static double _initialTemperature;

        // 
        private static int _horizontalSplitting;

        // 
        private static int _verticalSplitting;

        // 
        private static int _monitoredIndexRowNumber;

        // 
        private static int _monitoredIndexColumnNumber;

        // The index of the horizontal part that is monitored.
        private static int[] _monitoredIndex;

        // The differencial value of time that is incemented in the cooldown.
        private static double _timeDifference;

        // The number of dimensions of the heat transfer coefficients.
        private static int _dimensionRowNumber;

        // The number of dimensions of the heat transfer coefficients.
        private static int _dimensionColumnNumber;

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
        private static double _rX;

        // 
        private static double _rY;

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

        // 
        private static double _enviromentTemperature;

        // 
        private static bool _partlyCooling;

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
        private static double[] _particlePosition1;

        // 
        private static double[] _particlePosition2;

        // 
        private static double[] _particleVelocity1;

        // 
        private static double[] _particleVelocity2;

        // 
        private static double[] _particleBestPosition1;

        // 
        private static double[] _particleBestPosition2;

        // 
        private static double[] _particleFitness1;

        // 
        private static double[] _particleFitness2;

        // 
        private static double[] _particleBestFitness1;

        // 
        private static double[] _particleBestFitness2;

        // 
        private static int _globalBestSize;

        // 
        private static double[] _globalBestFitness1;

        // 
        private static double[] _globalBestFitness2;

        // 
        private static double[] _globalBestPosition1;

        // 
        private static double[] _globalBestPosition2;

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

        // 
        private static Stopwatch _stopwatch;

        // 
        private static string _exitReason;

        #endregion PSO

        // Pseudo random number generator.
        private static Random _random;

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
                if (_particleBestFitness1[i] < _globalBestFitness1[_globalBestSize - 1] &&
                    _particleBestFitness2[i] < _globalBestFitness2[_globalBestSize - 1])
                {
                    _globalBestFitness1[_globalBestSize - 1] = _particleBestFitness1[i];
                    _globalBestFitness2[_globalBestSize - 1] = _particleBestFitness2[i];

                    for (int j = 0; j < _dimensionRowNumber; j++)
                    {
                        _globalBestPosition1[(_globalBestSize - 1) + j] = _particlePosition1[(i * _dimensionRowNumber) + j];
                        _globalBestPosition2[(_globalBestSize - 1) + j] = _particlePosition2[(i * _dimensionRowNumber) + j];
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
                                string[] values = streamReader.ReadLine().Trim().Replace(".", ",").Split(' ');
                                _monitoredIndexRowNumber = Convert.ToInt32(values[0]);
                                _monitoredIndexColumnNumber = Convert.ToInt32(values[1]);

                                if (_monitoredIndexRowNumber > 0 && _monitoredIndexColumnNumber > 0)
                                {
                                    _monitoredIndex = new int[_monitoredIndexRowNumber * _monitoredIndexColumnNumber];

                                    for (int i = 0; i < _monitoredIndexRowNumber; i++)
                                    {
                                        string[] result = streamReader.ReadLine().Trim().Replace(".", ",").Split(' ');

                                        for (int j = 0; j < _monitoredIndexColumnNumber; j++)
                                        {
                                            _monitoredIndex[(i * _monitoredIndexRowNumber) + j] = Convert.ToInt32(result[j]);
                                        }
                                    }
                                }
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
                            else if (currentLine.ToLower().Contains("enviroment temperature"))
                            {
                                _enviromentTemperature = Convert.ToDouble(streamReader.ReadLine().Trim().Replace(".", ","));
                            }
                            else if (currentLine.ToLower().Contains("shielded sides"))
                            {
                                if (Convert.ToInt32(streamReader.ReadLine().Trim().Replace(".", ",")) == 1)
                                {
                                    _partlyCooling = true;
                                }
                                else
                                {
                                    _partlyCooling = false;
                                }
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
                                string[] result = streamReader.ReadLine().Split(' ');
                                // Number of HTC values.
                                _dimensionRowNumber = Convert.ToInt32(result[0]);
                                // Number of positions needed for a particle + temperature for the htc values.
                                _dimensionColumnNumber = Convert.ToInt32(result[1]);

                                if (_dimensionRowNumber > 0)
                                {
                                    _htcValues = new double[_dimensionRowNumber * _dimensionColumnNumber];

                                    for (int i = 0; i < (_dimensionRowNumber * _dimensionColumnNumber); i += _dimensionColumnNumber)
                                    {
                                        string[] values = (streamReader.ReadLine().Trim().Replace(".", ",")).Split(' ');

                                        for (int j = 0; j < _dimensionColumnNumber; j++)
                                        {
                                            _htcValues[i + j] = Convert.ToDouble(values[j]);
                                        }
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
            _rX = _timeDifference / (_dX * _dX);
            _rY = _timeDifference / (_dY * _dY);
            _stopwatch = new Stopwatch();
            _random = new Random();
            _exitReason = "";
            // Get particle swarm initial values.
            _particlePosition1 = new double[_particleSwarmSize * _dimensionRowNumber];
            _particlePosition2 = new double[_particleSwarmSize * _dimensionRowNumber];
            _particleVelocity1 = new double[_particleSwarmSize * _dimensionRowNumber];
            _particleVelocity2 = new double[_particleSwarmSize * _dimensionRowNumber];
            _particleBestPosition1 = new double[_particleSwarmSize * _dimensionRowNumber];
            _particleBestPosition2 = new double[_particleSwarmSize * _dimensionRowNumber];
            _particleFitness1 = new double[_particleSwarmSize];
            _particleFitness2 = new double[_particleSwarmSize];
            _particleBestFitness1 = new double[_particleSwarmSize];
            _particleBestFitness2 = new double[_particleSwarmSize];
            _globalBestSize = 1;
            _globalBestFitness1 = new double[_globalBestSize];
            _globalBestFitness2 = new double[_globalBestSize];
            _globalBestFitness1[_globalBestSize - 1] = double.MaxValue;
            _globalBestFitness2[_globalBestSize - 1] = double.MaxValue;
            _globalBestPosition1 = new double[_globalBestSize * _dimensionRowNumber];
            _globalBestPosition2 = new double[_globalBestSize * _dimensionRowNumber];
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        private static double[,] SetInitialCurrentTemperature()
        {
            //double[] returnValue = new double[_horizontalSplitting * _verticalSplitting];
            double[,] returnValue = new double[_horizontalSplitting, _verticalSplitting];

            for (int i = 0; i < _horizontalSplitting; i++)
            {
                for (int j = 0; j < _verticalSplitting; j++)
                {
                    //returnValue[(i * _verticalSplitting) + j] = _initialTemperature;
                    returnValue[i, j] = _initialTemperature;
                }
            }

            return returnValue;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        private static double[,] SetInitialPreviousTemperature()
        {
            //double[] returnValue = new double[_horizontalSplitting * _verticalSplitting];
            double[,] returnValue = new double[_horizontalSplitting, _verticalSplitting];

            for (int i = 0; i < _horizontalSplitting; i++)
            {
                for (int j = 0; j < _verticalSplitting; j++)
                {
                    //returnValue[(i * _verticalSplitting) + j] = _initialTemperature;
                    returnValue[i, j] = _initialTemperature;
                }
            }

            return returnValue;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        private static double[,] SetInitialG()
        {
            //double[] returnValue = new double[_horizontalSplitting * _verticalSplitting];
            double[,] returnValue = new double[_horizontalSplitting, _verticalSplitting];

            for (int i = 0; i < _horizontalSplitting; i++)
            {
                for (int j = 0; j < _verticalSplitting; j++)
                {
                    //returnValue[(i * _verticalSplitting) + j] = 0;
                    returnValue[i, j] = 0;
                }
            }

            return returnValue;
        }

        /// <summary>
        /// 
        /// </summary>
        private static void CalculateReferenceCooldownCurve()
        {
            _referenceCooldownCurveSizeRows = (int)(_simulationTime / _timeDifference) + 1;
            _referenceCooldownCurve = new double[_referenceCooldownCurveSizeRows * (_monitoredIndexRowNumber + 1)];
            double[,] currentTemperature = SetInitialCurrentTemperature();
            double[,] previousTemperature = SetInitialPreviousTemperature();
            double[,] g = SetInitialG();

            for (int i = 0; i < (_referenceCooldownCurveSizeRows * (_monitoredIndexRowNumber + 1)); i += (_monitoredIndexRowNumber + 1))
            {
                currentTemperature = CalculateCooldownCurve2D(false, -1, currentTemperature, previousTemperature, g);
                _referenceCooldownCurve[i + 0] = ((i / (_monitoredIndexRowNumber + 1)) * _timeDifference);

                for (int j = 1; j < (_monitoredIndexRowNumber + 1); j++)
                {
                    int index0 = _monitoredIndex[(j - 1) * _monitoredIndexColumnNumber];
                    int index1 = _monitoredIndex[((j - 1) * _monitoredIndexColumnNumber) + 1];
                    int index = (index0 * _horizontalSplitting) + index1;
                    _referenceCooldownCurve[i + j] = currentTemperature[index1, index0];
                }

                //for (int j = 0; j < (_horizontalSplitting * _verticalSplitting); j++)
                //{
                //    previousTemperature[j] = currentTemperature[j];
                //}

                for (int j = 0; j < _horizontalSplitting; j++)
                {
                    for (int k = 0; k < _verticalSplitting; k++)
                    {
                        previousTemperature[j, k] = currentTemperature[j, k];
                    }
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
        private static double[,] CalculateCooldownCurve2D(bool isParticle, int particleIndex, double[,] currentTemperature, double[,] previousTemperature, double[,] g)
        {
            double temperature = 0;
            double heatConductivity = 0;
            double alpha = 0;

            // Top side.
            for (int i = 1; i < _horizontalSplitting - 1; i++)
            {
                //temperature = previousTemperature[(i * _verticalSplitting) + _verticalSplitting - 1];
                //heatConductivity = GetHeatConductivity(temperature);
                //alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                //currentTemperature[(i * _verticalSplitting) + _verticalSplitting - 1] =
                //    2 * _rY * alpha * (previousTemperature[(i * _verticalSplitting) + _verticalSplitting - 2] -
                //    previousTemperature[(i * _verticalSplitting) + _verticalSplitting - 1] + (_dY / heatConductivity) *
                //    GetHeatTransferCoefficient(temperature, isParticle, particleIndex, i, (_verticalSplitting - 1)) *
                //    (_enviromentTemperature - previousTemperature[(i * _verticalSplitting) + _verticalSplitting - 1])) +
                //    _rX * alpha * (previousTemperature[((i - 1) * _verticalSplitting) + _verticalSplitting - 1] - 2 *
                //    previousTemperature[(i * _verticalSplitting) + _verticalSplitting - 1] +
                //    previousTemperature[((i + 1) * _verticalSplitting) + _verticalSplitting - 1]) +
                //    (_timeDifference / _dX) * (1 / (2 * _radius)) * alpha *
                //    Math.Abs(previousTemperature[((i + 1) * _verticalSplitting) + _verticalSplitting - 1] -
                //    previousTemperature[((i - 1) * _verticalSplitting) + _verticalSplitting - 1]) +
                //    (_timeDifference / heatConductivity) * alpha +
                //    previousTemperature[(i * _verticalSplitting) + _verticalSplitting - 1];
                temperature = previousTemperature[i, _verticalSplitting - 1];
                heatConductivity = GetHeatConductivity(temperature);
                alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                currentTemperature[i, _verticalSplitting - 1] =
                    2 * _rY * alpha * (previousTemperature[i, _verticalSplitting - 2] -
                    previousTemperature[i, _verticalSplitting - 1] + (_dY / heatConductivity) *
                    GetHeatTransferCoefficient(temperature, isParticle, particleIndex, i, _verticalSplitting - 1) *
                    (_enviromentTemperature - previousTemperature[i, _verticalSplitting - 1])) +
                    _rX * alpha * (previousTemperature[i - 1, _verticalSplitting - 1] -
                    2 * previousTemperature[i, _verticalSplitting - 1] +
                    previousTemperature[i + 1, _verticalSplitting - 1]) +
                    (_timeDifference / _dX) * (1 / (2 * _radius)) * alpha *
                    Math.Abs(previousTemperature[i + 1, _verticalSplitting - 1] -
                    previousTemperature[i - 1, _verticalSplitting - 1]) +
                    (_timeDifference / heatConductivity) * alpha +
                    previousTemperature[i, _verticalSplitting - 1];
            }

            // Bottom side.
            for (int i = 1; i < _horizontalSplitting - 1; i++)
            {
                //temperature = previousTemperature[i * _verticalSplitting];
                //heatConductivity = GetHeatConductivity(temperature);
                //alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                //currentTemperature[i * _verticalSplitting] =
                //    2 * _rY * alpha * (previousTemperature[(i * _verticalSplitting) + 1] -
                //    previousTemperature[i * _verticalSplitting] + (_dY / heatConductivity) *
                //    GetHeatTransferCoefficient(temperature, isParticle, particleIndex, i, 0) *
                //    (_enviromentTemperature - previousTemperature[i * _verticalSplitting])) + _rX * alpha *
                //    (previousTemperature[(i - 1) * _verticalSplitting] - 2 *
                //    previousTemperature[i * _verticalSplitting] + previousTemperature[(i + 1) * _verticalSplitting]) +
                //    (_timeDifference / _dX) * (1 / (2 * _radius)) * alpha *
                //    Math.Abs(previousTemperature[(i + 1) * _verticalSplitting] - previousTemperature[(i - 1) * _verticalSplitting]) +
                //    (_timeDifference / heatConductivity) * alpha + previousTemperature[i * _verticalSplitting];
                temperature = previousTemperature[i, 0];
                heatConductivity = GetHeatConductivity(temperature);
                alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                currentTemperature[i, 0] = 
                    2 * _rY * alpha * (previousTemperature[i, 1] - previousTemperature[i, 0] + (_dY / heatConductivity) * 
                    GetHeatTransferCoefficient(temperature, isParticle, particleIndex, i, 0) * 
                    (_enviromentTemperature - previousTemperature[i, 0])) + _rX * alpha *
                    (previousTemperature[i - 1, 0] - 2 * previousTemperature[i, 0] + previousTemperature[i + 1, 0]) +
                    (_timeDifference / _dX) * (1 / (2 * _radius)) * alpha *
                    Math.Abs(previousTemperature[i - 1, 0] - previousTemperature[i, 1]) +
                    (_timeDifference / heatConductivity) * alpha + previousTemperature[i, 0];
            }

            if (!_partlyCooling)
            {
                // Inner side.
                for (int i = 1; i < _verticalSplitting - 1; i++)
                {
                    //temperature = previousTemperature[i];
                    //heatConductivity = GetHeatConductivity(temperature);
                    //alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                    //currentTemperature[i] =
                    //    previousTemperature[i] + 4 * _rX * alpha * (previousTemperature[(1 * _verticalSplitting) + i] -
                    //    previousTemperature[i]) + _rY * alpha * (previousTemperature[i - 1] - 2 * previousTemperature[i] +
                    //    previousTemperature[i + 1]) + (_timeDifference / heatConductivity) * alpha;
                    temperature = previousTemperature[0, i];
                    heatConductivity = GetHeatConductivity(temperature);
                    alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                    currentTemperature[0, i] = 
                        previousTemperature[0, i] + 4 * _rX * alpha * (previousTemperature[1, i] - 
                        previousTemperature[0, i]) + _rY * alpha * (previousTemperature[0, i - 1] - 
                        2 * previousTemperature[0, i] + previousTemperature[0, i + 1]) +
                        (_timeDifference / heatConductivity) * alpha;
                }

                // Outer side.
                for (int i = 1; i < _verticalSplitting - 1; i++)
                {
                    //temperature = previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + i];
                    //heatConductivity = GetHeatConductivity(temperature);
                    //alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                    //currentTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + i] =
                    //    _rY * alpha * (previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (i - 1)] - 2 *
                    //    previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + i] +
                    //    previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (i + 1)]) + 2 * _rX * alpha *
                    //    (previousTemperature[((_horizontalSplitting - 2) * _verticalSplitting) + i] -
                    //    previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + i] + (_dX / heatConductivity) *
                    //    GetHeatTransferCoefficient(temperature, isParticle, particleIndex, (_horizontalSplitting - 1), i) *
                    //    (_enviromentTemperature - previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + i])) +
                    //    (_timeDifference / heatConductivity) * (1 / _radius) * alpha *
                    //    GetHeatTransferCoefficient(temperature, isParticle, particleIndex, (_horizontalSplitting - 1), i) *
                    //    (_enviromentTemperature - previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + i]) +
                    //    (_timeDifference / heatConductivity) * alpha + previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + i];
                    temperature = previousTemperature[(_horizontalSplitting - 1), i];
                    heatConductivity = GetHeatConductivity(temperature);
                    alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                    currentTemperature[_horizontalSplitting - 1, i] = 
                        _rY * alpha * (previousTemperature[_horizontalSplitting - 1, i - 1] - 2 * 
                        previousTemperature[_horizontalSplitting - 1, i] + 
                        previousTemperature[_horizontalSplitting - 1, i + 1]) + 2 * _rX * alpha *
                        (previousTemperature[_horizontalSplitting - 1 - 1, i] - 
                        previousTemperature[_horizontalSplitting - 1, i] + (_dX / heatConductivity) *
                        GetHeatTransferCoefficient(temperature, isParticle, particleIndex, _horizontalSplitting - 1, i) * 
                        (_enviromentTemperature - previousTemperature[_horizontalSplitting - 1, i])) +
                        (_timeDifference / heatConductivity) * (1 / _radius) * alpha *
                        GetHeatTransferCoefficient(temperature, isParticle, particleIndex, _horizontalSplitting - 1, i) * 
                        (_enviromentTemperature - previousTemperature[_horizontalSplitting - 1, i]) +
                        (_timeDifference / heatConductivity) * alpha + previousTemperature[_horizontalSplitting - 1, i];
                }
            }
            else
            {
                int wholeNumber = ((_verticalSplitting - 1) / 2);

                // Inner side.
                for (int i = 1; i <= wholeNumber; i++)
                {
                    //temperature = previousTemperature[i];
                    //heatConductivity = GetHeatConductivity(temperature);
                    //alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                    //currentTemperature[i] =
                    //    previousTemperature[i] + 4 * _rX * alpha * (previousTemperature[(1 * _verticalSplitting) + i - 1] -
                    //    previousTemperature[i]) + _rY * alpha * (previousTemperature[i - 1] - 2 * previousTemperature[i] +
                    //    previousTemperature[i + 1]) + (_timeDifference / heatConductivity) * alpha;
                    //temperature = previousTemperature[(_verticalSplitting - 1) - i];
                    //heatConductivity = GetHeatConductivity(temperature);
                    //alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                    //currentTemperature[(_verticalSplitting - 1) - i] =
                    //    previousTemperature[(_verticalSplitting - 1) - i] + 4 * _rX * alpha *
                    //    (previousTemperature[(1 * _verticalSplitting) - 1 + (_verticalSplitting - 1) - i] -
                    //    previousTemperature[(_verticalSplitting - 1) - i]) + _rY * alpha *
                    //    (previousTemperature[(_verticalSplitting - 1) - i - 1] - 2 * previousTemperature[(_verticalSplitting - 1) - i] +
                    //    previousTemperature[(_verticalSplitting - 1) - i + 1]) + (_timeDifference / heatConductivity) * alpha;
                    temperature = previousTemperature[0, i];
                    heatConductivity = GetHeatConductivity(temperature);
                    alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                    currentTemperature[0, i] = 
                        previousTemperature[0, i] + 4 * _rX * alpha * (previousTemperature[1, i] - 
                        previousTemperature[0, i]) + _rY * alpha * (previousTemperature[0, i - 1] - 2 * previousTemperature[0, i] + 
                        previousTemperature[0, i + 1]) + (_timeDifference / heatConductivity) * alpha;
                    temperature = previousTemperature[0, (_verticalSplitting - 1) - i];
                    heatConductivity = GetHeatConductivity(temperature);
                    alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                    currentTemperature[0, (_verticalSplitting - 1) - i] = 
                        previousTemperature[0, (_verticalSplitting - 1) - i] + 4 * _rX * alpha *
                        (previousTemperature[1, (_verticalSplitting - 1) - i] - 
                        previousTemperature[0, (_verticalSplitting - 1) - i]) + _rY * alpha * 
                        (previousTemperature[0, (_verticalSplitting - 1) - i - 1] - 
                        2 * previousTemperature[0, (_verticalSplitting - 1) - i] +
                        previousTemperature[0, (_verticalSplitting - 1) - i + 1]) + 
                        (_timeDifference / heatConductivity) * alpha;
                }

                // Outer side.
                for (int i = 1; i <= wholeNumber; i++)
                {
                    //temperature = previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + i];
                    //heatConductivity = GetHeatConductivity(temperature);
                    //alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                    //currentTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + i] =
                    //    _rY * alpha * (previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (i - 1)] - 2 *
                    //    previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + i] +
                    //    previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (i + 1)]) + 2 * _rX * alpha *
                    //    (previousTemperature[((_horizontalSplitting - 2) * _verticalSplitting) + i] -
                    //    previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + i] + (_dX / heatConductivity) *
                    //    1 * (_enviromentTemperature - previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + i])) +
                    //    (_timeDifference / heatConductivity) * (1 / _radius) * alpha *
                    //    1 * (_enviromentTemperature - previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + i]) +
                    //    (_timeDifference / heatConductivity) * alpha + previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + i];
                    //temperature = previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (_verticalSplitting - 1) - i];
                    //heatConductivity = GetHeatConductivity(temperature);
                    //alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                    //currentTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (_verticalSplitting - 1) - i] =
                    //    _rY * alpha * (previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (_verticalSplitting - 1) - i - 1] - 2 *
                    //    previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (_verticalSplitting - 1) - i] +
                    //    previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (_verticalSplitting - 1) - i + 1]) + 2 * _rX * alpha *
                    //    (previousTemperature[((_horizontalSplitting - 2) * _verticalSplitting) + (_verticalSplitting - 1) - i] -
                    //    previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (_verticalSplitting - 1) - i] + (_dX / heatConductivity) *
                    //    1 * (_enviromentTemperature - previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (_verticalSplitting - 1) - i])) +
                    //    (_timeDifference / heatConductivity) * (1 / _radius) * alpha *
                    //    1 * (_enviromentTemperature - previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (_verticalSplitting - 1) - i]) +
                    //    (_timeDifference / heatConductivity) * alpha + previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (_verticalSplitting - 1) - i];
                    temperature = previousTemperature[(_horizontalSplitting - 1), i];
                    heatConductivity = GetHeatConductivity(temperature);
                    alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                    currentTemperature[_horizontalSplitting - 1, i] = 
                        _rY * alpha * (previousTemperature[_horizontalSplitting - 1, i - 1] - 2 * 
                        previousTemperature[_horizontalSplitting - 1, i] + 
                        previousTemperature[_horizontalSplitting - 1, i + 1]) + 2 * _rX * alpha *
                        (previousTemperature[_horizontalSplitting - 2, i] - 
                        previousTemperature[_horizontalSplitting - 1, i] + (_dX / heatConductivity) *
                        1 * (_enviromentTemperature - previousTemperature[_horizontalSplitting - 1, i])) +
                        (_timeDifference / heatConductivity) * (1 / _radius) * alpha *
                        1 * (_enviromentTemperature - previousTemperature[_horizontalSplitting - 1, i]) +
                        (_timeDifference / heatConductivity) * alpha + previousTemperature[_horizontalSplitting - 1, i];
                    temperature = previousTemperature[(_horizontalSplitting - 1), (_verticalSplitting - 1) - i];
                    heatConductivity = GetHeatConductivity(temperature);
                    alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                    currentTemperature[_horizontalSplitting - 1, (_verticalSplitting - 1) - i] = 
                        _rY * alpha * (previousTemperature[_horizontalSplitting - 1, (_verticalSplitting - 1) - i - 1] - 2 *
                        previousTemperature[_horizontalSplitting - 1, (_verticalSplitting - 1) - i] + 
                        previousTemperature[_horizontalSplitting - 1, (_verticalSplitting - 1) - i + 1]) + 2 * _rX * alpha *
                        (previousTemperature[_horizontalSplitting - 2, (_verticalSplitting - 1) - i] - 
                        previousTemperature[_horizontalSplitting - 1, (_verticalSplitting - 1) - i] + (_dX / heatConductivity) *
                        1 * (_enviromentTemperature - previousTemperature[_horizontalSplitting - 1, (_verticalSplitting - 1) - i])) +
                        (_timeDifference / heatConductivity) * (1 / _radius) * alpha *
                        1 * (_enviromentTemperature - previousTemperature[_horizontalSplitting - 1, (_verticalSplitting - 1) - i]) +
                        (_timeDifference / heatConductivity) * alpha + previousTemperature[_horizontalSplitting - 1, (_verticalSplitting - 1) - i];
                }

                if ((_verticalSplitting - 1) % 2 != 0)
                {
                    int middle = ((_verticalSplitting - 1) / 2) + 1;

                    // Inner side.
                    //temperature = previousTemperature[middle];
                    //heatConductivity = GetHeatConductivity(temperature);
                    //alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                    //currentTemperature[middle] =
                    //    previousTemperature[middle] + 4 * _rX * alpha * (previousTemperature[(1 * _verticalSplitting) + middle - 1] -
                    //    previousTemperature[middle]) + _rY * alpha * (previousTemperature[middle - 1] - 2 *
                    //    previousTemperature[middle] + previousTemperature[middle + 1]) + (_timeDifference / heatConductivity) * alpha;
                    temperature = previousTemperature[0, middle];
                    heatConductivity = GetHeatConductivity(temperature);
                    alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                    currentTemperature[0, middle] = 
                        previousTemperature[0, middle] + 4 * _rX * alpha * (previousTemperature[1, middle] - 
                        previousTemperature[0, middle]) + _rY * alpha * (previousTemperature[0, middle - 1] - 2 * 
                        previousTemperature[0, middle] + previousTemperature[0, middle + 1]) + (_timeDifference / heatConductivity) * alpha;

                    // Outer side.
                    //temperature = previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + middle - 1];
                    //heatConductivity = GetHeatConductivity(temperature);
                    //alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                    //currentTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + middle - 1] =
                    //    _rY * alpha * (previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (middle - 2)] - 2 *
                    //    previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + middle - 1] +
                    //    previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + middle]) + 2 * _rX * alpha *
                    //    (previousTemperature[((_horizontalSplitting - 2) * _verticalSplitting) + middle - 1] -
                    //    previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + middle - 1] + (_dX / heatConductivity) *
                    //    1 * (_enviromentTemperature - previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + middle - 1])) +
                    //    (_timeDifference / heatConductivity) * (1 / _radius) * alpha *
                    //    1 * (_enviromentTemperature - previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + middle - 1]) +
                    //    (_timeDifference / heatConductivity) * alpha + previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + middle - 1];
                    temperature = previousTemperature[(_horizontalSplitting - 1), middle - 1];
                    heatConductivity = GetHeatConductivity(temperature);
                    alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                    currentTemperature[_horizontalSplitting - 1, middle] = 
                        _rY * alpha * (previousTemperature[_horizontalSplitting - 1, middle - 1] - 2 * 
                        previousTemperature[_horizontalSplitting - 1, middle] + 
                        previousTemperature[_horizontalSplitting - 1, middle + 1]) + 2 * _rX * alpha *
                        (previousTemperature[_horizontalSplitting - 2, middle] - 
                        previousTemperature[_horizontalSplitting - 1, middle] + (_dX / heatConductivity) *
                        1 * (_enviromentTemperature - previousTemperature[_horizontalSplitting - 1, middle])) +
                        (_timeDifference / heatConductivity) * (1 / _radius) * alpha *
                        1 * (_enviromentTemperature - previousTemperature[_horizontalSplitting - 1, middle]) +
                        (_timeDifference / heatConductivity) * alpha + previousTemperature[_horizontalSplitting - 1, middle];
                }
            }

            // Corners.
            // Inner bottom corner.
            //temperature = previousTemperature[0];
            //heatConductivity = GetHeatConductivity(temperature);
            //alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
            //currentTemperature[0] =
            //    previousTemperature[0] + 4 * _rX * alpha * (previousTemperature[1 * _verticalSplitting - 1] -
            //    previousTemperature[0]) + 2 * _rY * alpha * (previousTemperature[1] - previousTemperature[0] +
            //    (_dY / heatConductivity) * GetHeatTransferCoefficient(temperature, isParticle, particleIndex, 0, 0) *
            //    (_enviromentTemperature - previousTemperature[0])) + (_timeDifference / heatConductivity) * alpha;
            temperature = previousTemperature[0, 0];
            heatConductivity = GetHeatConductivity(temperature);
            alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
            currentTemperature[0, 0] = 
                previousTemperature[0, 0] + 4 * _rX * alpha * (previousTemperature[1, 0] - 
                previousTemperature[0, 0]) + 2 * _rY * alpha * (previousTemperature[0, 1] - previousTemperature[0, 0] + 
                (_dY / heatConductivity) * GetHeatTransferCoefficient(temperature, isParticle, particleIndex, 0, 0) *
                (_enviromentTemperature - previousTemperature[0, 0])) + (_timeDifference / heatConductivity) * alpha;
            // Inner top corner.
            //temperature = previousTemperature[_verticalSplitting - 1];
            //heatConductivity = GetHeatConductivity(temperature);
            //alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
            //currentTemperature[_verticalSplitting - 1] =
            //    previousTemperature[_verticalSplitting - 1] + 4 * _rX * alpha * (previousTemperature[(1 * _verticalSplitting) +
            //    (_verticalSplitting - 2)] - previousTemperature[_verticalSplitting - 1]) + 2 * _rY * alpha *
            //    (previousTemperature[_verticalSplitting - 2] - previousTemperature[_verticalSplitting - 1] + (_dY / heatConductivity) *
            //    GetHeatTransferCoefficient(temperature, isParticle, particleIndex, 0, _verticalSplitting - 1) *
            //    (_enviromentTemperature - previousTemperature[_verticalSplitting - 1])) + (_timeDifference / heatConductivity) * alpha;
            temperature = previousTemperature[0, _verticalSplitting - 1];
            heatConductivity = GetHeatConductivity(temperature);
            alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
            currentTemperature[0, _verticalSplitting - 1] = 
                previousTemperature[0, _verticalSplitting - 1] + 4 * _rX * alpha * (previousTemperature[1, _verticalSplitting - 1] - 
                previousTemperature[0, _verticalSplitting - 1]) + 2 * _rY * alpha * 
                (previousTemperature[0, _verticalSplitting - 2] - previousTemperature[0, _verticalSplitting - 1] + (_dY / heatConductivity) *
                GetHeatTransferCoefficient(temperature, isParticle, particleIndex, 0, _verticalSplitting - 1) * 
                (_enviromentTemperature - previousTemperature[0, _verticalSplitting - 1])) + (_timeDifference / heatConductivity) * alpha;
            // Outer bottom corner.
            //temperature = previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting)];
            //heatConductivity = GetHeatConductivity(temperature);
            //alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
            //currentTemperature[((_horizontalSplitting - 1) * _verticalSplitting)] =
            //    2 * _rY * alpha * ((_dY / heatConductivity) * GetHeatTransferCoefficient(temperature, isParticle, particleIndex, (_horizontalSplitting - 1), 0) *
            //    (_enviromentTemperature - previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting)]) -
            //    previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting)] + previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + 1]) +
            //    2 * _rX * alpha * ((_dX / heatConductivity) * GetHeatTransferCoefficient(temperature, isParticle, particleIndex, (_horizontalSplitting - 1), 0) *
            //    (_enviromentTemperature - previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting)]) -
            //    previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting)] + previousTemperature[((_horizontalSplitting - 2) * _verticalSplitting)]) +
            //    (_timeDifference / heatConductivity) * (1 / _radius) * alpha *
            //    GetHeatTransferCoefficient(temperature, isParticle, particleIndex, (_horizontalSplitting - 1), 0) *
            //    (_enviromentTemperature - previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting)]) +
            //    (_timeDifference / heatConductivity) * alpha + previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting)];
            temperature = previousTemperature[(_horizontalSplitting - 1), 0];
            heatConductivity = GetHeatConductivity(temperature);
            alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
            currentTemperature[_horizontalSplitting - 1, 0] = 
                2 * _rY * alpha * ((_dY / heatConductivity) * GetHeatTransferCoefficient(temperature, isParticle, particleIndex, _horizontalSplitting - 1, 0) *
                (_enviromentTemperature - previousTemperature[_horizontalSplitting - 1, 0]) - 
                previousTemperature[_horizontalSplitting - 1, 0] + previousTemperature[_horizontalSplitting - 1, 1]) +
                2 * _rX * alpha * ((_dX / heatConductivity) * GetHeatTransferCoefficient(temperature, isParticle, particleIndex, _horizontalSplitting - 1, 0) *
                (_enviromentTemperature - previousTemperature[_horizontalSplitting - 1, 0]) - 
                previousTemperature[_horizontalSplitting - 1, 0] + previousTemperature[_horizontalSplitting - 2, 0]) +
                (_timeDifference / heatConductivity) * (1 / _radius) * alpha * 
                GetHeatTransferCoefficient(temperature, isParticle, particleIndex, _horizontalSplitting - 1, 0) *
                (_enviromentTemperature - previousTemperature[_horizontalSplitting - 1, 0]) + 
                (_timeDifference / heatConductivity) * alpha + previousTemperature[_horizontalSplitting - 1, 0];
            // Outer top corner.
            //temperature = previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (_verticalSplitting - 1)];
            //heatConductivity = GetHeatConductivity(temperature);
            //alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
            //currentTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (_verticalSplitting - 1)] =
            //    2 * _rY * alpha * (previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (_verticalSplitting - 1) - 1] -
            //    previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (_verticalSplitting - 1)] + (_dY / heatConductivity) *
            //    GetHeatTransferCoefficient(temperature, isParticle, particleIndex, (_horizontalSplitting - 1), (_verticalSplitting - 1)) *
            //    (_enviromentTemperature - previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (_verticalSplitting - 1)])) + 2 * _rX * alpha *
            //    (previousTemperature[((_horizontalSplitting - 2) * _verticalSplitting) + (_verticalSplitting - 1)] -
            //    previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (_verticalSplitting - 1)] + (_dX / heatConductivity) *
            //    GetHeatTransferCoefficient(temperature, isParticle, particleIndex, (_horizontalSplitting - 1), (_verticalSplitting - 1)) *
            //    (_enviromentTemperature - previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (_verticalSplitting - 1)])) + (_timeDifference / heatConductivity) * (1 / _radius) * alpha *
            //    GetHeatTransferCoefficient(temperature, isParticle, particleIndex, (_horizontalSplitting - 1), (_verticalSplitting - 1)) *
            //    (_enviromentTemperature - previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (_verticalSplitting - 1)]) +
            //    (_timeDifference / heatConductivity) * alpha + previousTemperature[((_horizontalSplitting - 1) * _verticalSplitting) + (_verticalSplitting - 1)];
            temperature = previousTemperature[(_horizontalSplitting - 1), (_verticalSplitting - 1)];
            heatConductivity = GetHeatConductivity(temperature);
            alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
            currentTemperature[_horizontalSplitting - 1, _verticalSplitting - 1] = 
                2 * _rY * alpha * (previousTemperature[_horizontalSplitting - 1, _verticalSplitting - 2] -
                previousTemperature[_horizontalSplitting - 1, _verticalSplitting - 1] + (_dY / heatConductivity) *
                GetHeatTransferCoefficient(temperature, isParticle, particleIndex, _horizontalSplitting - 1, _verticalSplitting - 1) *
                (_enviromentTemperature - previousTemperature[_horizontalSplitting - 1, _verticalSplitting - 1])) + 2 * _rX * alpha * 
                (previousTemperature[_horizontalSplitting - 2, _verticalSplitting - 1] -
                previousTemperature[_horizontalSplitting - 1, _verticalSplitting - 1] + (_dX / heatConductivity) *
                GetHeatTransferCoefficient(temperature, isParticle, particleIndex, _horizontalSplitting - 1, _verticalSplitting - 1) *
                (_enviromentTemperature - previousTemperature[_horizontalSplitting - 1, _verticalSplitting - 1])) + (_timeDifference / heatConductivity) * (1 / _radius) * alpha *
                GetHeatTransferCoefficient(temperature, isParticle, particleIndex, _horizontalSplitting - 1, _verticalSplitting - 1) *
                (_enviromentTemperature - previousTemperature[_horizontalSplitting - 1, _verticalSplitting - 1]) + 
                (_timeDifference / heatConductivity) * alpha + previousTemperature[_horizontalSplitting - 1, _verticalSplitting - 1];

            // Inside.
            for (int i = 1; i < _horizontalSplitting - 1; i++)
            {
                for (int j = 1; j < _verticalSplitting - 1; j++)
                {
                    //temperature = previousTemperature[(i * _verticalSplitting) + j];
                    //heatConductivity = GetHeatConductivity(temperature);
                    //alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                    //currentTemperature[(i * _verticalSplitting) + j] =
                    //    _rY * alpha * (previousTemperature[(i * _verticalSplitting) + j - 1] - 2 *
                    //    previousTemperature[(i * _verticalSplitting) + j] + previousTemperature[(i * _verticalSplitting) + j]) +
                    //    _rX * alpha * (previousTemperature[((i - 1) * _verticalSplitting) + j] - 2 * previousTemperature[(i * _verticalSplitting) + j] +
                    //    previousTemperature[((i + 1) * _verticalSplitting) + j - 1]) + (_timeDifference / _dX) * (1 / (2 * _radius)) * alpha *
                    //    Math.Abs(previousTemperature[((i + 1) * _verticalSplitting) + j] - previousTemperature[((i - 1) * _verticalSplitting) + j]) +
                    //    (_timeDifference / heatConductivity) * alpha + previousTemperature[(i * _verticalSplitting) + j];
                    temperature = previousTemperature[i, j];
                    heatConductivity = GetHeatConductivity(temperature);
                    alpha = GetAlpha(heatConductivity, GetSpecificHeat(temperature));
                    currentTemperature[i, j] = 
                        _rY * alpha * (previousTemperature[i, j - 1] - 2 * 
                        previousTemperature[i, j] + previousTemperature[i, j + 1]) + 
                        _rX * alpha * (previousTemperature[i - 1, j] - 2 * previousTemperature[i, j] + 
                        previousTemperature[i + 1, j]) + (_timeDifference / _dX) * (1 / (2 * _radius)) * alpha * 
                        Math.Abs(previousTemperature[i + 1, j] - previousTemperature[i - 1, j]) + 
                        (_timeDifference / heatConductivity) * alpha + previousTemperature[i, j];
                }
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
        private static double GetHeatTransferCoefficient(double temperature, bool isParticle, int particleIndex, int x, int y)
        {
            if (!isParticle)
            {
                if (_dimensionRowNumber > 0)
                {
                    double temperature0 = 0;
                    double temperature1 = 0;
                    double[] heatTransferCoefficient0 = new double[_dimensionColumnNumber - 1];
                    double[] heatTransferCoefficient1 = new double[_dimensionColumnNumber - 1];
                    int i = 0;

                    // Fill with 0-s.
                    for (int j = 1; j < _dimensionColumnNumber; j++)
                    {
                        heatTransferCoefficient0[j - 1] = 0;
                        heatTransferCoefficient1[j - 1] = 0;
                    }

                    // Get smaller htc row.
                    while ((i < (_dimensionRowNumber * _dimensionColumnNumber)) &&
                        (_htcValues[i + 0] <= temperature))
                    {
                        temperature0 = _htcValues[i + 0];

                        for (int j = 1; j < _dimensionColumnNumber; j++)
                        {
                            heatTransferCoefficient0[j - 1] = _htcValues[i + j];
                        }

                        i += _dimensionColumnNumber;
                    }

                    // Larger htc row is awailable.
                    if (i < (_dimensionRowNumber * _dimensionColumnNumber))
                    {
                        temperature1 = _htcValues[i + 0];

                        for (int j = 1; j < _dimensionColumnNumber; j++)
                        {
                            heatTransferCoefficient1[j - 1] = _htcValues[i + j];
                        }

                        if (i == 0)
                        {
                            temperature0 = 0;

                            for (int j = 1; j < _dimensionColumnNumber; j++)
                            {
                                heatTransferCoefficient0[j - 1] = 0;
                            }
                        }

                        double htc0 = ((heatTransferCoefficient1[0] - heatTransferCoefficient0[0]) / (temperature1 - temperature0) * temperature -
                            ((heatTransferCoefficient1[0] - heatTransferCoefficient0[0]) / (temperature1 - temperature0) * temperature0 - heatTransferCoefficient0[0]));
                        double htc1 = ((heatTransferCoefficient1[1] - heatTransferCoefficient0[1]) / (temperature1 - temperature0) * temperature -
                            ((heatTransferCoefficient1[1] - heatTransferCoefficient0[1]) / (temperature1 - temperature0) * temperature0 - heatTransferCoefficient0[1]));
                        double result = (htc0 < htc1 ? htc0 : htc1) + (Math.Abs(htc0 - htc1) * ((y + 1.0) / _verticalSplitting));
                        
                        return result;
                    }
                    // Larger htc row is not awailable.
                    else
                    {
                        double htc0 = _htcValues[((_dimensionRowNumber - 1) * _dimensionColumnNumber) + 1];
                        double htc1 = _htcValues[((_dimensionRowNumber - 1) * _dimensionColumnNumber) + 2];
                        double result = (htc0 < htc1 ? htc0 : htc1) + (Math.Abs(htc0 - htc1) * ((y + 1.0) / _verticalSplitting));
                        
                        return result;
                    }
                }

                return 0;
            }
            else
            {
                if (_dimensionRowNumber > 0)
                {
                    double temperature0 = 0;
                    double temperature1 = 0;
                    double[] heatTransferCoefficient0 = new double[_dimensionColumnNumber - 1];
                    double[] heatTransferCoefficient1 = new double[_dimensionColumnNumber - 1];
                    int i = 0;

                    // Fill with 0-s.
                    for (int j = 1; j < _dimensionColumnNumber; j++)
                    {
                        heatTransferCoefficient0[j - 1] = 0;
                        heatTransferCoefficient1[j - 1] = 0;
                    }

                    while ((i < (_dimensionRowNumber * _dimensionColumnNumber)) &&
                        (_htcValues[i + 0] <= temperature))
                    {
                        heatTransferCoefficient0[0] = _particlePosition1[(particleIndex * _dimensionRowNumber) + (i / _dimensionColumnNumber)];
                        heatTransferCoefficient0[1] = _particlePosition2[(particleIndex * _dimensionRowNumber) + (i / _dimensionColumnNumber)];
                        temperature0 = _htcValues[i + 0];
                        i += _dimensionColumnNumber;
                    }

                    if (i < (_dimensionRowNumber * _dimensionColumnNumber))
                    {
                        heatTransferCoefficient1[0] = _particlePosition1[(particleIndex * _dimensionRowNumber) + (i / _dimensionColumnNumber)];
                        heatTransferCoefficient1[1] = _particlePosition2[(particleIndex * _dimensionRowNumber) + (i / _dimensionColumnNumber)];
                        temperature1 = _htcValues[i + 0];

                        if (i == 0)
                        {
                            temperature0 = 0;

                            for (int j = 1; j < _dimensionColumnNumber; j++)
                            {
                                heatTransferCoefficient0[j - 1] = 0;
                            }
                        }
                        
                        double r1 = ((heatTransferCoefficient1[0] - heatTransferCoefficient0[0]) / (temperature1 - temperature0) * temperature -
                        ((heatTransferCoefficient1[0] - heatTransferCoefficient0[0]) / (temperature1 - temperature0) * temperature0 - heatTransferCoefficient0[0]));
                        double r2 = ((heatTransferCoefficient1[1] - heatTransferCoefficient0[1]) / (temperature1 - temperature0) * temperature -
                        ((heatTransferCoefficient1[1] - heatTransferCoefficient0[1]) / (temperature1 - temperature0) * temperature0 - heatTransferCoefficient0[1]));
                        double result = (r1 < r2 ? r1 : r2) + (Math.Abs(r1 - r2) * ((y + 1.0) / _verticalSplitting));

                        return result;
                    }
                    else
                    {
                        double htc0 = _particlePosition1[(particleIndex * _dimensionRowNumber) + (_dimensionRowNumber - 1)];
                        double htc1 = _particlePosition2[(particleIndex * _dimensionRowNumber) + (_dimensionRowNumber - 1)];
                        double result = (htc0 < htc1 ? htc0 : htc1) + (Math.Abs(htc0 - htc1) * ((y + 1.0) / _verticalSplitting));

                        return result;
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

                        for (int j = 0; j < _monitoredIndexRowNumber; j++)
                        {
                            streamWriter.Write("temperature" + (j + 1) + " ");
                        }

                        streamWriter.Write("\r\n");

                        for (int i = 0; i < (_referenceCooldownCurveSizeRows * (_monitoredIndexRowNumber + 1)); i += (_monitoredIndexRowNumber + 1))
                        {
                            for (int j = 0; j < (_monitoredIndexRowNumber + 1); j++)
                            {
                                streamWriter.Write(_referenceCooldownCurve[i + j].ToString("000.000").Replace(",", ".") + " ");
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
        private static void WritePsoGlobalBestLogToFile()
        {
            if (_globalBestFitness1 != null && _globalBestPosition1 != null &&
                _globalBestFitness2 != null && _globalBestPosition2 != null &&
                _globalBestSize > 0)
            {
                using (Stream stream = File.Open("ParticleSwarmOptimizationGlobalBestLog.txt", FileMode.Create, FileAccess.Write))
                {
                    using (StreamWriter streamWriter = new StreamWriter(stream))
                    {
                        streamWriter.Write("fitness ");

                        for (int i = 0; i < _dimensionRowNumber; i++)
                        {
                            streamWriter.Write("htc" + (i + 1) + " ");
                        }

                        streamWriter.Write("fitness ");

                        for (int i = 0; i < _dimensionRowNumber; i++)
                        {
                            streamWriter.Write("htc" + (i + 1) + " ");
                        }

                        streamWriter.Write("\r\n");

                        for (int i = 0; i < _globalBestSize; i++)
                        {
                            streamWriter.Write(_globalBestFitness1[i].ToString("#.000").Replace(",", ".") + " ");

                            for (int j = 0; j < _dimensionRowNumber; j++)
                            {
                                streamWriter.Write(_globalBestPosition1[(i * _dimensionRowNumber) + j].ToString("00000.000").Replace(",", ".") + " ");
                            }

                            streamWriter.Write(_globalBestFitness2[i].ToString("#.000").Replace(",", ".") + " ");

                            for (int j = 0; j < _dimensionRowNumber; j++)
                            {
                                streamWriter.Write(_globalBestPosition2[(i * _dimensionRowNumber) + j].ToString("00000.000").Replace(",", ".") + " ");
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
        /// <param name="particleIndex"></param>
        private static void InitializeParticle(int particleIndex)
        {
            GetInitialPosition(particleIndex);
            GetInitialVelocity(particleIndex);
            double[] fitnessValues = ObjectiveFunction(particleIndex);
            _particleFitness1[particleIndex] = fitnessValues[0];
            _particleFitness2[particleIndex] = fitnessValues[1];
            _particleBestFitness1[particleIndex] = _particleFitness1[particleIndex];
            _particleBestFitness2[particleIndex] = _particleFitness2[particleIndex];

            for (int j = 0; j < _dimensionRowNumber; j++)
            {
                _particleBestPosition1[(particleIndex * _dimensionRowNumber) + j] = _particlePosition1[(particleIndex * _dimensionRowNumber) + j];
                _particleBestPosition2[(particleIndex * _dimensionRowNumber) + j] = _particlePosition2[(particleIndex * _dimensionRowNumber) + j];
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="particleIndex"></param>
        /// <returns></returns>
        public static double[] ObjectiveFunction(int particleIndex)
        {
            int calculatedCooldownCurveSizeRows = _referenceCooldownCurveSizeRows;
            double[] calculatedCooldownCurve = new double[calculatedCooldownCurveSizeRows * (_monitoredIndexRowNumber + 1)];
            double[,] currentTemperature = SetInitialCurrentTemperature();
            double[,] previousTemperature = SetInitialPreviousTemperature();
            double[,] g = SetInitialG();

            for (int i = 0; i < (calculatedCooldownCurveSizeRows * (_monitoredIndexRowNumber + 1)); i += (_monitoredIndexRowNumber + 1))
            {
                currentTemperature = CalculateCooldownCurve2D(true, particleIndex, currentTemperature, previousTemperature, g);
                calculatedCooldownCurve[i + 0] = ((i / (_monitoredIndexRowNumber + 1)) * _timeDifference);

                for (int j = 1; j < (_monitoredIndexRowNumber + 1); j++)
                {
                    int index0 = _monitoredIndex[(j - 1) * _monitoredIndexColumnNumber];
                    int index1 = _monitoredIndex[((j - 1) * _monitoredIndexColumnNumber) + 1];
                    int index = (index0 * _horizontalSplitting) + index1;
                    calculatedCooldownCurve[i + j] = currentTemperature[index1, index0];
                    //calculatedCooldownCurve[i + j] = currentTemperature[(_monitoredIndex[j * _monitoredIndexColumnNumber] * _verticalSplitting) + _monitoredIndex[(j * _monitoredIndexColumnNumber) + 1]];
                    //calculatedCooldownCurve[i + j] = currentTemperature[_monitoredIndex[(j - 1) * _monitoredIndexColumnNumber], _monitoredIndex[((j - 1) * _monitoredIndexColumnNumber) + 1]];
                }

                //for (int j = 0; j < (_horizontalSplitting * _verticalSplitting); j++)
                //{
                //    previousTemperature[j] = currentTemperature[j];
                //}

                for (int j = 0; j < _horizontalSplitting; j++)
                {
                    for (int k = 0; k < _verticalSplitting; k++)
                    {
                        previousTemperature[j, k] = currentTemperature[j, k];
                    }
                }
            }

            currentTemperature = null;
            previousTemperature = null;
            g = null;

            if (calculatedCooldownCurve != null && calculatedCooldownCurveSizeRows > 0)
            {
                if (_referenceCooldownCurveSizeRows == calculatedCooldownCurveSizeRows)
                {
                    double[] sum = new double[_monitoredIndexRowNumber];

                    for (int i = 0; i < _monitoredIndexRowNumber; i++)
                    {
                        sum[i] = 0;
                    }

                    for (int i = 0; i < _referenceCooldownCurveSizeRows; i++)
                    {
                        for (int j = 1; j < (_monitoredIndexRowNumber + 1); j++)
                        {
                            sum[j - 1] += Math.Pow(Convert.ToDouble(_referenceCooldownCurve[(i * (_monitoredIndexRowNumber + 1)) + j].ToString("000.0")) -
                                Convert.ToDouble(calculatedCooldownCurve[(i * (_monitoredIndexRowNumber + 1)) + j].ToString("000.0")), 2);
                        }
                    }

                    calculatedCooldownCurve = null;

                    return sum;
                }
                else
                {
                    calculatedCooldownCurve = null;

                    return new double[] { double.MaxValue, double.MaxValue };
                }
            }
            else
            {
                calculatedCooldownCurve = null;

                return new double[] { double.MaxValue, double.MaxValue };
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
                double mBest1 = GetSwarmAverageBestPosition(1);
                double mBest2 = GetSwarmAverageBestPosition(2);

                for (int i = 0; i < _particleSwarmSize; i++)
                {
                    OptimizePosition(i, mBest1, mBest2, beta);
                }

                double improovedFitness1 = 0;
                double improovedFitness2 = 0;
                double[] improovedPosition1 = new double[_dimensionRowNumber];
                double[] improovedPosition2 = new double[_dimensionRowNumber];

                for (int i = 0; i < _particleSwarmSize; i++)
                {
                    if (_particleFitness1[i] < _globalBestFitness1[_globalBestSize - 1] &&
                        _particleFitness2[i] < _globalBestFitness2[_globalBestSize - 1])
                    {
                        improovedFitness1 = _particleBestFitness1[i];
                        improovedFitness2 = _particleBestFitness2[i];

                        for (int j = 0; j < _dimensionRowNumber; j++)
                        {
                            improovedPosition1[j] = _particleBestPosition1[(i * _dimensionRowNumber) + j];
                            improovedPosition2[j] = _particleBestPosition2[(i * _dimensionRowNumber) + j];
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
                    double[] tempFitness1 = new double[_globalBestSize];
                    double[] tempFitness2 = new double[_globalBestSize];
                    double[] tempPosition1 = new double[_globalBestSize * _dimensionRowNumber];
                    double[] tempPosition2 = new double[_globalBestSize * _dimensionRowNumber];

                    for (int i = 0; i < _globalBestSize; i++)
                    {
                        tempFitness1[i] = _globalBestFitness1[i];
                        tempFitness2[i] = _globalBestFitness2[i];

                        for (int j = 0; j < _dimensionRowNumber; j++)
                        {
                            tempPosition1[(i * _dimensionRowNumber) + j] = _globalBestPosition1[(i * _dimensionRowNumber) + j];
                            tempPosition2[(i * _dimensionRowNumber) + j] = _globalBestPosition2[(i * _dimensionRowNumber) + j];
                        }
                    }

                    _globalBestFitness1 = null;
                    _globalBestFitness2 = null;
                    _globalBestPosition1 = null;
                    _globalBestPosition2 = null;
                    _globalBestSize++;
                    _globalBestFitness1 = new double[_globalBestSize];
                    _globalBestFitness2 = new double[_globalBestSize];
                    _globalBestPosition1 = new double[_globalBestSize * _dimensionRowNumber];
                    _globalBestPosition2 = new double[_globalBestSize * _dimensionRowNumber];

                    for (int i = 0; i < _globalBestSize; i++)
                    {
                        if (i < _globalBestSize - 1)
                        {
                            _globalBestFitness1[i] = tempFitness1[i];
                            _globalBestFitness2[i] = tempFitness2[i];

                            for (int j = 0; j < _dimensionRowNumber; j++)
                            {
                                _globalBestPosition1[(i * _dimensionRowNumber) + j] = tempPosition1[(i * _dimensionRowNumber) + j];
                                _globalBestPosition2[(i * _dimensionRowNumber) + j] = tempPosition2[(i * _dimensionRowNumber) + j];
                            }
                        }
                        else
                        {
                            _globalBestFitness1[i] = improovedFitness1;
                            _globalBestFitness2[i] = improovedFitness2;
                            string position = "New best found at iteration: " + _epoch + ", with value1 of: " + improovedFitness1 + ", at ";
                            string position2 = " and with value2 of: " + improovedFitness2 + ", at ";

                            for (int j = 0; j < _dimensionRowNumber; j++)
                            {
                                position += "position" + (j + 1) + ": " + improovedPosition1[j] + " ";
                                position2 += "position" + (j + 1) + ": " + improovedPosition2[j] + " ";
                                _globalBestPosition1[(i * _dimensionRowNumber) + j] = improovedPosition1[j];
                                _globalBestPosition2[(i * _dimensionRowNumber) + j] = improovedPosition2[j];
                            }

                            Console.WriteLine(position + position2 + "\r\n");
                        }
                    }

                    tempFitness1 = null;
                    tempFitness2 = null;
                    tempPosition1 = null;
                    tempPosition2 = null;

                    if (Math.Abs(_globalBestFitness1[_globalBestSize - 1]) <= _particleEpsilon &&
                        Math.Abs(_globalBestFitness2[_globalBestSize - 1]) <= _particleEpsilon)
                    {
                        _exitReason = "The particle swarm optimization reached an acceptable fitness value. The value was given at: " + _particleEpsilon;

                        return;
                    }
                }

                improovedPosition1 = null;
                improovedPosition2 = null;
                _epoch++;

                if (_globalBestSize >= 5)
                {
                    if ((_globalBestFitness1[_globalBestSize - 1] / _globalBestFitness1[_globalBestSize - 5]) < 0.003 ||
                        (_globalBestFitness2[_globalBestSize - 1] / _globalBestFitness2[_globalBestSize - 5]) < 0.003)
                    {
                        _exitReason = "The PSO global fitness value changed less then 0.003% over the last 5 global value refresh.";

                        return;
                    }
                }

                //if (((double)_globalBestSize / (double)_epoch) < 0.03)
                //{
                //    _exitReason = "The convergence of the " + (_particleOptimalisationType == 1 ? "Clerc" : "Quantum") + " PSO is too slow, it may not find the optimum";

                //    return;
                //}

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
                    int p = i - n;

                    if (p < 0)
                    {
                        p = _particleSwarmSize + p;
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
                        p = p - _particleSwarmSize;
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
        /// <param name="particleIndex"></param>
        /// <param name="mbest"></param>
        /// <param name="beta"></param>
        private static void OptimizePosition(int particleIndex, double mBest1, double mBest2, double beta)
        {
            switch (_particleOptimalisationType)
            {
                case 1:
                    UpdateClercPosition(particleIndex);
                    break;
                case 2:
                    UpdateQuantumPosition(particleIndex, mBest1, mBest2, beta);
                    break;
            }

            UpdateErrorValues(particleIndex);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="particleIndex"></param>
        private static void UpdateClercPosition(int particleIndex)
        {
            UpdateVelocity(particleIndex, GetBestLocalPositions(particleIndex, 1), GetBestLocalPositions(particleIndex, 2));
            UpdatePosition(particleIndex);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="particleIndex"></param>
        private static void UpdateErrorValues(int particleIndex)
        {
            double[] checkPosition1 = new double[_dimensionRowNumber];
            double[] checkPosition2 = new double[_dimensionRowNumber];

            for (int i = 0; i < _dimensionRowNumber; i++)
            {
                checkPosition1[i] = _particlePosition1[(particleIndex * _dimensionRowNumber) + i];
                checkPosition2[i] = _particlePosition2[(particleIndex * _dimensionRowNumber) + i];
            }

            if (CheckIsInRange(checkPosition1) && CheckIsInRange(checkPosition2))
            {
                double[] fitnessValues = ObjectiveFunction(particleIndex);
                _particleFitness1[particleIndex] = fitnessValues[0];
                _particleFitness2[particleIndex] = fitnessValues[1];

                if (_particleFitness1[particleIndex] < _particleBestFitness1[particleIndex])
                {
                    _particleBestFitness1[particleIndex] = _particleFitness1[particleIndex];

                    for (int i = 0; i < _dimensionRowNumber; i++)
                    {
                        _particleBestPosition1[(particleIndex * _dimensionRowNumber) + i] = _particlePosition1[(particleIndex * _dimensionRowNumber) + i];
                    }
                }

                if (_particleFitness2[particleIndex] < _particleBestFitness2[particleIndex])
                {
                    _particleBestFitness2[particleIndex] = _particleFitness2[particleIndex];

                    for (int i = 0; i < _dimensionRowNumber; i++)
                    {
                        _particleBestPosition2[(particleIndex * _dimensionRowNumber) + i] = _particlePosition2[(particleIndex * _dimensionRowNumber) + i];
                    }
                }
            }

            checkPosition1 = null;
            checkPosition2 = null;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="particleIndex"></param>
        private static void UpdatePosition(int particleIndex)
        {
            for (int i = 0; i < _dimensionRowNumber; i++)
            {
                _particlePosition1[(particleIndex * _dimensionRowNumber) + i] =
                    (_particlePosition1[(particleIndex * _dimensionRowNumber) + i] +
                    _particleVelocity1[(particleIndex * _dimensionRowNumber) + i]);
                _particlePosition2[(particleIndex * _dimensionRowNumber) + i] =
                    (_particlePosition2[(particleIndex * _dimensionRowNumber) + i] +
                    _particleVelocity2[(particleIndex * _dimensionRowNumber) + i]);
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="particleIndex"></param>
        /// <param name="mBest1"></param>
        /// <param name="mBest2"></param>
        /// <param name="beta"></param>
        private static void UpdateQuantumPosition(int particleIndex, double mBest1, double mBest2, double beta)
        {
            for (int i = 0; i < _dimensionRowNumber; i++)
            {
                double fi = _random.NextDouble();
                double p = fi * _particleBestPosition1[(particleIndex * _dimensionRowNumber) + i] +
                    (1 - fi) * _globalBestPosition1[(_globalBestSize - 1 - _dimensionRowNumber) + i];

                if (fi > 0.5)
                {
                    _particlePosition1[(particleIndex * _dimensionRowNumber) + i] = p - beta *
                        Math.Abs(mBest1 - _particlePosition1[(particleIndex * _dimensionRowNumber) + i]) * (-Math.Log10(fi));
                }
                else
                {
                    _particlePosition1[(particleIndex * _dimensionRowNumber) + i] = p + beta *
                        Math.Abs(mBest1 - _particlePosition1[(particleIndex * _dimensionRowNumber) + i]) * (-Math.Log10(fi));
                }

                if (_particlePosition1[(particleIndex * _dimensionRowNumber) + i] < _rangeMin)
                {
                    _particlePosition1[(particleIndex * _dimensionRowNumber) + i] = 2 *
                        _rangeMin - _particlePosition1[(particleIndex * _dimensionRowNumber) + i];
                }

                if (_particlePosition1[(particleIndex * _dimensionRowNumber) + i] > _rangeMax)
                {
                    _particlePosition1[(particleIndex * _dimensionRowNumber) + i] = 2 *
                        _rangeMax - _particlePosition1[(particleIndex * _dimensionRowNumber) + i];
                }

                fi = _random.NextDouble();
                p = fi * _particleBestPosition2[(particleIndex * _dimensionRowNumber) + i] +
                    (1 - fi) * _globalBestPosition2[(_globalBestSize - 1 - _dimensionRowNumber) + i];

                if (fi > 0.5)
                {
                    _particlePosition2[(particleIndex * _dimensionRowNumber) + i] = p - beta *
                        Math.Abs(mBest2 - _particlePosition2[(particleIndex * _dimensionRowNumber) + i]) * (-Math.Log10(fi));
                }
                else
                {
                    _particlePosition2[(particleIndex * _dimensionRowNumber) + i] = p + beta *
                        Math.Abs(mBest2 - _particlePosition2[(particleIndex * _dimensionRowNumber) + i]) * (-Math.Log10(fi));
                }

                if (_particlePosition2[(particleIndex * _dimensionRowNumber) + i] < _rangeMin)
                {
                    _particlePosition2[(particleIndex * _dimensionRowNumber) + i] = 2 *
                        _rangeMin - _particlePosition2[(particleIndex * _dimensionRowNumber) + i];
                }

                if (_particlePosition2[(particleIndex * _dimensionRowNumber) + i] > _rangeMax)
                {
                    _particlePosition2[(particleIndex * _dimensionRowNumber) + i] = 2 *
                        _rangeMax - _particlePosition2[(particleIndex * _dimensionRowNumber) + i];
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="particleIndex"></param>
        /// <param name="bestLocalPosition1"></param>
        /// <param name="bestLocalPosition2"></param>
        private static void UpdateVelocity(int particleIndex, double[] bestLocalPosition1, double[] bestLocalPosition2)
        {
            for (int i = 0; i < _dimensionRowNumber; i++)
            {
                // particle.velocity[i] = (w * paritcle.velocity[i]) + (constant1 * rand() * (particle.bestPosition[i] - particle.position[i])) + (constant2 * rand() * bestPosition[i] - particle.position[i])
                _particleVelocity1[(particleIndex * _dimensionRowNumber) + i] =
                    (_weight * _particleVelocity1[(particleIndex * _dimensionRowNumber) + i]) +
                    (_particleConstant1 * _random.NextDouble() *
                    (_particleBestPosition1[(particleIndex * _dimensionRowNumber) + i] -
                    _particlePosition1[(particleIndex * _dimensionRowNumber) + i])) +
                    (_particleConstant2 * _random.NextDouble() *
                    (bestLocalPosition1[i] - _particlePosition1[(particleIndex * _dimensionRowNumber) + i]));
                _particleVelocity2[(particleIndex * _dimensionRowNumber) + i] =
                    (_weight * _particleVelocity2[(particleIndex * _dimensionRowNumber) + i]) +
                    (_particleConstant1 * _random.NextDouble() *
                    (_particleBestPosition2[(particleIndex * _dimensionRowNumber) + i] -
                    _particlePosition2[(particleIndex * _dimensionRowNumber) + i])) +
                    (_particleConstant2 * _random.NextDouble() *
                    (bestLocalPosition2[i] - _particlePosition2[(particleIndex * _dimensionRowNumber) + i]));
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="position"></param>
        /// <returns></returns>
        private static bool CheckIsInRange(double[] position)
        {
            for (int i = 0; i < _dimensionRowNumber; i++)
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
        /// <param name="index"></param>
        /// <returns></returns>
        private static double[] GetBestLocalPositions(int particleIndex, int index)
        {
            double[] returnValue = new double[_dimensionRowNumber];
            int bestIndex = particleIndex;

            switch (index)
            {
                case 1:
                    for (int i = 0; i < _particleInformerNumber; i++)
                    {
                        if (_particleBestFitness1[_particleInformers[(particleIndex * _particleInformerNumber) + i]] < _particleBestFitness1[bestIndex])
                        {
                            bestIndex = _particleInformers[(particleIndex * _particleInformerNumber) + i];
                        }
                    }

                    for (int i = 0; i < _dimensionRowNumber; i++)
                    {
                        returnValue[i] = _particleBestPosition1[(bestIndex * _dimensionRowNumber) + i];
                    }
                    break;
                case 2:
                    for (int i = 0; i < _particleInformerNumber; i++)
                    {
                        if (_particleBestFitness2[_particleInformers[(particleIndex * _particleInformerNumber) + i]] < _particleBestFitness2[bestIndex])
                        {
                            bestIndex = _particleInformers[(particleIndex * _particleInformerNumber) + i];
                        }
                    }

                    for (int i = 0; i < _dimensionRowNumber; i++)
                    {
                        returnValue[i] = _particleBestPosition2[(bestIndex * _dimensionRowNumber) + i];
                    }
                    break;
            }

            return returnValue;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="particleIndex"></param>
        public static void GetInitialPosition(int particleIndex)
        {
            for (int i = 0; i < _dimensionRowNumber; i++)
            {
                // The given values are null.
                if (_particleInitialPosition == null)
                {
                    _particlePosition1[(particleIndex * _dimensionRowNumber) + i] = ((_rangeMax - _rangeMin) * _random.NextDouble() + _rangeMin);
                }
                // There is only one given value.
                else if (_particleInitialPositionSize == 1)
                {
                    _particlePosition1[(particleIndex * _dimensionRowNumber) + i] = _particleInitialPosition[0];
                }
                // There are as many given value as the position dimensions.
                else if (_particleInitialPositionSize == _dimensionRowNumber)
                {
                    _particlePosition1[(particleIndex * _dimensionRowNumber) + i] = _particleInitialPosition[i];
                }
                // The current position can be set from the given values.
                else if (i < _particleInitialPositionSize)
                {
                    _particlePosition1[(particleIndex * _dimensionRowNumber) + i] = _particleInitialPosition[i];
                }
                // The current position can't be set from the given values.
                else
                {
                    _particlePosition1[(particleIndex * _dimensionRowNumber) + i] = ((_rangeMax - _rangeMin) * _random.NextDouble() + _rangeMin);
                }
            }

            for (int i = 0; i < _dimensionRowNumber; i++)
            {
                // The given values are null.
                if (_particleInitialPosition == null)
                {
                    _particlePosition2[(particleIndex * _dimensionRowNumber) + i] = ((_rangeMax - _rangeMin) * _random.NextDouble() + _rangeMin);
                }
                // There is only one given value.
                else if (_particleInitialPositionSize == 1)
                {
                    _particlePosition2[(particleIndex * _dimensionRowNumber) + i] = _particleInitialPosition[0];
                }
                // There are as many given value as the position dimensions.
                else if (_particleInitialPositionSize == _dimensionRowNumber)
                {
                    _particlePosition2[(particleIndex * _dimensionRowNumber) + i] = _particleInitialPosition[i];
                }
                // The current position can be set from the given values.
                else if (i < _particleInitialPositionSize)
                {
                    _particlePosition2[(particleIndex * _dimensionRowNumber) + i] = _particleInitialPosition[i];
                }
                // The current position can't be set from the given values.
                else
                {
                    _particlePosition2[(particleIndex * _dimensionRowNumber) + i] = ((_rangeMax - _rangeMin) * _random.NextDouble() + _rangeMin);
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="particleIndex"></param>
        public static void GetInitialVelocity(int particleIndex)
        {
            for (int i = 0; i < _dimensionRowNumber; i++)
            {
                // The given values are null.
                if (_particleInitialVelocity == null)
                {
                    _particleVelocity1[(particleIndex * _dimensionRowNumber) + i] = ((_rangeMax - _rangeMin) * _random.NextDouble() + _rangeMin);
                }
                //  There is only one given value.
                else if (_particleInitialVelocitySize == 1)
                {
                    _particleVelocity1[(particleIndex * _dimensionRowNumber) + i] = _particleInitialVelocity[0];
                }
                // There are as many given value as the velocity dimensions.
                else if (_particleInitialVelocitySize == _dimensionRowNumber)
                {
                    _particleVelocity1[(particleIndex * _dimensionRowNumber) + i] = _particleInitialVelocity[i];
                }
                // The current velocity can be set from the given values.
                else if (i < _particleInitialVelocitySize)
                {
                    _particleVelocity1[(particleIndex * _dimensionRowNumber) + i] = _particleInitialVelocity[i];
                }
                // The current velocity can't be set from the given values.
                else
                {
                    _particleVelocity1[(particleIndex * _dimensionRowNumber) + i] = ((_rangeMax - _rangeMin) * _random.NextDouble() + _rangeMin);
                }
            }

            for (int i = 0; i < _dimensionRowNumber; i++)
            {
                // The given values are null.
                if (_particleInitialVelocity == null)
                {
                    _particleVelocity2[(particleIndex * _dimensionRowNumber) + i] = ((_rangeMax - _rangeMin) * _random.NextDouble() + _rangeMin);
                }
                //  There is only one given value.
                else if (_particleInitialVelocitySize == 1)
                {
                    _particleVelocity2[(particleIndex * _dimensionRowNumber) + i] = _particleInitialVelocity[0];
                }
                // There are as many given value as the velocity dimensions.
                else if (_particleInitialVelocitySize == _dimensionRowNumber)
                {
                    _particleVelocity2[(particleIndex * _dimensionRowNumber) + i] = _particleInitialVelocity[i];
                }
                // The current velocity can be set from the given values.
                else if (i < _particleInitialVelocitySize)
                {
                    _particleVelocity2[(particleIndex * _dimensionRowNumber) + i] = _particleInitialVelocity[i];
                }
                // The current velocity can't be set from the given values.
                else
                {
                    _particleVelocity2[(particleIndex * _dimensionRowNumber) + i] = ((_rangeMax - _rangeMin) * _random.NextDouble() + _rangeMin);
                }
            }
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
        /// <param name="index"></param>
        /// <returns></returns>
        private static double GetSwarmAverageBestPosition(int index)
        {
            double sum = 0;

            for (int i = 0; i < _particleSwarmSize; i++)
            {
                for (int j = 0; j < _dimensionRowNumber; j++)
                {
                    switch (index)
                    {
                        case 1:
                            sum += _particleBestPosition1[(i * _dimensionRowNumber) + j];
                            break;
                        case 2:
                            sum += _particleBestPosition2[(i * _dimensionRowNumber) + j];
                            break;
                    }
                }
            }

            return ((sum * 1.0) / _particleSwarmSize);
        }

        #endregion Methods
    }
}
