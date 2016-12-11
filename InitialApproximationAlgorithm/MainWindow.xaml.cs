using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using Window = System.Windows.Window;

namespace InitialApproximationAlgorithm
{
    public struct Circle
    {
        public double Center { get; set; }
        public double Radius { get; set; }
    }

    public class InitialСonditions
    {
        public Circle Circle { get; set; }
        public Func<Complex, Complex> OperatorCoeficientsFunction { get; set; }
        public Func<Complex, Complex> DerivedOperatorCoeficientsFunction { get; set; }
        public Func<Complex, Complex> SecondDerivedOperatorCoeficientsFunction { get; set; }
        public int NumberPartitions { get; set; }
        public double Eps { get; set; }
        public double L0 { get; set; }
    }

    public class Decomposition
    {
        public Matrix<Complex> L { get; set; }
        public Matrix<Complex> U { get; set; }
        public Matrix<Complex> V { get; set; }
        public Matrix<Complex> M { get; set; }
        public Matrix<Complex> W { get; set; }
        public Matrix<Complex> N { get; set; }

        public Decomposition(int demension)
        {
            L = Matrix<Complex>.Build.DenseIdentity(demension, demension);
            M = Matrix<Complex>.Build.DenseIdentity(demension, demension);
            N = Matrix<Complex>.Build.DenseIdentity(demension, demension);
            U = Matrix<Complex>.Build.Dense(demension, demension);
            V = Matrix<Complex>.Build.Dense(demension, demension);
            W = Matrix<Complex>.Build.Dense(demension, demension);
        }
    }

    public struct Result
    {
        public double InitialEigenValue { get; set; }
        public double ApproximatedEigenValue { get; set; }
        public int CountIteration { get; set; }
        public double LowerBound { get; set; }
        public double UpperBound { get; set; }
    }

    public partial class MainWindow : Window
    {
        Object lockMe = new Object();
        public MainWindow()
        {
            InitializeComponent();
        }

        List<Result> results = new List<Result>();

        private void CalculateButtonClick(object sender, RoutedEventArgs e)
        {
            InitialСonditions initialСonditions = new InitialСonditions
            {
                Circle = new Circle()
                {
                    Center = double.Parse(CenterTextBox.Text.Replace('.', ',')),
                    Radius = double.Parse(RadiusTextBox.Text.Replace('.', ','))
                },
                OperatorCoeficientsFunction = x => -(4 + x * x) * (4 + x * x),
                DerivedOperatorCoeficientsFunction = x => -4 * (4 + x * x) * x,
                SecondDerivedOperatorCoeficientsFunction = x => -16 - 12 * x * x,
                NumberPartitions = int.Parse(NumberPartitionsTextBox.Text),
                Eps = double.Parse(EpsTextBox.Text.Replace('.', ','))
            };

            var results = GetInitialApproximationWithAproximatedEigenValues(initialСonditions);
            AddRowsToEigenValuesTable(results);
        }

        private List<Result> GetInitialApproximationWithAproximatedEigenValues(InitialСonditions initialСonditions)
        {
            
            var sk0 = GetSkValue(initialСonditions, 0);
            int countEigenValues = (int)Math.Round(sk0);
            var initialEigenValues = new List<Complex>();
            var sk = new List<double>();
            if (countEigenValues > 0)
            {
                for (int j = 1; j <= countEigenValues; j++)
                {
                    var skj = GetSkValue(initialСonditions, j);
                    sk.Add(skj);
                    //double jPartition = j * 1.0 / countEigenValues;

                    //var jInitialLambda = initialСonditions.Circle.Center +
                    //GetSpectralRadius(initialСonditions.Circle.Radius, jPartition);

                    //initialEigenValues.Add(jInitialLambda);
                    initialСonditions.L0 = skj;

                    var result = SolveNewtonMethod(initialСonditions);
                    result.LowerBound = initialСonditions.Circle.Center - initialСonditions.Circle.Radius;
                    result.UpperBound = initialСonditions.Circle.Center + initialСonditions.Circle.Radius;

                    results.Add(result);
                }
            }
            else
            {
                InitialСonditions newInitialСonditions = new InitialСonditions
                {
                    Circle = new Circle()
                    {
                        Center = initialСonditions.Circle.Center + initialСonditions.Circle.Radius,
                        Radius = initialСonditions.Circle.Radius
                    },
                    OperatorCoeficientsFunction = x => -(4 + x * x) * (4 + x * x),
                    DerivedOperatorCoeficientsFunction = x => -4 * (4 + x * x) * x,
                    SecondDerivedOperatorCoeficientsFunction = x => -16 - 12 * x * x,
                    NumberPartitions = initialСonditions.NumberPartitions,
                    Eps = initialСonditions.Eps
                };

                results.Add(new Result() { ApproximatedEigenValue = 0, CountIteration = 0, InitialEigenValue = 0, LowerBound = initialСonditions.Circle.Center - initialСonditions.Circle.Radius, UpperBound = initialСonditions.Circle.Center + initialСonditions.Circle.Radius });

                GetInitialApproximationWithAproximatedEigenValues(newInitialСonditions);
            }

            #region Parallel
            /*Parallel.For(0, countEigenValues, (j) =>
            {
                double jPartition = j * 1.0 / countEigenValues;

                var jInitialLambda = initialСonditions.Circle.Center +
                                     GetSpectralRadius(initialСonditions.Circle.Radius, jPartition);

                initialСonditions.L0 = jInitialLambda;

                var result = SolveNewtonMethod(initialСonditions);

                lock (lockMe)
                {
                    results.Add(result);
                }
            });*/
            #endregion

            return results;
        }

        private double GetSkValue(InitialСonditions initialСonditions, int k)
        {
            Complex sum = Complex.Zero;

            var numberPartitions = initialСonditions.NumberPartitions;
            var radius = initialСonditions.Circle.Radius;
            var center = initialСonditions.Circle.Center;

            var xPartitions = GetPartitionForIntegrationLimit(initialСonditions);

            for (int j = 1; j <= numberPartitions; j++)
            {
                var spectralRadius = GetSpectralRadius(radius, j, numberPartitions);
                var jLambda = center + spectralRadius;

                var matrixD = GetMatrixForDecomposition(xPartitions, initialСonditions.OperatorCoeficientsFunction,
                    jLambda);

                var matrixB = GetMatrixForDecomposition(xPartitions, initialСonditions.DerivedOperatorCoeficientsFunction,
                    jLambda);

                var matrixC = GetMatrixForDecomposition(xPartitions, initialСonditions.SecondDerivedOperatorCoeficientsFunction,
                    jLambda);

                var decomposition = GetDecompositionMatrices(matrixD, matrixB, matrixC);

                var sumDividedDiagonalElements = GetSumForDivideVUDiagonal(decomposition);

                sum += Complex.Pow(jLambda, k) * spectralRadius * sumDividedDiagonalElements;
            }

            #region Comment
            /*var xPartitions = GetPartitionForIntegrationLimit(initialСonditions);

            for (int j = 0; j <= initialСonditions.NumberPartitions; j++)
            {
                var jPartition = j * 1.0 / initialСonditions.NumberPartitions;

                var spectralRadius = GetSpectralRadius(initialСonditions.Circle.Radius, jPartition);

                var jLambda = initialСonditions.Circle.Center + spectralRadius;

                var matrixD = GetMatrixForDecomposition(xPartitions, initialСonditions.OperatorCoeficientsFunction,
                    jLambda);

                var matrixB = GetMatrixForDecomposition(xPartitions, initialСonditions.DerivedOperatorCoeficientsFunction,
                    jLambda);

                var matrixC = GetMatrixForDecomposition(xPartitions, initialСonditions.SecondDerivedOperatorCoeficientsFunction,
                    jLambda);

                var decomposition = GetDecompositionMatrices(matrixD, matrixB, matrixC);

                var sumDividedDiagonalElements = GetSumForDivideVUDiagonal(decomposition);

                sum += spectralRadius * sumDividedDiagonalElements;
            }*/
            #endregion

            #region Parallel
            /*Parallel.For(1, initialСonditions.NumberPartitions, (j) =>
            {
                var jPartition = j * 1.0 / initialСonditions.NumberPartitions;
                var spectralRadius = GetSpectralRadius(initialСonditions.Circle.Radius, jPartition);
                var jLambda = initialСonditions.Circle.Center + spectralRadius;

                var matrixD = GetMatrixForDecomposition(xPartitions, initialСonditions.OperatorCoeficientsFunction,
                    jLambda.Magnitude);

                var matrixB = GetMatrixForDecomposition(xPartitions,
                    initialСonditions.DerivedOperatorCoeficientsFunction,
                    jLambda.Magnitude);

                var matrixC = GetMatrixForDecomposition(xPartitions,
                    initialСonditions.SecondDerivedOperatorCoeficientsFunction,
                    jLambda.Magnitude);

                var decomposition = GetDecompositionMatrices(matrixD, matrixB, matrixC);

                var deltaLambda = 1 / GetDeltaLambda(decomposition);

                lock (lockMe)
                {
                    sum += spectralRadius * deltaLambda;
                }
            });*/
            #endregion

            sum /= initialСonditions.NumberPartitions * initialСonditions.NumberPartitions;

            return sum.Magnitude;

        }

        private Vector<double> GetPartitionForIntegrationLimit(InitialСonditions initialСonditions)
        {
            var step = 1.0 / initialСonditions.NumberPartitions;
            return Vector<double>.Build.Dense(
                initialСonditions.NumberPartitions + 1,
                i => i * step);
            //var step = 2 * initialСonditions.Circle.Radius / initialСonditions.NumberPartitions;
            //var startPoint = initialСonditions.Circle.Center - initialСonditions.Circle.Radius;

            //return Vector<double>.Build.Dense(
            //    initialСonditions.NumberPartitions + 1,
            //    i => startPoint + i * step);
        }

        private Complex GetSpectralRadius(double radius, double j, int numberPartition)
        {
            return radius * Complex.Exp(Complex.ImaginaryOne * 2 * Math.PI * j / numberPartition);
        }

        private Matrix<Complex> GetMatrixForDecomposition(Vector<double> xPartitions, Func<Complex, Complex> function, Complex lambda)
        {
            int dimensional = xPartitions.Count;

            var matrixA = GetMatrixByFunction(xPartitions, function);
            var matrixLambdaI = Matrix<Complex>.Build.DenseDiagonal(dimensional, lambda);

            return matrixA - matrixLambdaI;
        }

        private Matrix<Complex> GetMatrixByFunction(Vector<double> xPartitions, Func<Complex, Complex> function)
        {
            var matrixDimensional = xPartitions.Count;
            Matrix<Complex> operatorMatrix = Matrix<Complex>.Build.Dense(matrixDimensional, matrixDimensional);
            var step = (xPartitions[matrixDimensional - 1] - xPartitions[0]) / (matrixDimensional - 1);
            var squaredStep = step * step;
            operatorMatrix[0, 0] = function(xPartitions[0]);
            operatorMatrix[matrixDimensional - 1, matrixDimensional - 1] = function(xPartitions[matrixDimensional - 1]);
            for (int i = 1; i < matrixDimensional - 1; i++)
            {

                operatorMatrix[i, i - 1] = 1.0 / squaredStep * function(xPartitions[i]);
                operatorMatrix[i, i] = -2.0 / squaredStep * function(xPartitions[i]);
                operatorMatrix[i, i + 1] = 1.0 / squaredStep * function(xPartitions[i]);
            }

            return operatorMatrix;
        }

        private Decomposition GetDecompositionMatrices(Matrix<Complex> matrixD, Matrix<Complex> matrixB, Matrix<Complex> matrixC)
        {
            int dimensional = matrixD.ColumnCount;

            Decomposition decomposition = new Decomposition(dimensional);

            for (int r = 0; r < dimensional; r++)
            {
                for (int k = r; k < dimensional; k++)
                {
                    var sumForU = Complex.Zero;
                    var sumForV = Complex.Zero;
                    var sumForW = Complex.Zero;

                    for (int j = 0; j < r; j++)
                    {
                        sumForU += decomposition.L[r, j] * decomposition.U[j, k];

                        sumForV += decomposition.M[r, j] * decomposition.U[j, k] +
                                   decomposition.L[r, j] * decomposition.V[j, k];

                        sumForW += decomposition.N[r, j] * decomposition.U[j, k] +
                                   2 * decomposition.M[r, j] * decomposition.V[j, k] +
                                   decomposition.L[r, j] * decomposition.W[j, k];
                    }

                    decomposition.U[r, k] = matrixD[r, k] - sumForU;
                    decomposition.V[r, k] = matrixB[r, k] - sumForV;
                    decomposition.W[r, k] = matrixC[r, k] - sumForW;
                }

                for (int i = r + 1; i < dimensional; i++)
                {
                    var sumForL = Complex.Zero;
                    var sumForM = Complex.Zero;
                    var sumForN = Complex.Zero;

                    for (int j = 0; j < r; j++)
                    {
                        sumForL += decomposition.L[i, j] * decomposition.U[j, r];

                        sumForM += decomposition.M[i, j] * decomposition.U[j, r] +
                                   decomposition.L[i, j] * decomposition.V[j, r];

                        sumForN += decomposition.N[i, j] * decomposition.U[j, r] +
                                   2 * decomposition.M[i, j] * decomposition.V[j, r] +
                                   decomposition.L[i, j] * decomposition.W[j, r];
                    }

                    decomposition.L[i, r] = (matrixD[i, r] - sumForL) / decomposition.U[r, r];

                    decomposition.M[i, r] = (matrixB[i, r] -
                                                sumForM -
                                                decomposition.L[i, r] * decomposition.V[r, r]) /
                                            decomposition.U[r, r];

                    decomposition.N[i, r] = (matrixC[i, r] -
                                                sumForN -
                                                2 * decomposition.M[i, r] * decomposition.V[r, r]
                                                - decomposition.L[i, r] * decomposition.W[r, r]) /
                                            decomposition.U[r, r];
                }
            }

            return decomposition;
        }

        private Complex GetSumForDivideVUDiagonal(Decomposition decomposition)
        {
            var sum = Complex.Zero;

            int dimension = decomposition.V.ColumnCount;

            for (int i = 0; i < dimension; i++)
            {
                sum += decomposition.V[i, i] / decomposition.U[i, i];
            }

            return sum;
        }

        private Result SolveNewtonMethod(InitialСonditions initialСonditions)
        {
            var eigenValues = new List<Complex>();
            eigenValues.Add(initialСonditions.L0);

            var xPartitions = GetPartitionForIntegrationLimit(initialСonditions);

            Complex prevEigenValue;
            Complex nextEigenValue;

            do
            {
                prevEigenValue = eigenValues.Last();

                var matrixD = GetMatrixForDecomposition(xPartitions, initialСonditions.OperatorCoeficientsFunction,
                    prevEigenValue);

                var matrixB = GetMatrixForDecomposition(xPartitions, initialСonditions.DerivedOperatorCoeficientsFunction,
                    prevEigenValue);

                var matrixC = GetMatrixForDecomposition(xPartitions, initialСonditions.SecondDerivedOperatorCoeficientsFunction,
                    prevEigenValue);

                var decomposition = GetDecompositionMatrices(matrixD, matrixB, matrixC);

                var deltaLambda = GetDeltaLambda(decomposition);

                nextEigenValue = prevEigenValue - deltaLambda;

                eigenValues.Add(nextEigenValue);
            } while ((nextEigenValue - prevEigenValue).Magnitude > initialСonditions.Eps/* && eigenValues.Count < 500*/);

            return new Result()
            {
                ApproximatedEigenValue = eigenValues.Last().Real,
                CountIteration = eigenValues.Count,
                InitialEigenValue = eigenValues.First().Real
            };
        }

        private Complex GetDeltaLambda(Decomposition decomposition)
        {
            Complex sum = GetSumForDivideVUDiagonal(decomposition);

            var deltaLambda = 1.0 / sum;

            return deltaLambda;
        }

        private void AddRowsToEigenValuesTable(List<Result> results)
        {
            var resultForTable = results.Select((x, i) => new { N = i, Interval = string.Format("[{0}, {1}]", x.LowerBound, x.UpperBound), x.ApproximatedEigenValue, x.InitialEigenValue, x.CountIteration });
            EigenValuesTable.ItemsSource = resultForTable;
        }
    }
}
