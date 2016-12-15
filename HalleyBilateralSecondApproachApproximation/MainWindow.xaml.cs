using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows;
using MathNet.Numerics.LinearAlgebra;

namespace HalleyBilateralSecondApproachApproximation
{
    public struct IntegrationLimits
    {
        public double UpperLimit { get; set; }
        public double LowerLimit { get; set; }
    }
    public class InitialСonditions
    {
        public IntegrationLimits IntegrationLimits { get; set; }
        public Func<double, double> OperatorCoeficientsFunction { get; set; }
        public Func<double, double> DerivedOperatorCoeficientsFunction { get; set; }
        public Func<double, double> SecondDerivedOperatorCoeficientsFunction { get; set; }
        public int NumberPartitions { get; set; }
        public double Eps { get; set; }
        public double L0 { get; set; }
    }
    public class Decomposition
    {
        public Matrix<double> L { get; set; }
        public Matrix<double> U { get; set; }
        public Matrix<double> V { get; set; }
        public Matrix<double> M { get; set; }
        public Matrix<double> W { get; set; }
        public Matrix<double> N { get; set; }

        public Decomposition(int demension)
        {
            L = Matrix<double>.Build.DenseIdentity(demension, demension);
            M = Matrix<double>.Build.DenseIdentity(demension, demension);
            N = Matrix<double>.Build.DenseIdentity(demension, demension);
            U = Matrix<double>.Build.Dense(demension, demension);
            V = Matrix<double>.Build.Dense(demension, demension);
            W = Matrix<double>.Build.Dense(demension, demension);
        }
    }
    public struct Result
    {
        public double Mu { get; set; }
        public double Lambda { get; set; }
        public double Nu { get; set; }
    }
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
        }

        private void CalculateButtonClick(object sender, RoutedEventArgs e)
        {
            InitialСonditions initialСonditions = new InitialСonditions
            {
                IntegrationLimits = new IntegrationLimits()
                {
                    LowerLimit = double.Parse(LowerLimitTextBox.Text),
                    UpperLimit = double.Parse(UpperLimitTextBox.Text)
                },
                OperatorCoeficientsFunction = x => -Math.Pow(4 + x * x, 2),
                DerivedOperatorCoeficientsFunction = x => 1,
                SecondDerivedOperatorCoeficientsFunction = x => 0,
                NumberPartitions = int.Parse(NumberPartitionsTextBox.Text),
                Eps = double.Parse(EpsTextBox.Text.Replace('.', ',')),
                L0 = double.Parse(Lambda0TextBox.Text.Replace('.', ','))
            };

            var eigenValues = SolveBilateralSecondApproachMethod(initialСonditions);
            AddRowToEigenValuesTable(eigenValues);
        }

        private Matrix<double> GetMatrixByFunction(Vector<double> xPartitions, Func<double, double> function)
        {
            var matrixDimensional = xPartitions.Count;
            Matrix<double> operatorMatrix = Matrix<double>.Build.Dense(matrixDimensional, matrixDimensional);
            var step = (xPartitions[matrixDimensional - 1] - xPartitions[0]) / (matrixDimensional - 1);
            var squaredStep = step * step;
            operatorMatrix[0, 0] = function(xPartitions[0]);
            operatorMatrix[matrixDimensional - 1, matrixDimensional - 1] = function(xPartitions[matrixDimensional - 1]);
            for (int i = 1; i < matrixDimensional - 1; i++)
            {
                operatorMatrix[i, i - 1] = -1.0 /*/ squaredStep*/ * function(xPartitions[i]);
                operatorMatrix[i, i] = 2.0 /*/ squaredStep*/ * function(xPartitions[i]);
                operatorMatrix[i, i + 1] = -1.0 /*/ squaredStep*/ * function(xPartitions[i]);
            }

            return operatorMatrix;
        }

        private Matrix<double> GetMatrixForDecomposition(Vector<double> xPartitions, Func<double, double> function, double lambda)
        {
            int dimensional = xPartitions.Count;

            var matrixA = GetMatrixByFunction(xPartitions, function);

            var step = (xPartitions[dimensional - 1] - xPartitions[0]) / (dimensional - 1);
            var squaredStep = step * step;

            var matrixLambdaI = Matrix<double>.Build.DenseDiagonal(dimensional, lambda * squaredStep);

            return matrixA + matrixLambdaI;
        }

        private double GetDeltaLambdaMu(double determinant, double determinantDerivative, double determinantSecondDerivative)
        {
            var numerator = determinant * determinantDerivative;
            var denominator = Math.Pow(determinantDerivative, 2.0) - determinant * determinantSecondDerivative;

            return numerator / denominator;
        }        

        private Vector<double> GetPartitionForIntegrationLimit(InitialСonditions initialСonditions)
        {
            var step = (initialСonditions.IntegrationLimits.UpperLimit - initialСonditions.IntegrationLimits.LowerLimit) / initialСonditions.NumberPartitions;
            return Vector<double>.Build.Dense(
                initialСonditions.NumberPartitions + 1,
                i => initialСonditions.IntegrationLimits.LowerLimit + i * step);
        }

        private Decomposition GetDecompositionMatrices(Matrix<double> matrixD, Matrix<double> matrixB, Matrix<double> matrixC)
        {
            int dimensional = matrixD.ColumnCount;

            Decomposition decomposition = new Decomposition(dimensional);

            for (int r = 0; r < dimensional; r++)
            {
                for (int k = r; k < dimensional; k++)
                {
                    double sumForU = 0;
                    double sumForV = 0;
                    double sumForW = 0;

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
                    double sumForL = 0;
                    double sumForM = 0;
                    double sumForN = 0;

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

            /*Console.WriteLine('D');
            Console.WriteLine(matrixD);
            Console.WriteLine('L');
            Console.WriteLine(decomposition.L);
            Console.WriteLine('U');
            Console.WriteLine(decomposition.U);*/

            //bool a = isEqual(matrixD, decomposition.L, decomposition.U);

            return decomposition;
        }

        private Matrix<double> SpectralOperatorSecondDerivative(int dimensional)
        {
            return Matrix<double>.Build.DenseDiagonal(dimensional, 0.0);
        }

        private Matrix<double> DerivedSpectralOperator(Vector<double> xPartitions)
        {
            int dimensional = xPartitions.Count;

            var step = (xPartitions[dimensional - 1] - xPartitions[0]) / (dimensional - 1);
            var squaredStep = step * step;

            var matrixLambdaI = Matrix<double>.Build.DenseDiagonal(dimensional, squaredStep);

            return matrixLambdaI;
        }

        private List<Result> SolveBilateralSecondApproachMethod(InitialСonditions initialСonditions)
        {
            var eigenValues = new List<Result>();

            var xPartitions = GetPartitionForIntegrationLimit(initialСonditions);

            var lambda = initialСonditions.L0;

            var matrixD = GetMatrixForDecomposition(xPartitions, initialСonditions.OperatorCoeficientsFunction,
                    lambda);

            var matrixB = DerivedSpectralOperator(xPartitions);

            var matrixC = SpectralOperatorSecondDerivative(xPartitions.Count);

            var decomposition = GetDecompositionMatrices(matrixD, matrixB, matrixC);

            var determinant = Determinant(decomposition.U);
            var determinantDerivative = DeterminantDerivative(decomposition.V, decomposition.U);
            var determinantSecondDerivative = DeterminantSecondDerivative(decomposition.W, decomposition.V, decomposition.U);

            var approximation = lambda;
            lambda = approximation - HalleyRelation(decomposition);

            double mu;
            if (determinant * determinantSecondDerivative > 0.0)
            {
                mu = approximation - determinant / determinantDerivative;
            }
            else
            {
                mu = approximation - GetDeltaLambdaMu(determinant, determinantDerivative, determinantSecondDerivative);
            }

            var nu = lambda;
            eigenValues.Add(new Result() {Mu = mu, Lambda = lambda, Nu = nu });

            while (Math.Abs(lambda - approximation) > initialСonditions.Eps)
            {
                matrixD = GetMatrixForDecomposition(xPartitions, initialСonditions.OperatorCoeficientsFunction,
                    mu);
                decomposition = GetDecompositionMatrices(matrixD, matrixB, matrixC);
                mu = mu - HalleyRelation(decomposition);

                matrixD = GetMatrixForDecomposition(xPartitions, initialСonditions.OperatorCoeficientsFunction,
                    nu);
                decomposition = GetDecompositionMatrices(matrixD, matrixB, matrixC);
                nu = nu - HalleyRelation(decomposition);

                approximation = lambda;
                lambda = (mu + nu) / 2.0;

                eigenValues.Add(new Result() { Mu = mu, Lambda = lambda, Nu = nu });
            };

            return eigenValues;
        }

        public static double Determinant(Matrix<double> u)
        {
            var result = 1.0;

            for (var i = 0; i < u.RowCount; ++i)
            {
                result *= u[i, i];
            }

            return result;
        }

        public static double DeterminantDerivative(Matrix<double> v, Matrix<double> u)
        {
            var summation = 0.0;

            for (var k = 0; k < v.RowCount; ++k)
            {
                var product = 1.0;

                for (var i = 0; i < u.RowCount; ++i)
                {
                    if (i != k)
                    {
                        product *= u[i, i];
                    }
                }

                summation += v[k, k] * product;
            }

            return summation;
        }

        public static double DeterminantSecondDerivative(Matrix<double> w, Matrix<double> v, Matrix<double> u)
        {
            var leftSummation = 0.0;

            for (var k = 0; k < w.RowCount; ++k)
            {
                var product = 1.0;

                for (var i = 0; i < u.RowCount; ++i)
                {
                    if (i != k)
                    {
                        product *= u[i, i];
                    }
                }

                leftSummation += w[k, k] * product;
            }

            var rightSummation = 0.0;

            for (var k = 0; k < v.RowCount; ++k)
            {
                var externalSummation = 0.0;

                for (var j = 0; j < v.RowCount; ++j)
                {
                    if (j != k)
                    {
                        var internalProduct = 1.0;

                        for (var i = 0; i < u.RowCount; ++i)
                        {
                            if (i != k && i != j)
                            {
                                internalProduct *= u[i, i];
                            }
                        }

                        externalSummation += v[j, j] * internalProduct;
                    }
                }

                rightSummation += v[k, k] * externalSummation;
            }

            return leftSummation + rightSummation;
        }

        private void AddRowToEigenValuesTable(List<Result> eigenValues)
        {
            /*var eigenValue = eigenValues.Mu.Last();

            List<string> leftApproximation = eigenValues.Mu.Select(x => x.ToString()).ToList();
            List<string> rightApproximation = eigenValues.Nu.Select(x => x.ToString()).ToList();

            leftApproximation.Add("");
            rightApproximation.Add("");

            List<string> eigenValueList = Enumerable.Repeat("", leftApproximation.Count).ToList();
            eigenValueList[eigenValueList.Count() - 1] = eigenValue.ToString();*/

            var resultsForTable = eigenValues.Select((x, i) => new { N = i + 1, Mu = x.Mu, Nu = x.Nu, ApproximatedEigenValue = x.Lambda });

            EigenValuesTable.ItemsSource = resultsForTable;
        }

        private double HalleyRelation(Decomposition decomposition)
        {
            double sum = 0;
            double sum2 = 0;

            int dimension = decomposition.V.ColumnCount;

            for (int i = 0; i < dimension; i++)
            {
                sum += decomposition.V[i, i] / decomposition.U[i, i];

                var sum3 = 0.0;
                for (int j = 0; j < dimension; j++)
                {
                    if (j != i)
                    {
                        sum3 += decomposition.V[j, j] / decomposition.U[j, j];
                    }
                }

                sum2 += Math.Pow(decomposition.V[i, i] / decomposition.U[i, i], 2) - 0.5 * decomposition.W[i, i] / decomposition.U[i, i] + 0.5 * decomposition.V[i, i] / decomposition.U[i, i] * sum3;
            }

            double deltaLambda = sum / sum2;

            return deltaLambda;
        }
    }
}
