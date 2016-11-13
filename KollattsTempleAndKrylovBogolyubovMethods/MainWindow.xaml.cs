using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Windows;

namespace KollattsTempleAndKrylovBogolyubovMethods
{
    using System.Linq;

    public struct IntegrationLimits
    {
        public double UpperLimit { get; set; }
        public double LowerLimit { get; set; }
    }
    public struct BoundaryConditions
    {
        public double LeftCondition { get; set; }
        public double RightCondition { get; set; }
    }
    public class InitialСonditions
    {
        public IntegrationLimits IntegrationLimits { get; set; }
        public Func<double, double> InitialApproximation { get; set; }
        public Func<double, double> OperatorCoeficientsFunction { get; set; }
        public BoundaryConditions BoundaryConditions { get; set; }
        public int NumberPartitions { get; set; }
        public double Eps { get; set; }
        public double L2 { get; set; }
    }
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow
    {
        private struct Bound
        {
            public double UpperBound { get; set; }
            public double LowerBound { get; set; }
        }
        private class Result
        {
            public List<Vector<double>> EigenVectors { get; }
            public List<double> SchwartzConstants { get; }
            public List<double> SchwartzRelations { get; }
            public List<Bound> Bounds { get; }
            public double EigenValue { get; set; }

            public int NumberSteps => Bounds.Count;

            public Result()
            {
                EigenVectors = new List<Vector<double>>();
                SchwartzConstants = new List<double>();
                SchwartzRelations = new List<double>();
                Bounds = new List<Bound>();
            }

        }

        public MainWindow()
        {
            InitializeComponent();
        }

        Result _kollattsTempleResult = new Result();
        Result _krylovBogolyubovResult = new Result();

        private void Solve(InitialСonditions initialСonditions)
        {
            _kollattsTempleResult = SolveKollattsTemple(initialСonditions);
            _krylovBogolyubovResult = SolveKrylovBogolyubov(initialСonditions);
        }

        private Result SolveKrylovBogolyubov(InitialСonditions initialСonditions)
        {
            Vector<double> xPartition = GetPartitionForIntegrationLimit(initialСonditions);
            Matrix<double> operatorMatrix = GetOperatorMatrix(xPartition, initialСonditions.OperatorCoeficientsFunction);

            var krylovBogolyubovResult = new Result();

            Vector<double> initialEigenVector = GetFunctionValuesOnPartition(xPartition, initialСonditions.InitialApproximation);
            Vector<double> operatorCoeficients = GetFunctionValuesOnPartition(xPartition, (x) => - 1.0 / initialСonditions.OperatorCoeficientsFunction(x));

            var schwartzConstant = GetSchwartzConstant(xPartition, initialEigenVector, initialEigenVector, operatorCoeficients);

            krylovBogolyubovResult.SchwartzConstants.Add(schwartzConstant);
            krylovBogolyubovResult.EigenVectors.Add(initialEigenVector);

            double krylovBogolyubovDifference = 2 * initialСonditions.Eps;
            int numberStep = 0;

            do
            {
                Vector<double> approxtimationVector = GetApproximationVector(krylovBogolyubovResult.EigenVectors[numberStep], /*operatorCoeficients, */initialСonditions.BoundaryConditions);

                Vector<double> nextEigenVector = TridiagonalMatrixAlgorithm(Matrix<double>.Build.DenseOfMatrix(operatorMatrix), Vector<double>.Build.DenseOfVector(approxtimationVector));

                krylovBogolyubovResult.EigenVectors.Add(nextEigenVector);

                var nextSchwartzConstant = GetSchwartzConstant(xPartition, initialEigenVector, nextEigenVector, operatorCoeficients);

                var nextSchwartzRelation = GetSchwartzRelation(krylovBogolyubovResult.SchwartzConstants.Last(), nextSchwartzConstant);

                krylovBogolyubovResult.SchwartzConstants.Add(nextSchwartzConstant);
                krylovBogolyubovResult.SchwartzRelations.Add(nextSchwartzRelation);

                int countSchwartzRelations = krylovBogolyubovResult.SchwartzRelations.Count;

                if (countSchwartzRelations > 1)
                {
                    double prevSchwartzRelation = krylovBogolyubovResult.SchwartzRelations[countSchwartzRelations - 2];

                    krylovBogolyubovDifference = 2 * Math.Sqrt((prevSchwartzRelation - nextSchwartzRelation) * nextSchwartzRelation);

                    krylovBogolyubovResult.Bounds.Add(new Bound()
                    {
                        LowerBound = nextSchwartzRelation - 0.5 * krylovBogolyubovDifference,
                        UpperBound = nextSchwartzRelation + 0.5 * krylovBogolyubovDifference
                    });
                }

                numberStep++;
            }
            while (Math.Abs(krylovBogolyubovDifference) > initialСonditions.Eps);

            krylovBogolyubovResult.EigenValue = krylovBogolyubovResult.SchwartzRelations.Last()
                                                + 0.5 * krylovBogolyubovDifference;

            return krylovBogolyubovResult;
        }

        private Result SolveKollattsTemple(InitialСonditions initialСonditions)
        {
            Vector<double> xPartition = GetPartitionForIntegrationLimit(initialСonditions);
            Matrix<double> operatorMatrix = GetOperatorMatrix(xPartition, initialСonditions.OperatorCoeficientsFunction);

            var kollattsTempleResult = new Result();

            Vector<double> initialEigenVector = GetFunctionValuesOnPartition(xPartition, initialСonditions.InitialApproximation);
            Vector<double> operatorCoeficients = GetFunctionValuesOnPartition(xPartition, (x) => -1.0 / initialСonditions.OperatorCoeficientsFunction(x));

            var schwartzConstant = GetSchwartzConstant(xPartition, initialEigenVector, initialEigenVector, operatorCoeficients);

            kollattsTempleResult.SchwartzConstants.Add(schwartzConstant);
            kollattsTempleResult.EigenVectors.Add(initialEigenVector);

            double kollattsTempleDifference = 2 * initialСonditions.Eps;
            int numberStep = 0;

            do
            {
                Vector<double> approxtimationVector = GetApproximationVector(kollattsTempleResult.EigenVectors[numberStep], /*operatorCoeficients,*/ initialСonditions.BoundaryConditions);

                Vector<double> nextEigenVector = TridiagonalMatrixAlgorithm(Matrix<double>.Build.DenseOfMatrix(operatorMatrix), Vector<double>.Build.DenseOfVector(approxtimationVector));

                kollattsTempleResult.EigenVectors.Add(nextEigenVector);

                var nextSchwartzConstant = GetSchwartzConstant(xPartition, initialEigenVector, nextEigenVector, operatorCoeficients);

                var nextSchwartzRelation = GetSchwartzRelation(kollattsTempleResult.SchwartzConstants.Last(), nextSchwartzConstant);
                kollattsTempleResult.SchwartzConstants.Add(nextSchwartzConstant);

                kollattsTempleResult.SchwartzRelations.Add(nextSchwartzRelation);

                int countSchwartzRelations = kollattsTempleResult.SchwartzRelations.Count;

                if (countSchwartzRelations > 1)
                {
                    double prevSchwartzRelation = kollattsTempleResult.SchwartzRelations[countSchwartzRelations - 2];

                    kollattsTempleDifference = (prevSchwartzRelation - nextSchwartzRelation) / ((initialСonditions.L2 / nextSchwartzRelation) - 1.0);

                    kollattsTempleResult.Bounds.Add(new Bound()
                    {
                        LowerBound = nextSchwartzRelation - kollattsTempleDifference,
                        UpperBound = nextSchwartzRelation
                    });
                }

                numberStep++;
            }
            while (Math.Abs(kollattsTempleDifference) > initialСonditions.Eps);

            kollattsTempleResult.EigenValue = kollattsTempleResult.SchwartzRelations.Last()
                                                - 0.5 * kollattsTempleDifference;

            return kollattsTempleResult;
        }

        private double GetSchwartzRelation(double prevSchwartzConstant, double nextSchwartzConstant)
        {
            return prevSchwartzConstant / nextSchwartzConstant;
        }

        private Vector<double> TridiagonalMatrixAlgorithm(Matrix<double> operatorMatrix, Vector<double> approxtimationVector)
        {
            int matrixDimension = approxtimationVector.Count;
            Vector<double> result = Vector<double>.Build.Dense(matrixDimension, 1);

            operatorMatrix[1, 0] /= operatorMatrix[0, 0];
            operatorMatrix[0, 0] /= operatorMatrix[0, 0];

            for (int i = 1; i < matrixDimension; i++)
            {
                double id = 1 / (operatorMatrix[i, i] - operatorMatrix[i, i - 1] * operatorMatrix[i - 1, i]);
                if (i != matrixDimension - 1)
                {
                    operatorMatrix[i + 1, i] *= id;
                }
                approxtimationVector[i] = (approxtimationVector[i] - approxtimationVector[i - 1] * operatorMatrix[i - 1, i]) * id;
            }

            result[matrixDimension - 1] = approxtimationVector[matrixDimension - 1];

            for (int i = matrixDimension - 2; i >= 0; i--)
            {
                result[i] = approxtimationVector[i] - operatorMatrix[i + 1, i] * result[i + 1];
            }

            return result;
        }

        private Matrix<double> GetOperatorMatrix(Vector<double> xPartitions, Func<double, double> function)
        {
            var matrixDimensional = xPartitions.Count;
            Matrix<double> operatorMatrix = Matrix<double>.Build.Dense(matrixDimensional, matrixDimensional);
            var step = (xPartitions[matrixDimensional - 1] - xPartitions[0]) / (matrixDimensional - 1);
            var squaredStep = step * step;
            operatorMatrix[0, 0] = function(xPartitions[0]);
            operatorMatrix[matrixDimensional - 1, matrixDimensional - 1] = function(xPartitions[matrixDimensional - 1]);
            for (int i = 1; i < matrixDimensional - 1; i++)
            {

                operatorMatrix[i - 1, i] = 1.0 / squaredStep * function(xPartitions[i]);
                operatorMatrix[i, i] = -2.0 / squaredStep * function(xPartitions[i]);
                operatorMatrix[i + 1, i] = 1.0 / squaredStep * function(xPartitions[i]);
            }

            return operatorMatrix;
        }

        private Vector<double> GetFunctionValuesOnPartition(Vector<double> partition, Func<double, double> function)
        {
            return Vector<double>.Build.Dense(partition.Count, i => function(partition[i]));
        }

        private Vector<double> GetApproximationVector(Vector<double> eigenVector, /*Vector<double> operatorCoeficiesnts,*/ BoundaryConditions boundaryConditions)
        {
            var approxtimationVector = Vector<double>.Build.Dense(
                eigenVector.Count,
                i => eigenVector[i]/* * operatorCoeficiesnts[i]*/);

            approxtimationVector[0] = boundaryConditions.LeftCondition;
            approxtimationVector[approxtimationVector.Count - 1] = boundaryConditions.RightCondition;

            return approxtimationVector;
        }

        private Vector<double> GetPartitionForIntegrationLimit(InitialСonditions initialСonditions)
        {
            var step = (initialСonditions.IntegrationLimits.UpperLimit - initialСonditions.IntegrationLimits.LowerLimit) / initialСonditions.NumberPartitions;
            return Vector<double>.Build.Dense(
                initialСonditions.NumberPartitions + 1,
                i => initialСonditions.IntegrationLimits.LowerLimit + i * step);
        }

        public double GetSchwartzConstant(Vector<double> xPartition, Vector<double> eigenVector1, Vector<double> eigenVector2, Vector<double> operatorCoeficiesnts)
        {
            Vector<double> subintegralVunstionVector = Vector<double>.Build.Dense(xPartition.Count, i => eigenVector1[i] * eigenVector2[i] * operatorCoeficiesnts[i]);
            return IntegrateFunctionVectorOnPartition(xPartition, subintegralVunstionVector);
        }

        private double IntegrateFunctionVectorOnPartition(Vector<double> xPartition, Vector<double> subintegralVunstionVector)
        {
            double valueOfIntegral = 0;

            for (var i = 1; i < xPartition.Count; i++)
            {
                valueOfIntegral += 0.5 * (subintegralVunstionVector[i] + subintegralVunstionVector[i - 1]) * (xPartition[i] - xPartition[i - 1]);
            }

            return valueOfIntegral;
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
                BoundaryConditions = new BoundaryConditions()
                {
                    LeftCondition = double.Parse(LeftConditionTextBox.Text),
                    RightCondition = double.Parse(RightConditionTextBox.Text)
                },
                InitialApproximation = x => 1,
                OperatorCoeficientsFunction = x => -Math.Pow(4 + x * x, 2),
                NumberPartitions = int.Parse(NumberPartitionsTextBox.Text),
                Eps = double.Parse(EpsTextBox.Text.Replace('.', ',')),
                L2 = 32 * Math.PI * Math.PI
            };

            Solve(initialСonditions);
            ShowKollattsTempleResult(_kollattsTempleResult, initialСonditions.L2);
            ShowKrylovBogolyubovResult(_krylovBogolyubovResult);
        }

        private void ShowKollattsTempleResult(Result resultKT, double l2)
        {
            IterationNumberKT.Content = resultKT.NumberSteps;
            ResultKT.Content = resultKT.EigenValue;
            SecondEigenValue.Content = l2;


            var constants = resultKT.SchwartzConstants.Select((x, i) => new
            {
                N = i,
                Constant = x
            }).ToList();

            schwartzValuesTableKT.ItemsSource = constants;

            var bounds = resultKT.Bounds.Select((x, i) => new
            {
                N = i,
                LowerBound = x.LowerBound,
                UpperBound = x.UpperBound
            }).ToList();

            BoundsTableKT.ItemsSource = bounds;
        }

        private void ShowKrylovBogolyubovResult(Result resultKB)
        {
            IterationNumberKB.Content = resultKB.NumberSteps;
            ResultKB.Content = resultKB.EigenValue;

            var constants = resultKB.SchwartzConstants.Select((x, i) => new
            {
                N = i,
                Constant = x
            }).ToList();

            schwartzValuesTableKB.ItemsSource = constants;

            var bounds = resultKB.Bounds.Select((x, i) => new
            {
                N = i,
                LowerBound = x.LowerBound,
                UpperBound = x.UpperBound
            }).ToList();

            BoundsTableKB.ItemsSource = bounds;
        }
    }
}
