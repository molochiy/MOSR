using System;
using System.Collections.Generic;
using System.Linq;
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

namespace WeinsteinMethod
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
        }

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

        private Vector<double> GetPartitionForIntegrationLimit(InitialСonditions initialСonditions)
        {
            var step = (initialСonditions.IntegrationLimits.UpperLimit - initialСonditions.IntegrationLimits.LowerLimit) / initialСonditions.NumberPartitions;
            return Vector<double>.Build.Dense(
                initialСonditions.NumberPartitions + 1,
                i => initialСonditions.IntegrationLimits.LowerLimit + i * step);
        }

        public static Vector<double> GetVectorOnPartition(Vector<double> xPartitions, int k)
        {
            Vector<double> vector = Vector<double>.Build.Dense(xPartitions.Count);

            var step = (xPartitions[xPartitions.Count - 1] - xPartitions[0]) / (xPartitions.Count - 1);

            return Vector<double>.Build.Dense(
                xPartitions.Count,
                i => Math.Sin((k + 1) * Math.PI * i * step));
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
                /*BoundaryConditions = new BoundaryConditions()
                {
                    LeftCondition = double.Parse(LeftConditionTextBox.Text),
                    RightCondition = double.Parse(RightConditionTextBox.Text)
                },*/
                InitialApproximation = x => 1,
                OperatorCoeficientsFunction = x => -Math.Pow(4 + x * x, 2),
                NumberPartitions = int.Parse(NumberPartitionsTextBox.Text),
                Eps = double.Parse(EpsTextBox.Text.Replace('.', ',')),
                L2 = 32 * Math.PI * Math.PI
            };

            var partition = GetPartitionForIntegrationLimit(initialСonditions);
            var matrixA = GetOperatorMatrix(partition, initialСonditions.OperatorCoeficientsFunction);
            var matrixB = GetOperatorMatrix(partition, x => -Math.Pow(3.9 + x * x, 2));
            var matrixC = matrixA - matrixB;

            var matricesA = new List<Matrix<double>>();
            var matricesC = new List<Matrix<double>>();
            var vectorsG = new List<Vector<double>>();

            var g0 = GetVectorOnPartition(partition, 0);

            vectorsG.Add((matrixC * g0).Normalize(g0));
            matricesC.Add(g0.ToScalarProductMatrix());

            matricesA.Add(matrixB + matricesC[0]);

            List<double> eigenValues = new List<double>();
            List<double> eigenValuesA = new List<double>();
            List<double> eigenValuesB = new List<double>();

            matricesA[0].Evd().EigenValues.Select(x => x.Real).ToList<double>();

            eigenValues.Add(matricesA[0].Evd().EigenValues.Select(x => x.Real).ToList<double>().Select(x => x).Where(x => x > 0).Min() - 41);
            eigenValuesA = matrixA.Evd().EigenValues.Select(x => x.Real).ToList<double>().Select(x => x).Where(x => x > 0).ToList<double>();
            eigenValuesB = matrixB.Evd().EigenValues.Select(x => x.Real).ToList<double>().Select(x => x).Where(x => x > 0).ToList<double>();
            eigenValuesA.Sort();
            eigenValuesB.Sort();

            for (int i = 1; i < initialСonditions.NumberPartitions; i++)
            {
                var gi = GetVectorOnPartition(partition, i);
                vectorsG.Add((matrixC * gi).Normalize(gi));
                matricesC.Add(matricesC[i - 1] + vectorsG[i].ToScalarProductMatrix());
                matricesA.Add(matrixB + matricesC[i]);
                eigenValues.Add(matricesA[i].Evd().EigenValues.Select(x => x.Real).ToList<double>().Select(x => x).Where(x => x > 0).Min() - 41);
                //.Select(x => x).Where(x => x > 2).Min() + 0.148);
            }
            eigenValues.Sort();
        }
    }

    public static class VectorExtensions
    {
        public static Vector<double> Normalize(this Vector<double> vector, Vector<double> basis)
        {
            double coef = 0;
            for (int i = 0; i < vector.Count; i++)
            {
                coef += vector[i] * basis[i];
            }

            return Vector<double>.Build.Dense(vector.Count, i => vector[i] / Math.Sqrt(coef));
        }

        public static Matrix<double> ToScalarProductMatrix(this Vector<double> vector)
        {
            Matrix<double> matrix = Matrix<double>.Build.Dense(vector.Count, vector.Count);

            for (int i = 0; i < vector.Count; i++)
            {
                for (int j = 0; j < vector.Count; j++)
                {
                    matrix[i, j] = vector[i] * vector[j];
                }
            }

            return matrix;
        }
    }
}
