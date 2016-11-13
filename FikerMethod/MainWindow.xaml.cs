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
using MathNet.Numerics.LinearAlgebra;

namespace FikerMethod
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
        public int NumberPartitions { get; set; }
        public double Eps { get; set; }
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
                NumberPartitions = int.Parse(NumberPartitionsTextBox.Text),
            };

            var partition = GetPartitionForIntegrationLimit(initialСonditions);
            var operatorMatrix = GetOperatorMatrix(partition, initialСonditions.OperatorCoeficientsFunction);
            var gMatrix = operatorMatrix.Inverse();
            var eigen = gMatrix.Evd();
            var allEigenValues = eigen.EigenValues.Map(c => 1.0 / c.Real);
            var eigenValues = allEigenValues.Where(x => x > 0).Select((x, i) => new
            {
                N = i + 1,
                EigenValue = x
            });

            EigenValuesTable.ItemsSource = eigenValues;
        }
    }
}
