using System;

namespace wavemodel
{
    class Program
    {
        static void Main(/*string[] args*/)
        {
            PDESolver pdeSolver = new PDESolver();
            pdeSolver.Init(0.0003, 300.0 , 300.0 ,300.0, 60, 60, 60);

            pdeSolver.SetCoefficients(1000, 1500);
            pdeSolver.InitMur1st();
            double graphStep = 0.01;
            double graphTime = graphStep;

            while (pdeSolver.GettCurrent() < 0.5)
            {
                pdeSolver.CalcNextStep();
                pdeSolver.SaveToCSV();
                Console.WriteLine(pdeSolver.GettCurrent());

                if (pdeSolver.GettCurrent() > graphTime)
                {
                    pdeSolver.SaveCurrentValuesP("outFixedP0" +  Math.Round(graphTime, 2) + ".kr");
                    // pdeSolver.SaveCurrentValuesVx("outFixedVx" + Math.Round(graphTime, 2) + ".kr");
                    // pdeSolver.SaveCurrentValuesVy("outFixedVy" + Math.Round(graphTime, 2) + ".kr");
                    graphTime += graphStep;
                    pdeSolver.FlushCSV();
                }
            }

            pdeSolver.Close();

            /*
            while (pdeSolver.GettCurrent() < 100)
            {
                pdeSolver.CalcNextStep();
            }
            pdeSolver.SaveCurrentValues("out.kr");
            */
        }
    }
}
