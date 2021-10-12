using System;

namespace wavemodel
{
    class Program
    {
        static void Main(/*string[] args*/)
        {
            PDESolver pdeSolver = new PDESolver();
            pdeSolver.SetCoefficients(1000, 1500);
            pdeSolver.Init(0.0005f, 3000, 3000, 300, 300);

            float graphStep = 0.01f;
            float graphTime = graphStep;

            while (pdeSolver.GettCurrent() < 0.5)
            {
                pdeSolver.CalcNextStep();
                Console.WriteLine(pdeSolver.GettCurrent());

                if (pdeSolver.GettCurrent() > graphTime)
                {
                    pdeSolver.SaveCurrentValuesP0("outFixedP0" + Math.Round(graphTime, 2) + ".kr");
                    // pdeSolver.SaveCurrentValuesVx("outFixedVx" + Math.Round(graphTime, 2) + ".kr");
                    // pdeSolver.SaveCurrentValuesVy("outFixedVy" + Math.Round(graphTime, 2) + ".kr");

                    graphTime += graphStep;
                }
            }

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
