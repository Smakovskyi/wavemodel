using System;

namespace wavemodel
{
    class Program
    {
        static void Main(/*string[] args*/)
        {
            PDESolver pdeSolver = new PDESolver();
            pdeSolver.Init(0.0003, 300, 300, 300, 300);
            pdeSolver.SetCoefficients(1000, 1500);
            double graphStep = 0.02;
            double graphTime = graphStep;

            while (pdeSolver.GettCurrent() < 2)
            {
                pdeSolver.CalcNextStep();
                Console.WriteLine(pdeSolver.GettCurrent());

                if (pdeSolver.GettCurrent() > graphTime)
                {
                    pdeSolver.SaveCurrentValues("outFixed" + graphTime + ".kr");
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
