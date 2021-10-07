using System;

namespace wavemodel
{
    class Program
    {
        static void Main(string[] args)
        {
            PDESolver pdeSolver = new PDESolver();
            pdeSolver.init(0.001, 500, 500, 500, 500);
            pdeSolver.setCoefficients(1000, 1500);
            double graphStep = 0.01;
            double graphTime = graphStep;
            while (pdeSolver.gettCurrent() < 1)
            {
                pdeSolver.calcNextStep();
                Console.WriteLine(pdeSolver.gettCurrent());
                if (pdeSolver.gettCurrent() > graphTime)
                {
                    pdeSolver.saveCurrentValues("outFixed" + graphTime + ".kr");
                    
                    graphTime += graphStep;
                }
            }
            while (pdeSolver.gettCurrent() < 50)
            {
                pdeSolver.calcNextStep();
            }
            pdeSolver.saveCurrentValues("out.kr");
        }
    }
}
