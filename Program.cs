using System;

namespace wavemodel
{
    class Program
    {
        static void Main(string[] args)
        {
            PDESolver pdeSolver = new PDESolver();
            pdeSolver.init(0.01, 1, 10, 100, 1000);
            pdeSolver.setCoefficients(0.5, 0.5);
            while (pdeSolver.gettCurrent() < 50)
            {
                pdeSolver.calcNextStep();
            }
            pdeSolver.saveCurrentValues("out.kr");
        }
    }
}
