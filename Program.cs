using System;

namespace wavemodel
{
    class Program
    {
        static void Main(/*string[] args*/)
        {
            int Nx = 300;
            int Ny = 300;
            int Nz = 300;

            float [,,] velocity = new float[Nx, Ny, Nz];
            for(int i = 0; i < Nz; i++)
                for(int j = 0; j < Nz; j++)
                    for(int k = 0; k < Nz; k++)
                        velocity[i, j, k] = 1500 + k;

            PDESolver pdeSolver = new PDESolver();
            pdeSolver.Init(0.0003f, 300, 300, 300, Nx, Ny, Nz);
            pdeSolver.SetCoefficients(1000, velocity);

            double graphStep = 0.01;
            double graphTime = graphStep;

            while(pdeSolver.GettCurrent() < 0.5)
            {
                pdeSolver.CalcNextStep();
                pdeSolver.SaveToCSV();
                Console.WriteLine(pdeSolver.GettCurrent());

                if(pdeSolver.GettCurrent() > graphTime)
                {
                    pdeSolver.SaveCurrentValuesP("outFixedP0" + Math.Round(graphTime, 2) + ".kr");
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
