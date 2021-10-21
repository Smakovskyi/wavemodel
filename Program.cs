using System;

namespace wavemodel
{
    class Program
    {
        static void Main(/*string[] args*/)
        {
            int Step = 10;
            int Nx = 300 / Step;
            int Ny = 300 / Step;
            int Nz = 200 / Step;

            float [] velocity = new float[Nz];
            
            for(int k = 0; k < Nz; k++)
                        velocity[k] = 1500 + k;

            PDESolver pdeSolver = new PDESolver();
            pdeSolver.Init(0.0001f, Nx, Ny, Nz);
            pdeSolver.SetCoefficients(1000, velocity);

            float graphStep = 0.01f;
            float graphTime = graphStep;


            while(pdeSolver.GettCurrent() < 0.5)
            {
                pdeSolver.CalcNextStep();
                pdeSolver.SaveToCSV();
                Console.WriteLine(pdeSolver.GettCurrent());

                if(pdeSolver.GettCurrent() > graphTime)
                {
                    pdeSolver.SaveCurrentValuesP("outFixedP0" + Math.Round(graphTime, 2) + ".kr");
                    graphTime += graphStep;
                    pdeSolver.FlushCSV();
                }
            }
            pdeSolver.Close();
        }
    }
}
