using System;
using System.IO;
namespace wavemodel
{
    public class PDESolver
    {
        double[,] uCurrent;
        double[,] uPrevious;
        double[,] uNext;
        double tCurrent;
        double dt;
        double dx;
        double dy;
        double ax;
        double ay;
        double lX;
        double lY;
        int Nx;
        int Ny;

        public double muX(double y, double t)
        {
            return 0.0;
        }

        public double nuX(double y, double t)
        {
            return 0.0;
        }

        public double muY(double x, double t)
        {
            if (x < 0.4 * lX && x >= 0.2 * lX)
            {
                return Math.Cos(200 * t);
            }
            return 0.0;
        }

        public double nuY(double x, double t)
        {
            return 0.0;
        }

        public double fi(double x, double y)
        {
            return 0.0;
        }

        public double psi(double x, double y)
        {
            return 0.0;
        }

        double f(double x, double y, double t)
        {
            /*if (x < 0.6 * lX && x >= 0.4 * lX
                && y < 0.6 * lY && y >= 0.4 * lY)
            {
                return Math.Cos(20 * t);
            }*/
            return 0;
        }

        double getAy(double x)
        {
            return ay*(1 - 0.3 * x/lX);
        }
        public void init(double dt, double lX, double lY, int Nx, int Ny)
        {
            this.Nx = Nx;
            this.Ny = Ny;
            this.dt = dt;
            this.lX = lX;
            this.dx = lX / Nx;
            this.lY = lY;
            this.dy = lY / Ny;
            initGrid();
            tCurrent = dt;
        }

        public void setCoefficients(double ax, double ay)
        {
            this.ax = ax;
            this.ay = ay;
        }

        void initGrid()
        {
            uCurrent = new double[Nx + 1,Ny + 1];
            uPrevious = new double[Nx + 1,Ny + 1];
            uNext = new double[Nx + 1,Ny + 1];

            for (int i = 0; i <= Nx; i++)
            {
                for (int j = 0; j <= Ny; j++)
                {
                    uPrevious[i,j] = fi(i * dx, j * dy);
                    uCurrent[i,j] = uPrevious[i,j] + psi(i * dx, j * dy) * dt;
                }
            }
            tCurrent = dt;
        }

        public void calcNextStep()
        {
            fillBoundaries();
            calculateInnerValues();
            tCurrent += dt;
            reassignArrays();
        }

        private void reassignArrays()
        {
            double[,] temp;
            temp = uPrevious;
            uPrevious = uCurrent;
            uCurrent = uNext;
            uNext = temp;
        }

        private void calculateInnerValues()
        {
            double gammaX = ax * ax * dt * dt / (dx * dx);
            //double gammaY = ay * ay * dt * dt / (dy * dy);
            for (int i = 1; i < Nx; i++)
            {
                double ay = getAy(i*dx);
                double gammaY = ay * ay * dt * dt / (dy * dy);
                for (int j = 1; j < Ny; j++)
                {
                    
                    uNext[i,j] = (2 - 2 * gammaX - 2 * gammaY) * uCurrent[i,j]
                    - uPrevious[i,j] + gammaX * (uCurrent[i - 1,j] + uCurrent[i + 1,j])
                    + gammaY * (uCurrent[i,j - 1] + uCurrent[i,j + 1])
                    + dt * dt * f(i * dx, j * dy, tCurrent);
                }
            }
        }

        private void fillBoundaries()
        {
            for (int i = 0; i <= Nx; i++)
            {
                double x = dx * i;
                uNext[i,0] = muY(x, tCurrent);
                uNext[i,Ny] = nuY(x, tCurrent);
            }

            for (int j = 1; j < Ny; j++)
            {
                double y = dy * j;
                uNext[0,j] = muX(y, tCurrent);
                uNext[Nx,j] = nuX(y, tCurrent);
            }

        }

        public void saveCurrentValues(String fileName)
        {
            using(StreamWriter outWriter = 
                new StreamWriter(fileName, false, System.Text.Encoding.Default))
            {
                
                for (int i = 0; i <= Nx; i++)
                {
                    for (int j = 0; j <= Ny; j+=10)
                    {
                        outWriter.WriteLine(dx * i + " " + dy * j + " " + uCurrent[i,j]);
                    }

                }

            }
        }

        public double gettCurrent()
        {
            return tCurrent;
        }
    }
}
