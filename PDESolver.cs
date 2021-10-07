using System;
using System.IO;
namespace wavemodel
{
    public class PDESolver
    {
        double[,] Vx;
        double[,] Vy;
        double[,] P;

        double tCurrent;
        double dt;
        double dx;
        double dy;
        double ro;
        double velocity;
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

        double f(int i, int j, double t)
        {
            if (i == Nx / 2 && j == Ny / 2 && t < 0.2)
            {
                return  Math.Cos(Math.PI * 2 * t);
            }
            return 0;
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

        public void setCoefficients(double ro, double velocity)
        {
            this.ro = ro;
            this.velocity = velocity;
        }

        void initGrid()
        {

            Vx = new double[Nx + 1, Ny + 1];
            Vy = new double[Nx + 1, Ny + 1];
            P = new double[Nx + 1, Ny + 1];

            for (int i = 0; i <= Nx; i++)
            {
                for (int j = 0; j <= Ny; j++)
                {
                    Vx[i, j] = 0;
                    Vy[i, j] = 0;
                    P[i, j] = 0;
                }
            }
            tCurrent = 0;
        }

        public void calcNextStep()
        {
            fillBoundaries();
            updateV();
            updateP();
            tCurrent += dt;
        }

       

        private void updateV()
        {
            double dt_over_dx = dt / dx;
            double dt_over_dy = dt / dy;


            for (int i = 1; i < Nx; i++)
            {
                for (int j = 1; j < Ny; j++)
                {
                   Vx[i,j] += - dt_over_dx / ro * (P[i,j] - P[i - 1,j]);
                }
            }

            for (int i = 1; i < Nx; i++)
            {
                for (int j = 1; j < Ny; j++)
                {
                    Vy[i, j] += - dt_over_dy / ro * (P[i, j] - P[i, j - 1]);
                }
            }
        }

        private void updateP()
        {
            double dt_ro_c2 = dt * ro * velocity *10/** velocity*/;
            for (int i = 1; i < Nx; i++)
            {
                
                for (int j = 1; j < Ny; j++)
                {

                    P[i,j] += -dt_ro_c2 * ((Vx[i+1,j] - Vx[i,j])/dx  + (Vy[i, j+1] - Vy[i, j]))/dy + dt*this.f(i , j , tCurrent);
                }
            }
        }

        private void fillBoundaries()
        {
            for (int i = 0; i <= Nx; i++)
            {
                double x = dx * i;
                Vx[i,0] = muY(x, tCurrent);
                Vx[i,Ny] = nuY(x, tCurrent);
                Vy[i, 0] = muY(x, tCurrent);
                Vy[i, Ny] = nuY(x, tCurrent);
                P[i, 0] = muY(x, tCurrent);
                P[i, Ny] = nuY(x, tCurrent);
            }

            for (int j = 1; j < Ny; j++)
            {
                double y = dy * j;
                Vx[0,j] = muX(y, tCurrent);
                Vx[Nx,j] = nuX(y, tCurrent);
                Vy[0, j] = muX(y, tCurrent);
                Vy[Nx, j] = nuX(y, tCurrent);
                P[0, j] = muX(y, tCurrent);
                P[Nx, j] = nuX(y, tCurrent);
            }

        }

        public void saveCurrentValues(String fileName)
        {
            using(StreamWriter outWriter = 
                new StreamWriter(fileName, false, System.Text.Encoding.Default))
            {
                
                for (int i = 0; i <= Nx; i+=1)
                {
                    for (int j = 0; j <= Ny; j+=1)
                    {
                        outWriter.WriteLine((dx * i + " " + dy * j + " " + Math.Abs(P[i,j])).Replace(',','.'));
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
