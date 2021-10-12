
using System;
using System.IO;

namespace wavemodel
{
    public class PDESolver
    {
        #region var

        double[,] Vx;
        double[,] Vy;
        double[,] P;

        double tCurrent;
        double dt;
        double dx;
        double dy;
        double ro;
        double velocity;
        // double lX;
        // double lY;
        int Nx;
        int Ny;

        #endregion

        #region unused 1

        //x = 0 border
        static double MuX(double y, double t)
        {
            return 0.0;
        }
        //x = x_max border
        public double NuX(double y, double t)
        {
            return 0.0;
        }
        //y = 0
        public double MuY(double x, double t)
        {
            return 0.0;
        }
        //y = y_max
        public double NuY(double x, double t)
        {
            return 0.0;
        }
        public double Fi(double x, double y)
        {
            return 0.0;
        }
        public double Psi(double x, double y)
        {
            return 0.0;
        }

        #endregion

        #region init

        public void Init(double dt, double lX, double lY, int Nx, int Ny)
        {
            this.Nx = Nx;
            this.Ny = Ny;
            this.dt = dt;

            // this.lX = lX;
            this.dx = lX / Nx;
            // this.lY = lY;
            this.dy = lY / Ny;

            InitGrid();
        }
        void InitGrid()
        {
            Vx = new double[Nx + 1, Ny];
            Vy = new double[Nx, Ny + 1];
            P = new double[Nx, Ny];

            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    P[i, j] = 0;
                }
            }
            for (int i = 0; i <= Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    Vx[i, j] = 0;
                }
            }
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j <= Ny; j++)
                {
                    Vy[i, j] = 0;
                }
            }
            tCurrent = 0;
        }

        double F(int i, int j, double t)
        {
            if ((Math.Abs(i - Nx / 2) <= 2) && (Math.Abs(j - Ny / 2) <= 2) && t < 0.2)
            {
                return Math.Cos(Math.PI * 2 * t);
            }
            return 0;
        }
        public void SetCoefficients(double ro, double velocity)
        {
            this.ro = ro;
            this.velocity = velocity;
        }

        #endregion

        public void CalcNextStep()
        {

            UpdateV();
            //FillBoundariesV();
            UpdateP();
            //FillBoundariesP();
            
            tCurrent += dt;
        }

        private void MurBoundaries()
        {

        }

        public void FillBoundariesV()
        {
            for (int j = 0; j < Ny; j++)
            {
                for (int i = 0; i < 10; i++)
                {
                    Vx[i, j] *= Math.Exp(-(10-i));
                    Vy[i, j] *= Math.Exp(-(10-i));
                    P[i, j] *= Math.Exp(-(10 - i));
                }
            }
        }


        public void FillBoundariesP()
        {
            for(int j=0; j <= Ny; j++)
            {
                P[0, j] = P[1, j];
            }            
        }
        private void FillBoundaries()
        {
            //for(int i = 0; i <= Nx; i++)
            //{
            //     double x = dx * i;

            //    Vx[i, 0] = MuY(x, tCurrent);
            //    Vy[i, 0] = MuY(x, tCurrent);
            //    P[i, 0] = MuY(x, tCurrent);
            //    P[i, 1] = MuY(x, tCurrent);
            //    Vx[i, 1] = MuY(x, tCurrent);
            //    Vy[i, 1] = MuY(x, tCurrent);

            //    //Vx[i, Ny] = NuY(x, tCurrent);
            //    //Vy[i, Ny] = NuY(x, tCurrent);
            //    Vx[i, Ny] = 0;
            //    Vy[i, Ny] = 0;
            //    P[i, Ny] = P[i, Ny-1];
            //}

            //for (int j = 0; j <= Ny; j++)
            //{
            //     double y = dy * j;

            //    Vx[0, j] = MuX(y, tCurrent);
            //    Vx[Nx, j] = NuX(y, tCurrent);

            //    Vy[0, j] = MuX(y, tCurrent);
            //    Vy[Nx, j] = NuX(y, tCurrent);

            //    P[0, j] = MuX(y, tCurrent);
            //    P[Nx, j] = NuX(y, tCurrent);
            //}
           /* for (int j = 0; j <= Ny; j++)
            {
                P[0, j] = 0;
                P[1, j] = 0;
                Vx[0, j] = 0;
                Vx[1, j] = 0;
                Vy[0, j] = 0;
                Vy[1, j] = 0;
            }*/

        }

        //

        private void UpdateV()
        {
            // dt_over_dx = (dt / dy) / ro;
            // dt_over_dy = (dt / dy) / ro;

            double dt_dx_ro = dt / (dx * ro);
            double dt_dy_ro = dt / (dy * ro);

            for (int i = 1; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    Vx[i, j] -= dt_dx_ro * (P[i, j] - P[i - 1, j]);
                }
            }
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 1; j < Ny; j++)
                {
                    Vy[i, j] -= dt_dy_ro * (P[i, j] - P[i, j - 1]);
                }
            }
        }

        private void UpdateP()
        {
            //double dt_ro_c2x = dt * ro * velocity * 10 / dx; /*velocity*/
            //double dt_ro_c2y = dt * ro * velocity * 10 / dy; /*velocity*/

            double dt_dx_ro = (dt / dx) * ro * velocity * velocity;
            double dt_dy_ro = (dt / dy) * ro * velocity * velocity;

            for(int i = 0; i < Nx; i++)
                for(int j = 0; j < Ny; j++)
                {
                    P[i, j] -=
                        dt_dx_ro * (Vx[i + 1, j] - Vx[i, j]) +
                        dt_dy_ro * (Vy[i, j + 1] - Vy[i, j]) -
                        dt * this.F(i, j, tCurrent);
                }
        }

        #region save

        public void SaveCurrentValuesP(String fileName)
        {
            using StreamWriter outWriter = new StreamWriter(fileName, false, System.Text.Encoding.Default);

            for(int i = 0; i < Nx; i ++)
                for(int j = 0; j < Ny; j ++)
                    outWriter.WriteLine((dx * i + " " + dy * j + " " + P[i, j]).Replace(',', '.'));
        }
        public void SaveCurrentValuesVx(String fileName)
        {
            using StreamWriter outWriter = new StreamWriter(fileName, false, System.Text.Encoding.Default);

            for (int i = 0; i <= Nx; i += 1)
                for (int j = 0; j < Ny; j += 1)
                    outWriter.WriteLine((dx * i + " " + dy * j + " " + Vx[i, j]).Replace(',', '.'));
        }
        public void SaveCurrentValuesVy(String fileName)
        {
            using StreamWriter outWriter = new StreamWriter(fileName, false, System.Text.Encoding.Default);

            for (int i = 0; i < Nx; i += 1)
                for (int j = 0; j <= Ny; j += 1)
                    outWriter.WriteLine((dx * i + " " + dy * j + " " + Vy[i, j]).Replace(',', '.'));
        }
        public double GettCurrent()
        {
            return tCurrent;
        }

        #endregion
    }
}
