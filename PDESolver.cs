
using System;
using System.IO;

namespace wavemodel
{
    public class PDESolver
    {
        #region var

        double[,] Vx;
        double[,] Vy;
        double[,] P0;

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

        static double MuX(/*double y, double t*/)
        {
            return 0.0;
        }
        public double NuX(/*double y, double t*/)
        {
            return 0.0;
        }
        public double MuY(/*double x, double t*/)
        {
            return 0.0;
        }
        public double NuY(/*double x, double t*/)
        {
            return 0.0;
        }
        public double Fi(/*double x, double y*/)
        {
            return 0.0;
        }
        public double Psi(/*double x, double y*/)
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
            Vx = new double[Nx + 1, Ny + 1];
            Vy = new double[Nx + 1, Ny + 1];
            P0 = new double[Nx + 1, Ny + 1];

            for(int i = 0; i <= Nx; i++)
                for(int j = 0; j <= Ny; j++)
                {
                    Vx[i, j] = 0;
                    Vy[i, j] = 0;
                    P0[i, j] = 0;
                }
            tCurrent = 0;
        }

        double F(int i, int j, double t)
        {
            if( (Math.Abs( i - Nx / 2) <= 2) && (Math.Abs(j - Ny / 2) <= 2) && t < 0.2)
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
            FillBoundaries();
            UpdateV();
            UpdateP();

            tCurrent += dt;
        }
        private void FillBoundaries()
        {
            for(int i = 0; i <= Nx; i++)
            {
                // double x = dx * i;

                Vx[i, 0] = MuY(/*x, tCurrent*/);
                Vx[i, Ny] = NuY(/*x, tCurrent*/);

                Vy[i, 0] = MuY(/*x, tCurrent*/);
                Vy[i, Ny] = NuY(/*x, tCurrent*/);

                P0[i, 0] = MuY(/*x, tCurrent*/);
                P0[i, Ny] = NuY(/*x, tCurrent*/);
            }

            for(int j = 0; j <= Ny; j++)
            {
                // double y = dy * j;

                Vx[0, j] = MuX(/*y, tCurrent*/);
                Vx[Nx, j] = NuX(/*y, tCurrent*/);

                Vy[0, j] = MuX(/*y, tCurrent*/);
                Vy[Nx, j] = NuX(/*y, tCurrent*/);

                P0[0, j] = MuX(/*y, tCurrent*/);
                P0[Nx, j] = NuX(/*y, tCurrent*/);
            }
        }

        //

        private void UpdateV()
        {
            // dt_over_dx = (dt / dy) / ro;
            // dt_over_dy = (dt / dy) / ro;

            double dt_dx_ro = dt / (dx * ro);
            double dt_dy_ro = (dt / dy) * (1 / (ro * 1));
                
            for(int i = 1; i < Nx; i++)
                for(int j = 1; j < Ny; j++)
                    Vx[i, j] -= dt_dx_ro * (P0[i, j] -   P0[i - 1, j]);

            for(int i = 1; i < Nx; i++)
                for(int j = 1; j < Ny; j++)
                    Vy[i, j] -= dt_dy_ro * (P0[i, j] -  P0[i, j - 1]);
        }

        private void UpdateP()
        {
            //double dt_ro_c2x = dt * ro * velocity * 10 / dx; /*velocity*/
            //double dt_ro_c2y = dt * ro * velocity * 10 / dy; /*velocity*/

            double dt_dx_ro = (dt / dx) * ro * velocity * velocity;
            double dt_dy_ro = (dt / dy) * ro * velocity * velocity;

            for(int i = 1; i < Nx; i++)
                for(int j = 1; j < Ny; j++)
                {
                    P0[i, j] -=
                        dt_dx_ro * (Vx[i + 1, j] - Vx[i, j]) +
                        dt_dy_ro * (Vy[i, j + 1] - Vy[i, j]) -
                        dt * this.F(i, j, tCurrent);
                }
        }

        #region save

        public void SaveCurrentValues(String fileName)
        {
            using StreamWriter outWriter = new StreamWriter(fileName, false, System.Text.Encoding.Default);

            for(int i = 0; i <= Nx; i += 1)
                for(int j = 0; j <= Ny; j += 1)
                    outWriter.WriteLine((dx * i + " " + dy * j + " " + Math.Abs(P0[i, j])).Replace(',', '.'));
        }
        public double GettCurrent()
        {
            return tCurrent;
        }

        #endregion
    }
}
