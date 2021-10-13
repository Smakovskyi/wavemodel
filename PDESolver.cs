
using System;
using System.IO;

namespace wavemodel
{
    public class PDESolver
    {
        #region var

        double[,,] Vx;
        double[,,] Vy;
        double[,,] Vz;
        double[,,] P;
        double[,,] Mur_X1;

        double tCurrent;
        double dt;
        double dx;
        double dy;
        double dz;
        double ro;
        double velocity;
        
        int Nx, Ny, Nz;

        #endregion

        #region unused 1

        public double Psi(double x, double y)
        {
            return 0.0;
        }

        #endregion

        #region init

        public void Init(double dt, double lX, double lY, double lZ, int Nx, int Ny, int Nz)
        {
            this.Nx = Nx;
            this.Ny = Ny;
            this.Nz = Nz;
            this.dt = dt;

            this.dx = lX / Nx;
            this.dy = lY / Ny;
            this.dz = lZ / Nz;

            InitGrid();
        }
        void InitGrid()
        {
            Vx = new double[Nx + 1, Ny, Nz];
            Vy = new double[Nx, Ny + 1, Nz];
            Vz = new double[Nx, Ny, Nz + 1];
            P = new double[Nx, Ny, Nz];

            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for(int k = 0; k < Nz; k++)
                    {
                        P[i, j, k] = 0;
                    }
                    
                }
            }
            for (int i = 0; i <= Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        Vx[i, j, k] = 0;
                    }
                }
            }
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j <= Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        Vy[i, j, k] = 0;
                    }
                }
            }
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k <= Nz; k++)
                    {
                        Vz[i, j, k] = 0;
                    }
                }
            }
            tCurrent = 0;
        }

        double F(int i, int j, int  k, double t)
        {
            if ( (Math.Abs(i - Nx / 2) <= 2) 
                 && (Math.Abs(j - Ny / 2) <= 2)
                 && (Math.Abs(k - Nz / 2) <= 2) 
                 && t < 0.2)
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
            MurBoundaries();   
            tCurrent += dt;
        }

        public void InitMur1st()
        {
            Mur_X1 = new double[4 ,Ny ,Nz];
            for (int j = 0; j < Ny; j++) 
            {
                for (int k = 0; j < Nz; j++)
                {
                    Mur_X1[0, j, k] = 0.0;
                    Mur_X1[1, j, k] = 0.0;
                    //Mur_X1[2,j] = 0.0;
                    //Mur_X1[3,j] = 0.0;
                }
            }
        }

        private void MurBoundaries()
        {
            double reflectionCoefficient = 0.2;

            for (int j = 1; j < Ny - 1; j++) 
            { 
                for (int k = 0; k < Nz; k++) 
                {

                        P[0, j, k] = reflectionCoefficient * P[0, j, k] +
                                  (1 - reflectionCoefficient) * (Mur_X1[1, j, k] + (velocity * dt - dx)
                                                 / (velocity * dt + dx) * (P[1, j, k] - Mur_X1[0, j, k]));
                }
            }

            for (int j = 0; j < Ny; j++)
            {
                for (int k = 0; k < Nz; k++) 
                {
                    Mur_X1[0, j, k] = P[0, j, k];
                    Mur_X1[1, j, k] = P[1, j, k];
                }
            }
        }
        
        //

        private void UpdateV()
        {
            // dt_over_dx = (dt / dy) / ro;
            // dt_over_dy = (dt / dy) / ro;

            double dt_dx_ro = dt / (dx * ro);
            double dt_dy_ro = dt / (dy * ro);
            double dt_dz_ro = dt / (dz * ro);

            for (int i = 1; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        Vx[i, j, k] -= dt_dx_ro * (P[i, j, k] - P[i - 1, j , k]);
                    }
                }
            }
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 1; j < Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        Vy[i, j, k] -= dt_dy_ro * (P[i, j, k] - P[i, j - 1, k]);
                    }
                }
            }

            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 1; k < Nz; k++)
                    {
                        Vz[i, j, k] -= dt_dz_ro * (P[i, j, k] - P[i , j, k - 1]);
                    }
                }
            }
        }

        private void UpdateP()
        {
            

            double dt_dx_ro = (dt / dx) * ro * velocity * velocity;
            double dt_dy_ro = (dt / dy) * ro * velocity * velocity;
            double dt_dz_ro = dt / (dz * ro);

            for (int i = 0; i < Nx; i++)
                for(int j = 0; j < Ny; j++)
                    for (int k = 0; k < Nz; k++)
                    {
                        P[i, j, k] -=
                        dt_dx_ro * (Vx[i + 1, j, k] - Vx[i, j, k]) +
                        dt_dy_ro * (Vy[i, j + 1, k] - Vy[i, j, k]) +
                        dt_dz_ro * (Vz[i, j , k + 1] - Vz[i, j, k])
                        - dt * this.F(i, j, k, tCurrent);
                    }
        }

        #region save

        public void SaveCurrentValuesP(String fileName)
        {
            using StreamWriter outWriter = new StreamWriter(fileName, false, System.Text.Encoding.Default);

            for(int i = 0; i < Nx; i ++)
                for(int j = 0; j < Ny; j ++)
                    outWriter.WriteLine((dx * i + " " + dy * j + " " + P[i, j , Nz / 2]).Replace(',', '.'));
        }
        /*
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
        }*/
        public double GettCurrent()
        {
            return tCurrent;
        }

        #endregion
    }
}
