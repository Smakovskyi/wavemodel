
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
        double[,,] MurX;
        double[,,] MurY;
        double[,,] MurZ;

        double tCurrent;
        double dt;
        double dx;
        double dy;
        double dz;
        double ro;
        double velocity;

        StreamWriter datWriter;

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

            datWriter = new StreamWriter("P_all.csv", false, System.Text.Encoding.Default);

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


        public void Close()
        {
            datWriter.Close();
        }


        public void SaveToCSV()
        {
            datWriter.WriteLine( ("" + P[Nx/2, Ny/2, Nz-4]+";").Replace(',', '.'));
        }

        public void FlushCSV()
        {
            datWriter.Flush();
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
            MurX = new double[4 ,Ny ,Nz];
            for (int j = 0; j < Ny; j++) 
            {
                for (int k = 0; j < Nz; j++)
                {
                    MurX[0, j, k] = 0.0;
                    MurX[1, j, k] = 0.0;
                    MurX[2, j, k] = 0.0;
                    MurX[3, j, k] = 0.0;
                }
            }
            MurY = new double[Nx, 4, Nz];
            for(int i=0; i<Nx; i++)
            {
                for(int k=0; k<Nz; k++)
                {
                    MurY[i, 0,  k] = 0.0;
                    MurY[i, 1,  k] = 0.0;
                    MurY[i, 2,  k] = 0.0;
                    MurY[i, 3,  k] = 0.0;
                }
            }
            MurZ = new double[Nx, Ny, 4];
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    MurZ[i, j, 0] = 0.0;
                    MurZ[i, j, 1] = 0.0;
                    MurZ[i, j, 2] = 0.0;
                    MurZ[i, j, 3] = 0.0;
                }
            }
        }

        private void MurBoundaries()
        {
            double reflectionCoefficient = 0.5;
            // x==0 boundary, partial reflection 
            for (int j = 1; j < Ny - 1; j++)
            {
                for (int k = 0; k < Nz; k++)
                {
                    P[0, j, k] = reflectionCoefficient * P[0, j, k] +
                              (1 - reflectionCoefficient) * (MurX[1, j, k] + (velocity * dt - dx)
                                             / (velocity * dt + dx) * (P[1, j, k] - MurX[0, j, k]));
                }
            }
            // x == Lx, Mur absorbing
            for (int j = 1; j < Ny - 1; j++)
            {
                for (int k = 0; k < Nz; k++)
                {
                    P[Nx - 1, j, k] = MurX[2, j, k] + (velocity * dt - dx) / (velocity * dt + dx) * (P[Nx - 2, j, k] - MurX[3, j, k]);
                }
            }
            // y == 0
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int k = 0; k < Nz; k++)                     // 
                {                                                // dx - maybe dy ?????             
                    P[i, 0, k] = MurY[i, 1, k] + (velocity * dt - dx) / (velocity * dt + dx) * (P[i, 1, k] - MurY[i, 0, k]);
                }
            }
            // y == Ly
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int k = 0; k < Nz; k++)
                {                                               //dx
                    P[i, Ny - 1, k] = MurY[i, 2, k] + (velocity * dt - dx) / (velocity * dt + dx) * (P[i, Ny - 2, k] - MurY[i, 3, k]);
                }
            }
            // z == 0
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int j = 1; j < Ny - 1; j++)
                {                                                 //dx
                    P[i, j, 0] = MurZ[i, j, 1] + (velocity * dt - dx) / (velocity * dt + dx) * (P[i, j, 1] - MurZ[i, j, 0]);
                }

            }
            // z == Lz
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int j = 1; j < Ny - 1; j++)
                {                                                    ///????       dx
                    P[i, j, Nz - 1] = MurZ[i, j, 2] + (velocity * dt - dx) / (velocity * dt + dx) * (P[i, j ,Nz - 2] - MurZ[i, j, 3]);
                }

            }
            //
            MurCopy();
        }


        private void MurCopy()
        {
            for (int j = 0; j < Ny; j++)
            {
                for (int k = 0; k < Nz; k++)
                {
                    MurX[0, j, k] = P[0, j, k];
                    MurX[1, j, k] = P[1, j, k];
                    MurX[2, j, k] = P[Nx - 2, j, k];
                    MurX[3, j, k] = P[Nx - 1, j, k];
                }
            }

            for (int i = 0; i < Nx; i++)
            {
                for (int k = 0; k < Nz; k++)
                {
                    MurY[i,0,k] = P[i, 0, k];
                    MurY[i,1,k] = P[i, 1, k];
                    MurY[i,2,k] = P[i, Ny - 2, k];
                    MurY[i,3,k] = P[i, Ny - 1, k];
                }
            }
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j< Ny; j++)
                {
                    MurZ[i, j, 0] = P[i, j, 0];
                    MurZ[i, j, 1] = P[i, j, 1];
                    MurZ[i, j, 2] = P[i, j, Nz - 2];
                    MurZ[i, j, 3] = P[i, j, Nz - 1];
                }
            }
        }
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
            double dt_dz_ro = (dt / dz) * ro * velocity * velocity;

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
