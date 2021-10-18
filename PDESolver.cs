
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
        double[,,] P0;

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

/*
        //x = 0 border
        static long MuX(int x, float t)
        {
            return 0;
        }
        //x = x_max border
        public long NuX(int x, float t)
        {
            return 0;
        }
        //y = 0
        public long MuY(int y, float t)
        {
            return 0;
        }
        //y = y_max
        public long NuY(int y, float t)
        {
            return 0;
        }

        public long Fi(int x, int y)
        {
            return 0;
        }
*/
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
            #region init P0, Vx, Vy

            P0 = new double[Nx, Ny, Nz];
            for (int i = 0; i < Nx; i++)
                for (int j = 0; j < Ny; j++)
                    for(int k = 0; k < Nz; k++)
                        P0[i, j, k] = 0;

            Vx = new double[Nx + 1, Ny, Nz];
            for (int i = 0; i <= Nx; i++)
                for (int j = 0; j < Ny; j++)
                    for (int k = 0; k < Nz; k++)
                        Vx[i, j, k] = 0;

            Vy = new double[Nx, Ny + 1, Nz];
            for(int i = 0; i < Nx; i++)
                for(int j = 0; j <= Ny; j++)
                    for(int k = 0; k < Nz; k++)
                        Vy[i, j, k] = 0;

            Vz = new double[Nx, Ny, Nz + 1];
            for(int i = 0; i < Nx; i++)
                for(int j = 0; j < Ny; j++)
                    for(int k = 0; k <= Nz; k++)
                        P0[i, j, k] = 0;

            #endregion

            #region init Mur 1st

            MurX = new double[4, Ny, Nz];
            for(int j = 0; j < Ny; j++)
                for(int k = 0; j < Nz; j++)
                {
                    MurX[0, j, k] = 0.0;
                    MurX[1, j, k] = 0.0;
                    MurX[2, j, k] = 0.0;
                    MurX[3, j, k] = 0.0;
                }

            MurY = new double[Nx, 4, Nz];
            for(int i = 0; i < Nx; i++)
                for(int k = 0; k < Nz; k++)
                {
                    MurY[i, 0, k] = 0.0;
                    MurY[i, 1, k] = 0.0;
                    MurY[i, 2, k] = 0.0;
                    MurY[i, 3, k] = 0.0;
                }

            MurZ = new double[Nx, Ny, 4];
            for(int i = 0; i < Nx; i++)
                for(int j = 0; j < Ny; j++)
                {
                    MurZ[i, j, 0] = 0.0;
                    MurZ[i, j, 1] = 0.0;
                    MurZ[i, j, 2] = 0.0;
                    MurZ[i, j, 3] = 0.0;
                }

            #endregion

            tCurrent = 0;
        }


        public void Close()
        {
            datWriter.Close();
        }


        public void SaveToCSV()
        {
            datWriter.WriteLine(("" + P0[Nx / 2, Ny / 2, Nz - 4] + ";").Replace(',', '.'));
        }

        public void FlushCSV()
        {
            datWriter.Flush();
        }

        double F(int i, int j, int k, double t)
        {
            if( (Math.Abs(i - Nx / 2) <= 2) &&
                (Math.Abs(j - Ny / 2) <= 2) &&
                (Math.Abs(k - Nz / 2) <= 2) && t < 0.2)
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
            UpdateP();
            MurBoundaries();  
            
            tCurrent += dt;
        }

        #region Boundaries

        private void MurBoundaries()
        {
            double reflectionCoefficient = 0.2;

            #region X повне згасання сигналу

            // x == 0
            for(int j = 1; j < Ny - 1; j++)
                for(int k = 0; k < Nz; k++)

                    P0[0, j, k] =
                        MurX[1, j, k] +
                        (velocity * dt - dx) / (velocity * dt + dx) * (P0[1, j, k] - MurX[0, j, k]);

            // x == Lx
            for(int j = 1; j < Ny - 1; j++)
                for(int k = 0; k < Nz; k++)

                    P0[Nx - 1, j, k] =
                        MurX[2, j, k] +
                        (velocity * dt - dx) / (velocity * dt + dx) * (P0[Nx - 2, j, k] - MurX[3, j, k]);


            #endregion

            #region Y повне згасання сигналу

            // y == 0
            for(int i = 1; i < Nx - 1; i++)
                for(int k = 0; k < Nz; k++)

                    P0[i, 0, k] =
                        MurY[i, 1, k] +
                        (velocity * dt - dy) / (velocity * dt + dy) * (P0[i, 1, k] - MurY[i, 0, k]);

            // y == Ly
            for(int i = 1; i < Nx - 1; i++)
                for(int k = 0; k < Nz; k++)

                    P0[i, Ny - 1, k] =
                        MurY[i, 2, k] +
                        (velocity * dt - dy) / (velocity * dt + dy) * (P0[i, Ny - 2, k] - MurY[i, 3, k]);


            #endregion

            #region Z відбиття від поверхні і дна

            // z == 0
            for(int i = 1; i < Nx - 1; i++)
                for(int j = 1; j < Ny - 1; j++)

                    P0[i, j, 0] =
                        reflectionCoefficient * P0[i, j, 0] +
                        (1 - reflectionCoefficient) * (MurZ[i, j, 1] +
                        (velocity * dt - dz) / (velocity * dt + dz) * (P0[i, j, 1] - MurZ[i, j, 0]));

            // z == Lz
            for(int i = 1; i < Nx - 1; i++)
                for(int j = 1; j < Ny - 1; j++)

                    P0[i, j, Nz - 1] =
                        reflectionCoefficient * P0[i, j, Nz - 1] +
                        (1 - reflectionCoefficient) * (MurZ[i, j, 2] +
                        (velocity * dt - dz) / (velocity * dt + dz) * (P0[i, j, Nz - 2] - MurZ[i, j, 3]));

            #endregion

            MurCopy();
        }

        private void MurCopy()
        {
            for(int j = 0; j < Ny; j++)
                for(int k = 0; k < Nz; k++)
                {
                    MurX[0, j, k] = P0[0, j, k];
                    MurX[1, j, k] = P0[1, j, k];
                    MurX[2, j, k] = P0[Nx - 2, j, k];
                    MurX[3, j, k] = P0[Nx - 1, j, k];
                }

            for(int i = 0; i < Nx; i++)
                for(int k = 0; k < Nz; k++)
                {
                    MurY[i, 0, k] = P0[i, 0, k];
                    MurY[i, 1, k] = P0[i, 1, k];
                    MurY[i, 2, k] = P0[i, Ny - 2, k];
                    MurY[i, 3, k] = P0[i, Ny - 1, k];
                }

            for(int i = 0; i < Nx; i++)
                for(int j = 0; j < Ny; j++)
                {
                    MurZ[i, j, 0] = P0[i, j, 0];
                    MurZ[i, j, 1] = P0[i, j, 1];
                    MurZ[i, j, 2] = P0[i, j, Nz - 2];
                    MurZ[i, j, 3] = P0[i, j, Nz - 1];
                }
        }

        #endregion

        #region Update P, V

        private void UpdateV()
        {
            double dt_dx_ro = dt / (dx * ro);
            double dt_dy_ro = dt / (dy * ro);
            double dt_dz_ro = dt / (dz * ro);

            for (int i = 1; i < Nx; i++)
                for (int j = 0; j < Ny; j++)
                    for (int k = 0; k < Nz; k++)
                        Vx[i, j, k] -= dt_dx_ro * (P0[i, j, k] - P0[i - 1, j , k]);

            for (int i = 0; i < Nx; i++)
                for (int j = 1; j < Ny; j++)
                    for (int k = 0; k < Nz; k++)
                        Vy[i, j, k] -= dt_dy_ro * (P0[i, j, k] - P0[i, j - 1, k]);

            for (int i = 0; i < Nx; i++)
                for (int j = 0; j < Ny; j++)
                    for (int k = 1; k < Nz; k++)
                        Vz[i, j, k] -= dt_dz_ro * (P0[i, j, k] - P0[i , j, k - 1]);

        }

        //

        private void UpdateP()
        {
            double dt_dx_ro = (dt / dx) * ro * velocity * velocity;
            double dt_dy_ro = (dt / dy) * ro * velocity * velocity;
            double dt_dz_ro = (dt / dz) * ro * velocity * velocity;

            for (int i = 0; i < Nx; i++)
                for(int j = 0; j < Ny; j++)
                    for (int k = 0; k < Nz; k++)
                    {
                        P0[i, j, k] -=
                        dt_dx_ro * (Vx[i + 1, j, k] - Vx[i, j, k]) +
                        dt_dy_ro * (Vy[i, j + 1, k] - Vy[i, j, k]) +
                        dt_dz_ro * (Vz[i, j, k + 1] - Vz[i, j, k]) - dt * this.F(i, j, k, tCurrent);
                    }
        }

        #endregion

        #region save

        public void SaveCurrentValuesP(String fileName)
        {
            using StreamWriter outWriter = new StreamWriter(fileName, false, System.Text.Encoding.Default);

            for(int i = 0; i < Nx; i++)
                for(int j = 0; j < Ny; j++)
                    outWriter.WriteLine((dx * i + " " + dy * j + " " + P0[i, j, Nz / 2]).Replace(',', '.'));
        }

        public double GettCurrent()
        {
            return tCurrent;
        }

        #endregion
    }
}
