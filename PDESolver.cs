
using System;
using System.IO;

namespace wavemodel
{
    public class PDESolver
    {
        #region var

        StreamWriter datWriter;

        float[,,] P0;
        float[,,] Vx; int dx; int Nx;
        float[,,] Vy; int dy; int Ny;
        float[,,] Vz; int dz; int Nz;

        float[,,] MurX;
        float[,,] MurY;
        float[,,] MurZ;
        float [,,] velocityMur;

        float tCurrent;
        float dt;

        float dt_dx_ro;
        float dt_dy_ro;
        float dt_dz_ro;
        float [,,] dt_dy_vvro;
        float [,,] dt_dz_vvro;
        float [,,] dt_dx_vvro;

        #endregion

        #region unused

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
        public double Psi(double x, double y)
        {
            return 0.0;
        }
        */

        #endregion

        #region init

        public void Init(float dt, int lX, int lY, int lZ, int Nx, int Ny, int Nz)
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
        public void SetCoefficients(int ro, float[,,] velocity)
        {
            velocityMur = new float[Nx, Ny, Nz];
            for(int i = 0; i < Nx; i++)
                for(int j = 0; j < Ny; j++)
                    for(int k = 0; k < Nz; k++)
                        velocityMur[i, j, k] = (velocity[i, j, k] * dt - dz) / (velocity[i, j, k] * dt + dz);

            dt_dx_ro = dt / (dx * ro);
            dt_dy_ro = dt / (dy * ro);
            dt_dz_ro = dt / (dz * ro);

            dt_dx_vvro = new float[Nx, Ny, Nz];
            dt_dy_vvro = new float[Nx, Ny, Nz];
            dt_dz_vvro = new float[Nx, Ny, Nz];
            for(int i = 0; i < Nz; i++)
                for(int j = 0; j < Nz; j++)
                    for(int k = 0; k < Nz; k++)
                    {
                        dt_dx_vvro[i, j, k] = (dt / dx) * velocity[i, j, k] * velocity[i, j, k] * ro;
                        dt_dy_vvro[i, j, k] = (dt / dy) * velocity[i, j, k] * velocity[i, j, k] * ro;
                        dt_dz_vvro[i, j, k] = (dt / dz) * velocity[i, j, k] * velocity[i, j, k] * ro;
                    }
        }

        void InitGrid()
        {
            #region init P0, Vx, Vy

            P0 = new float[Nx, Ny, Nz];
            for (int i = 0; i < Nx; i++)
                for (int j = 0; j < Ny; j++)
                    for(int k = 0; k < Nz; k++)
                        P0[i, j, k] = 0.0f;

            Vx = new float[Nx + 1, Ny, Nz];
            for (int i = 0; i <= Nx; i++)
                for (int j = 0; j < Ny; j++)
                    for (int k = 0; k < Nz; k++)
                        Vx[i, j, k] = 0.0f;

            Vy = new float[Nx, Ny + 1, Nz];
            for(int i = 0; i < Nx; i++)
                for(int j = 0; j <= Ny; j++)
                    for(int k = 0; k < Nz; k++)
                        Vy[i, j, k] = 0.0f;

            Vz = new float[Nx, Ny, Nz + 1];
            for(int i = 0; i < Nx; i++)
                for(int j = 0; j < Ny; j++)
                    for(int k = 0; k <= Nz; k++)
                        Vz[i, j, k] = 0.0f;

            #endregion

            #region init Mur 1st

            MurX = new float[4, Ny, Nz];
            for(int j = 0; j < Ny; j++)
                for(int k = 0; j < Nz; j++)
                {
                    MurX[0, j, k] = 0.0f;
                    MurX[1, j, k] = 0.0f;
                    MurX[2, j, k] = 0.0f;
                    MurX[3, j, k] = 0.0f;
                }

            MurY = new float[Nx, 4, Nz];
            for(int i = 0; i < Nx; i++)
                for(int k = 0; k < Nz; k++)
                {
                    MurY[i, 0, k] = 0.0f;
                    MurY[i, 1, k] = 0.0f;
                    MurY[i, 2, k] = 0.0f;
                    MurY[i, 3, k] = 0.0f;
                }

            MurZ = new float[Nx, Ny, 4];
            for(int i = 0; i < Nx; i++)
                for(int j = 0; j < Ny; j++)
                {
                    MurZ[i, j, 0] = 0.0f;
                    MurZ[i, j, 1] = 0.0f;
                    MurZ[i, j, 2] = 0.0f;
                    MurZ[i, j, 3] = 0.0f;
                }

            #endregion

            tCurrent = 0.0f;
        }

        //

        float F(int i, int j, int  k, double t)
        {
            if( (Math.Abs(i - Nx / 2) <= 2) &&
                (Math.Abs(j - Ny / 2) <= 2) &&
                (Math.Abs(k - Nz / 2) <= 2) && t < 0.2)
            {
                return (float)Math.Cos(Math.PI * 2f * t);
            }
            return 0.0f;
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
            float reflectionCoefficient = 0.2f;

            #region X повне згасання сигналу

            // x == 0
            for(int j = 1; j < Ny - 1; j++)
                for(int k = 0; k < Nz; k++)

                    P0[0, j, k] =
                        MurX[1, j, k] +
                        (P0[1, j, k] - MurX[0, j, k]) * velocityMur[0, j, k];

            // x == Lx
            for(int j = 1; j < Ny - 1; j++)
                for(int k = 0; k < Nz; k++)

                    P0[Nx - 1, j, k] =
                        MurX[2, j, k] +
                        (P0[Nx - 2, j, k] - MurX[3, j, k]) * velocityMur[Nx - 1, j, k];


            #endregion

            #region Y повне згасання сигналу

            // y == 0
            for(int i = 1; i < Nx - 1; i++)
                for(int k = 0; k < Nz; k++)

                    P0[i, 0, k] =
                        MurY[i, 1, k] +
                        (P0[i, 1, k] - MurY[i, 0, k]) * velocityMur[i, 0, k];

            // y == Ly
            for(int i = 1; i < Nx - 1; i++)
                for(int k = 0; k < Nz; k++)

                    P0[i, Ny - 1, k] =
                        MurY[i, 2, k] +
                        (P0[i, Ny - 2, k] - MurY[i, 3, k]) * velocityMur[i, Ny - 1, k];


            #endregion

            #region Z відбиття від поверхні і дна

            // z == 0
            for(int i = 1; i < Nx - 1; i++)
                for(int j = 1; j < Ny - 1; j++)

                    P0[i, j, 0] =
                        reflectionCoefficient * P0[i, j, 0] + (1 - reflectionCoefficient) *
                            (MurZ[i, j, 1] +
                            (P0[i, j, 1] - MurZ[i, j, 0]) * velocityMur[i, j, 0]);

            // z == Lz
            for(int i = 1; i < Nx - 1; i++)
                for(int j = 1; j < Ny - 1; j++)

                    P0[i, j, Nz - 1] =
                        reflectionCoefficient * P0[i, j, Nz - 1] + (1 - reflectionCoefficient) *
                            (MurZ[i, j, 2] +
                            (P0[i, j, Nz - 2] - MurZ[i, j, 3]) * velocityMur[i, j, Nz - 1]);

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
            for (int i = 0; i < Nx; i++)
                for(int j = 0; j < Ny; j++)
                    for (int k = 0; k < Nz; k++)
                    {
                        P0[i, j, k] -=
                        dt_dx_vvro[i, j, k] * (Vx[i + 1, j, k] - Vx[i, j, k]) +
                        dt_dy_vvro[i, j, k] * (Vy[i, j + 1, k] - Vy[i, j, k]) +
                        dt_dz_vvro[i, j, k] * (Vz[i, j, k + 1] - Vz[i, j, k]) - dt * this.F(i, j, k, tCurrent);
                    }
        }

        #endregion

        #region save

        public void SaveToCSV()
        {
            datWriter.WriteLine(("" + P0[Nx / 2, Ny / 2, Nz - 4] + ";").Replace(',', '.'));
        }
        public void FlushCSV()
        {
            datWriter.Flush();
        }
        public void Close()
        {
            datWriter.Close();
        }

        //

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
