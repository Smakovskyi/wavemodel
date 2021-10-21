
using System;
using System.IO;

namespace wavemodel
{
    public class PDESolver
    {
        const float Step = 10;

        #region var

        int Step = 50;
        StreamWriter datWriter;

        float[,,] P0;
        float[,,] Vx; int Nx;
        float[,,] Vy; int Ny;
        float[,,] Vz; int Nz;

        float[,,] MurX;
        float[,,] MurY;
        float[,,] MurZ;
        float[] velocityMur;

        float dt;
        float tCurrent;

        float dt_dx_ro;
        float [] dt_dx_vvro;

        #endregion

        #region init

        public void Init(float dt, int Nx, int Ny, int Nz)
        {
            this.dt = dt;

            this.Nx = Nx;
            this.Ny = Ny;
            this.Nz = Nz;

            datWriter = new StreamWriter("P_all.csv", false, System.Text.Encoding.Default);

            InitGrid();
        }
        public void SetCoefficients(float ro, float[] velocity)
        {
            velocityMur = new float[Nz];
            for(int k = 0; k < Nz; k++)
                velocityMur[k] = (velocity[k] * dt - Step) / (velocity[k] * dt + Step);

            dt_dx_ro = dt / (Step * ro);
            dt_dx_vvro = new float[Nz];
            
            for(int k = 0; k < Nz; k++)
                dt_dx_vvro[k] = (dt / Step) * velocity[k] * velocity[k] * ro;
        }

        void InitGrid()
        {
            #region init P0, Vx, Vy

            P0 = new float[Nx, Ny, Nz];
            for (int i = 0; i < Nx; i++)
                for (int j = 0; j < Ny; j++)
                    for(int h = 0; h < Nz; h++)
                        P0[i, j, h] = 0.0f;

            Vx = new float[Nx + 1, Ny, Nz];
            for (int i = 0; i <= Nx; i++)
                for (int j = 0; j < Ny; j++)
                    for (int h = 0; h < Nz; h++)
                        Vx[i, j, h] = 0.0f;

            Vy = new float[Nx, Ny + 1, Nz];
            for(int i = 0; i < Nx; i++)
                for(int j = 0; j <= Ny; j++)
                    for(int h = 0; h < Nz; h++)
                        Vy[i, j, h] = 0.0f;

            Vz = new float[Nx, Ny, Nz + 1];
            for(int i = 0; i < Nx; i++)
                for(int j = 0; j < Ny; j++)
                    for(int h = 0; h <= Nz; h++)
                        Vz[i, j, h] = 0.0f;

            #endregion

            #region init Mur 1st

            MurX = new float[4, Ny, Nz];
            for(int j = 0; j < Ny; j++)
                for(int h = 0; h < Nz; h++)
                {
                    MurX[0, j, h] = 0.0f;
                    MurX[1, j, h] = 0.0f;
                    MurX[2, j, h] = 0.0f;
                    MurX[3, j, h] = 0.0f;
                }

            MurY = new float[Nx, 4, Nz];
            for(int i = 0; i < Nx; i++)
                for(int h = 0; h < Nz; h++)
                {
                    MurY[i, 0, h] = 0.0f;
                    MurY[i, 1, h] = 0.0f;
                    MurY[i, 2, h] = 0.0f;
                    MurY[i, 3, h] = 0.0f;
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

        float F(int i, int j, int  k, float t)
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
            
            for(int j = 0; j < Ny; j++)
                for(int h = 0; h < Nz; h++)
                {
                    P0[0, j, h] = P0[1, j, h];
                    P0[Nx - 1, j, h] = P0[Nx - 1, j, h];
                }

            tCurrent += dt;
        }

        #region Boundaries

        private void MurBoundaries()
        {
            float reflectionCoefficient = 0.2f;

            #region X повне згасання сигналу

            /*
            for(int j = 1; j < Ny - 1; j++)
                for(int h = 0; h < Nz; h++)
                {
                    P0[0, j, h] =
                        MurX[1, j, h] +
                        (P0[1, j, h] - MurX[0, j, h]) * velocityMur[h];

                    P0[0, j, k] =
                        MurX[1, j, k] +
                        (P0[1, j, k] - MurX[0, j, k]) * velocityMur[k];

            // x == Lx
            for(int j = 1; j < Ny - 1; j++)
                for(int k = 0; k < Nz; k++)

                    P0[Nx - 1, j, k] =
                        MurX[2, j, k] +
                        (P0[Nx - 2, j, k] - MurX[3, j, k]) * velocityMur[k];


            #endregion

            #region Y повне згасання сигналу

            // y == 0
            for(int i = 1; i < Nx - 1; i++)
                for(int h = 0; h < Nz; h++)

                    P0[i, 0, k] =
                        MurY[i, 1, k] +
                        (P0[i, 1, k] - MurY[i, 0, k]) * velocityMur[k];

            // y == Ly
            for(int i = 1; i < Nx - 1; i++)
                for(int h = 0; h < Nz; h++)

                    P0[i, Ny - 1, k] =
                        MurY[i, 2, k] +
                        (P0[i, Ny - 2, k] - MurY[i, 3, k]) * velocityMur[k];


            #endregion

            #region Z відбиття від поверхні і дна

            // z == 0
            for(int i = 1; i < Nx - 1; i++)
                for(int j = 1; j < Ny - 1; j++)

                    P0[i, j, 0] =
                        reflectionCoefficient * P0[i, j, 0] + (1 - reflectionCoefficient) *
                            (MurZ[i, j, 1] +
                            (P0[i, j, 1] - MurZ[i, j, 0]) * velocityMur[0]);

            // z == Lz
            for(int i = 1; i < Nx - 1; i++)
                for(int j = 1; j < Ny - 1; j++)

                    P0[i, j, Nz - 1] =
                        reflectionCoefficient * P0[i, j, Nz - 1] + (1 - reflectionCoefficient) *
                            (MurZ[i, j, 2] +
                            (P0[i, j, Nz - 2] - MurZ[i, j, 3]) * velocityMur[Nz - 1]);

            #endregion

            MurCopy();
        }

        private void MurCopy()
        {
            for(int j = 0; j < Ny; j++)
                for(int h = 0; h < Nz; h++)
                {
                    MurX[0, j, h] = P0[0, j, h];
                    MurX[1, j, h] = P0[1, j, h];
                    MurX[2, j, h] = P0[Nx - 2, j, h];
                    MurX[3, j, h] = P0[Nx - 1, j, h];
                }

            for(int i = 0; i < Nx; i++)
                for(int h = 0; h < Nz; h++)
                {
                    MurY[i, 0, h] = P0[i, 0, h];
                    MurY[i, 1, h] = P0[i, 1, h];
                    MurY[i, 2, h] = P0[i, Ny - 2, h];
                    MurY[i, 3, h] = P0[i, Ny - 1, h];
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
                    for (int h = 0; h < Nz; h++)
                        Vx[i, j, h] -= dt_dx_ro * (P0[i, j, h] - P0[i - 1, j, h]);

            for (int i = 0; i < Nx; i++)
                for (int j = 1; j < Ny; j++)
                    for (int h = 0; h < Nz; h++)
                        Vy[i, j, h] -= dt_dx_ro * (P0[i, j, h] - P0[i, j - 1, h]);

            for (int i = 0; i < Nx; i++)
                for (int j = 0; j < Ny; j++)
                    for (int h = 1; h < Nz; h++)
                        Vz[i, j, h] -= dt_dx_ro * (P0[i, j, h] - P0[i , j, h - 1]);
        }

        //

        private void UpdateP()
        {
            for (int i = 0; i < Nx; i++)
                for(int j = 0; j < Ny; j++)
                    for (int h = 0; h < Nz; h++)
                    {
                        P0[i, j, k] -=
                            dt_dx_vvro[k] *
                                (
                                Vx[i + 1, j, h] - Vx[i, j, h] +
                                Vy[i, j + 1, h] - Vy[i, j, h] +
                                Vz[i, j, h + 1] - Vz[i, j, h]
                                ) - dt * this.F(i, j, h, tCurrent);
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
                    outWriter.WriteLine((Step * i + " " + Step * j + " " + P0[i, j, Nz / 2]).Replace(',', '.'));
        }
        public float GettCurrent()
        {
            return tCurrent;
        }

        #endregion
    }
}
