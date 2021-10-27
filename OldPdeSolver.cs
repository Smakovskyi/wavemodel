using System;
using System.IO;

namespace wavemodel
{
    public class OldPDESolver
    {
        #region var

        int Step = 100;
        StreamWriter datWriter;

        float[,,] P0;
        float[,,] Vx; int Nx;
        float[,,] Vy; int Ny;
        float[,,] Vz; int Nz;

        float[,,] MurX;
        float[,,] MurY;
        float[,,] MurZ;
        float[] velocityMur;
        float[,] bathymetry;
        int[,] maxK;
        //angles for normal vector to the bottom
        float[,] cosA2;
        float[,] cosB2;
        float[,] cosY2;

        float dt;
        float tCurrent;

        float dt_dx_ro;
        float[] dt_dx_vvro;
        float reflectionCoefficient = 0.2f;

        #endregion

        #region init

        public void Init(float dt, int Nx, int Ny, int Nz, float[,] bathymetry)
        {
            this.dt = dt;

            this.Nx = Nx;
            this.Ny = Ny;
            this.Nz = Nz;

            this.bathymetry = bathymetry;

            datWriter = new StreamWriter("P_all.csv", false, System.Text.Encoding.Default);

            InitMaxK();
            InitGrid();
            initAngles();
            MurInit();

        }



        public void SetCoefficients(float ro, float[] velocity)
        {
            velocityMur = new float[Nz];
            for (int k = 0; k < Nz; k++)
                velocityMur[k] = (velocity[k] * dt - Step) / (velocity[k] * dt + Step);

            dt_dx_ro = dt / (Step * ro);
            dt_dx_vvro = new float[Nz];

            for (int k = 0; k < Nz; k++)
                dt_dx_vvro[k] = (dt / Step) * velocity[k] * velocity[k] * ro;
        }

        void InitGrid()
        {
            #region init P0, Vx, Vy

            P0 = new float[Nx, Ny, Nz];
            for (int i = 0; i < Nx; i++)
                for (int j = 0; j < Ny; j++)
                    for (int k = 0; k < Nz; k++)
                        P0[i, j, k] = 0.0f;

            Vx = new float[Nx + 1, Ny, Nz];
            for (int i = 0; i <= Nx; i++)
                for (int j = 0; j < Ny; j++)
                    for (int k = 0; k < Nz; k++)
                        Vx[i, j, k] = 0.0f;

            Vy = new float[Nx, Ny + 1, Nz];
            for (int i = 0; i < Nx; i++)
                for (int j = 0; j <= Ny; j++)
                    for (int k = 0; k < Nz; k++)
                        Vy[i, j, k] = 0.0f;

            Vz = new float[Nx, Ny, Nz + 1];
            for (int i = 0; i < Nx; i++)
                for (int j = 0; j < Ny; j++)
                    for (int k = 0; k <= Nz; k++)
                        Vz[i, j, k] = 0.0f;

            #endregion



            tCurrent = 0.0f;
        }

        private void MurInit()
        {
            #region init Mur 1st

            MurX = new float[4, Ny, Nz];
            for (int j = 0; j < Ny; j++)
                for (int k = 0; k < Nz; k++)
                {
                    MurX[0, j, k] = 0.0f;
                    MurX[1, j, k] = 0.0f;
                    MurX[2, j, k] = 0.0f;
                    MurX[3, j, k] = 0.0f;
                }

            MurY = new float[Nx, 4, Nz];
            for (int i = 0; i < Nx; i++)
                for (int k = 0; k < Nz; k++)
                {
                    MurY[i, 0, k] = 0.0f;
                    MurY[i, 1, k] = 0.0f;
                    MurY[i, 2, k] = 0.0f;
                    MurY[i, 3, k] = 0.0f;
                }

            MurZ = new float[Nx, Ny, 4];
            for (int i = 0; i < Nx; i++)
                for (int j = 0; j < Ny; j++)
                {
                    MurZ[i, j, 0] = 0.0f;
                    MurZ[i, j, 1] = 0.0f;
                    MurZ[i, j, 2] = 0.0f;
                    MurZ[i, j, 3] = 0.0f;
                }

            #endregion

        }

        void initAngles()
        {
            cosA2 = new float[Nx + 1, Ny + 1];
            cosB2 = new float[Nx + 1, Ny + 1];
            cosY2 = new float[Nx + 1, Ny + 1];
            initAnglesForBottom();
        }

        // HACK: interpolation mass. re-do to bath
        void initAnglesForBottom()
        {
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    float dfx = (bathymetry[i + 1, j] - bathymetry[i, j]) / Step;
                    float dfy = (bathymetry[i, j + 1] - bathymetry[i, j]) / Step;
                    float len = (float)Math.Sqrt(dfx * dfx + dfy * dfy + 1);
                    cosA2[i, j] = dfx / len;
                    cosB2[i, j] = dfy / len;
                    cosY2[i, j] = -1 / len;
                }
                cosA2[i, Ny] = cosA2[i, Ny - 1];
                cosB2[i, Ny] = cosB2[i, Ny - 1];
                cosY2[i, Ny] = cosY2[i, Ny - 1];
            }
            for (int j = 0; j <= Ny; j++)
            {
                cosA2[Nx, j] = cosA2[Nx - 1, j];
                cosB2[Nx, j] = cosB2[Nx - 1, j];
                cosY2[Nx, j] = cosY2[Nx - 1, j];
            }
        }
        private void InitMaxK()
        {
            maxK = new int[Nx + 1, Ny + 1];

            for (int i = 0; i <= Nx; i++)
            {
                for (int j = 0; j <= Ny; j++)
                {
                    double Deepth = bathymetry[i, j];
                    maxK[i, j] = (int)(Deepth / Step);
                }
            }
        }

        //

        float F(int i, int j, int k, float t)
        {
            if ((Math.Abs(i - Nx / 2) <= 2) &&
                (Math.Abs(j - Ny / 2) <= 2) &&
                (Math.Abs(k - Nz / 2) <= 2) && t < 0.2)
            {
                return (float)Math.Cos(Math.PI * 2 * t);
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


            #region X повне згасання сигналу

            // x == 0
            int ii = 0;
            for (int j = 1; j < Ny - 1; j++)
                for (int k = 0; k < maxK[ii, j]; k++)

                    P0[0, j, k] =
                        MurX[1, j, k] +
                        (P0[1, j, k] - MurX[0, j, k]) * velocityMur[k];

            // x == Lx
            ii = Nx;
            for (int j = 1; j < Ny - 1; j++)
                for (int k = 0; k < maxK[ii, j]; k++)

                    P0[Nx - 1, j, k] =
                        MurX[2, j, k] +
                        (P0[Nx - 2, j, k] - MurX[3, j, k]) * velocityMur[k];

            #endregion

            #region Y повне згасання сигналу

            // y == 0
            int jj = 0;
            for (int i = 1; i < Nx - 1; i++)
                for (int k = 0; k < maxK[i, jj]; k++)

                    P0[i, 0, k] =
                        MurY[i, 1, k] +
                        (P0[i, 1, k] - MurY[i, 0, k]) * velocityMur[k];

            // y == Ly
            jj = Ny;
            for (int i = 1; i < Nx - 1; i++)
                for (int k = 0; k < maxK[i, jj]; k++)

                    P0[i, Ny - 1, k] =
                        MurY[i, 2, k] +
                        (P0[i, Ny - 2, k] - MurY[i, 3, k]) * velocityMur[k];


            #endregion

            #region Z відбиття від поверхні і дна
            TopReflection();
            BottomReflection();



            #endregion

            MurCopy();
        }

        private void TopReflection()
        {
            // z == 0
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int j = 1; j < Ny - 1; j++)
                {
                    int kMax = maxK[i, j];
                    if (kMax > 2)
                    {
                        P0[i, j, 0] = reflectionCoefficient * P0[i, j, 0] +
                                (1 - reflectionCoefficient) *
                                (MurZ[i, j, 1] + (P0[i, j, 1] - MurZ[i, j, 0]) * velocityMur[0]);
                    }
                }
            }

        }

        private void BottomReflection()
        {

            for (int i = 1; i < Nx - 1; i++)
                for (int j = 1; j < Ny - 1; j++)
                {
                    int kMax = maxK[i, j];
                    if (kMax > 2)
                    {
                        P0[i, j, kMax - 1] =
                            /*reflectionCoefficient * P0[i, j, kMax - 1] + (1 - reflectionCoefficient) **/
                            //P0[i, j, Nz - 2]                     P0[i, j, Nz - 1]
                            (MurZ[i, j, 2] + (P0[i, j, kMax - 2] - MurZ[i, j, 3]) * velocityMur[kMax - 1]);
                    }
                }
        }

        private void MurCopy()
        {
            for (int j = 0; j < Ny; j++)
                for (int k = 0; k < Nz; k++)
                {
                    MurX[0, j, k] = P0[0, j, k];
                    MurX[1, j, k] = P0[1, j, k];
                    MurX[2, j, k] = P0[Nx - 2, j, k];
                    MurX[3, j, k] = P0[Nx - 1, j, k];
                }

            for (int i = 0; i < Nx; i++)
                for (int k = 0; k < Nz; k++)
                {
                    MurY[i, 0, k] = P0[i, 0, k];
                    MurY[i, 1, k] = P0[i, 1, k];
                    MurY[i, 2, k] = P0[i, Ny - 2, k];
                    MurY[i, 3, k] = P0[i, Ny - 1, k];
                }

            for (int i = 0; i < Nx; i++)
                for (int j = 0; j < Ny; j++)
                {
                    MurZ[i, j, 0] = P0[i, j, 0];
                    MurZ[i, j, 1] = P0[i, j, 1];
                    int kMax = maxK[i, j];
                    if (kMax > 2)
                    {
                        MurZ[i, j, 2] = P0[i, j, kMax - 2];
                        MurZ[i, j, 3] = P0[i, j, kMax - 1];
                    }
                }
        }

        #endregion

        #region Update P, V

        private void UpdateV()
        {
            for (int i = 1; i < Nx; i++)
                for (int j = 0; j < Ny; j++)
                {
                    int MaxK = maxK[i, j];
                    for (int k = 0; k < MaxK; k++)
                        Vx[i, j, k] -= dt_dx_ro * (P0[i, j, k] - P0[i - 1, j, k]);
                }


            for (int i = 0; i < Nx; i++)
                for (int j = 1; j < Ny; j++)
                {
                    int MaxK = maxK[i, j];
                    for (int k = 0; k < MaxK; k++)
                        Vy[i, j, k] -= dt_dx_ro * (P0[i, j, k] - P0[i, j - 1, k]);
                }


            for (int i = 0; i < Nx; i++)
                for (int j = 0; j < Ny; j++)
                {
                    int MaxK = maxK[i, j];
                    for (int k = 1; k < MaxK; k++)
                        Vz[i, j, k] -= dt_dx_ro * (P0[i, j, k] - P0[i, j, k - 1]);
                }
        }

        //

        private void UpdateP()
        {
            for (int i = 0; i < Nx; i++)
                for (int j = 0; j < Ny; j++)
                {

                    int MaxK = maxK[i, j];
                    for (int k = 0; k < MaxK; k++)
                    {
                        P0[i, j, k] -=
                            dt_dx_vvro[k] *
                                (
                                Vx[i + 1, j, k] - Vx[i, j, k] +
                                Vy[i, j + 1, k] - Vy[i, j, k] +
                                Vz[i, j, k + 1] - Vz[i, j, k]
                                ) - dt * this.F(i, j, k, tCurrent);
                    }
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

            for (int i = 0; i < Nx; i++)
                for (int j = 0; j < Ny; j++)
                    outWriter.WriteLine((Step * i + " " + Step * j + " " + P0[i, j, Nz / 2]).Replace(',', '.'));
        }
        public float GettCurrent()
        {
            return tCurrent;
        }

        #endregion
    }
}

