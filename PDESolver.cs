
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

        /* https://github.com/nagataniyoshiki/FDTD_processing/blob/master/acoustic_2d_fdtd_processing2/acoustic_2d_fdtd_processing2.pde */
        double[,] Mur_X1;   // Mur's 2nd-Order Absorption Layer
        double[,] Mur_X2;   // Mur's 2nd-Order Absorption Layer
        double[,] Mur_Y1;   // Mur's 2nd-Order Absorption Layer
        double[,] Mur_Y2;   // Mur's 2nd-Order Absorption Layer

        double tCurrent;
        double dt;
        double dx;
        double dy;
        double ro;
        double velocity;

        int Nx;
        int Ny;

        double dt_dx_ro = 0;
        double dt_dy_ro = 0;

        double dt_dx_vro = 0;
        double dt_dy_vro = 0;

        #endregion

        #region unused 1

        //x = 0 border
        static double MuX(double y, double t)
        {
            return 0;
        }
        //x = x_max border
        public double NuX(double y, double t)
        {
            return 0;
        }
        //y = 0
        public double MuY(double x, double t)
        {
            return 0;
        }
        //y = y_max
        public double NuY(double x, double t)
        {
            return 0;
        }
        public double Fi(double x, double y)
        {
            return 0;
        }
        public double Psi(double x, double y)
        {
            return 0;
        }

        #endregion

        #region init

        public void Init(double dt, double lX, double lY, int Nx, int Ny)
        {
            this.Nx = (int)((double)Nx / lX);
            this.Ny = (int)((double)Ny / lY);
            this.dt = dt;

            this.dx = lX;
            this.dy = lY;

            InitGrid();
        }
        void InitGrid()
        {
            #region init P0, Vx, Vy

            P0 = new double[Nx + 1, Ny + 1];
            Vx = new double[Nx + 1, Ny + 1];
            Vy = new double[Nx + 1, Ny + 1];

            for (int i = 0; i <= Nx; i++)
                for (int j = 0; j < Ny; j++)
                    Vx[i, j] = 0;

            for(int i = 0; i < Nx; i++)
                for(int j = 0; j <= Ny; j++)
                    Vy[i, j] = 0;

            for(int i = 0; i < Nx; i++)
                for(int j = 0; j < Ny; j++)
                    P0[i, j] = 0;


            dt_dx_ro = (dt / dx) * ro;
            dt_dy_ro = (dt / dy) * ro;

            dt_dx_vro = dt_dx_ro * velocity * velocity;
            dt_dy_vro = dt_dy_ro * velocity * velocity;

            #endregion

            #region init Mur

            Mur_X1 = new double[4, Ny + 1];
            Mur_X2 = new double[4, Ny + 1];
            Mur_Y1 = new double[Nx + 1, 4];
            Mur_Y2 = new double[Nx + 1, 4];

            for(int i = 0; i <= Nx; i++)
            {
                Mur_Y1[i, 0] = 0;
                Mur_Y1[i, 1] = 0;
                Mur_Y1[i, 2] = 0;
                Mur_Y1[i, 3] = 0;
                Mur_Y2[i, 0] = 0;
                Mur_Y2[i, 1] = 0;
                Mur_Y2[i, 2] = 0;
                Mur_Y2[i, 3] = 0;
            }

            for(int j = 0; j <= Ny; j++)
            {
                Mur_X1[0, j] = 0;
                Mur_X1[1, j] = 0;
                Mur_X1[2, j] = 0;
                Mur_X1[3, j] = 0;
                Mur_X2[0, j] = 0;
                Mur_X2[1, j] = 0;
                Mur_X2[2, j] = 0;
                Mur_X2[3, j] = 0;
            }

            #endregion

            tCurrent = 0;
        }

        double F(int i, int j, double t)
        {
            if( (Math.Abs( i - Nx / 2) <= 2) && (Math.Abs(j - Ny / 2) <= 2) && t < 0.2)
                return (double)Math.Cos(Math.PI * 2 * t);

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

            FillBoundaries_null();
            //FillBoundaries_1ndOrder();
            //FillBoundaries_2ndOrder();

            tCurrent += dt;
        }

        #region Boundaries

        private void FillBoundaries_null()
        {
            for(int j = 0; j <= Nx; j++)
            {
                P0[0, j] = 0;
                P0[Nx, j] = 0;
                Vx[0, j] = 0;
                Vx[Nx, j] = 0;
                Vy[0, j] = 0;
                Vy[Nx, j] = 0;
            }

            for(int i = 0; i <= Ny; i++)
            {
                P0[i, 0] = 0;
                P0[i, Ny] = 0;
                Vx[i, 0] = 0;
                Vx[i, Ny] = 0;
                Vy[i, 0] = 0;
                Vy[i, Ny] = 0;
            }
        }

        private void FillBoundaries_2ndOrder()
        {
            int i = 0, j = 0;

            #region Mur's 2nd Order Absorption

            for(i = 2; i < Nx - 2; i++)
            {
                P0[i, 0] =
                            -Mur_Y2[i, 1]
                            + (velocity * dt - dx) / (velocity * dt + dx) * (P0[i, 1] + Mur_Y2[i, 0])
                            + (2 * dx) / (velocity * dt + dx) * (Mur_Y1[i, 0] + Mur_Y1[i, 1])
                            + (dx * velocity * velocity * dt * dt) / (2 * dx * dx * (velocity * dt + dx))
                                * ( + Mur_Y1[i + 1, 0] - 2 * Mur_Y1[i, 0]
                                    + Mur_Y1[i - 1, 0] + Mur_Y1[i + 1, 1]
                                    - 2 * Mur_Y1[i, 1] + Mur_Y1[i - 1, 1]);

                P0[i, Ny - 1] =
                            -Mur_Y2[i, 2]
                            + (velocity * dt - dx) / (velocity * dt + dx) * (P0[i, Ny - 2] + Mur_Y2[i, 3])
                            + (2 * dx) / (velocity * dt + dx) * (Mur_Y1[i, 3] + Mur_Y1[i, 2])
                            + (dx * velocity * velocity * dt * dt) / (2 * dx * dx * (velocity * dt + dx))
                                * ( + Mur_Y1[i + 1, 3] - 2 * Mur_Y1[i, 3]
                                    + Mur_Y1[i - 1, 3] + Mur_Y1[i + 1, 2]
                                    - 2 * Mur_Y1[i, 2] + Mur_Y1[i - 1, 2]);
            }

            for(j = 2; j < Ny - 2; j++)
            {
                P0[0, j] =
                            -Mur_X2[1, j]
                            + (velocity * dt - dx) / (velocity * dt + dx) * (P0[1, j] + Mur_X2[0, j])
                            + (2 * dx) / (velocity * dt + dx) * (Mur_X1[0, j] + Mur_X1[1, j])
                            + (dx * velocity * velocity * dt * dt) / (2 * dx * dx * (velocity * dt + dx))
                                * ( + Mur_X1[0, j + 1] - 2 * Mur_X1[0, j]
                                    + Mur_X1[0, j - 1] + Mur_X1[1, j + 1]
                                    - 2 * Mur_X1[1, j] + Mur_X1[1, j - 1]);

                P0[Nx - 1, j] =
                            -Mur_X2[2, j]
                            + (velocity * dt - dx) / (velocity * dt + dx) * (P0[Nx - 2, j] + Mur_X2[3, j])
                            + (2 * dx) / (velocity * dt + dx) * (Mur_X1[3, j] + Mur_X1[2, j])
                            + (dx * velocity * velocity * dt * dt) / (2 * dx * dx * (velocity * dt + dx))
                                * ( + Mur_X1[3, j + 1] - 2 * Mur_X1[3, j]
                                    + Mur_X1[3, j - 1] + Mur_X1[2, j + 1]
                                    - 2 * Mur_X1[2, j] + Mur_X1[2, j - 1]);
            }

            #endregion

            #region Mur's 1st Order Absorption for 4 corners

            i = 1;

            P0[i, 0] = Mur_Y1[i, 1] + (velocity * dt - dx) / (velocity * dt + dx) * (P0[i, 1] - Mur_Y1[i, 0]);
            P0[i, Ny - 1] = Mur_Y1[i, 2] + (velocity * dt - dx) / (velocity * dt + dx) * (P0[i, Ny - 2] - Mur_Y1[i, 3]);


            i = Nx - 2;

            P0[i, 0] = Mur_Y1[i, 1] + (velocity * dt - dx) / (velocity * dt + dx) * (P0[i, 1] - Mur_Y1[i, 0]);
            P0[i, Ny - 1] = Mur_Y1[i, 2] + (velocity * dt - dx) / (velocity * dt + dx) * (P0[i, Ny - 2] - Mur_Y1[i, 3]);


            j = 1;

            P0[0, j] = Mur_X1[1, j] + (velocity * dt - dx) / (velocity * dt + dx) * (P0[1, j] - Mur_X1[0, j]);
            P0[Nx - 1, j] = Mur_X1[2, j] + (velocity * dt - dx) / (velocity * dt + dx) * (P0[Nx - 2, j] - Mur_X1[3, j]);


            j = Ny - 2;

            P0[0, j] = Mur_X1[1, j] + (velocity * dt - dx) / (velocity * dt + dx) * (P0[1, j] - Mur_X1[0, j]);
            P0[Nx - 1, j] = Mur_X1[2, j] + (velocity * dt - dx) / (velocity * dt + dx) * (P0[Nx - 2, j] - Mur_X1[3, j]);

            #endregion

            #region Copy Previous Values

            for(i = 0; i < Nx; i++)
            {
                // Copy 1st Old Values to 2nd Old Values
                Mur_Y2[i, 0] = Mur_Y1[i, 0];
                Mur_Y2[i, 1] = Mur_Y1[i, 1];
                Mur_Y2[i, 2] = Mur_Y1[i, 2];
                Mur_Y2[i, 3] = Mur_Y1[i, 3];

                // Copy Present Values
                Mur_Y1[i, 0] = P0[i, 0];
                Mur_Y1[i, 1] = P0[i, 1];
                Mur_Y1[i, 2] = P0[i, Ny - 2];
                Mur_Y1[i, 3] = P0[i, Ny - 1];
            }

            for(j = 0; j < Ny; j++)
            {
                // Copy 1st Old Values to 2nd Old Values
                Mur_X2[0, j] = Mur_X1[0, j];
                Mur_X2[1, j] = Mur_X1[1, j];
                Mur_X2[2, j] = Mur_X1[2, j];
                Mur_X2[3, j] = Mur_X1[3, j];

                // Copy Present Values
                Mur_X1[0, j] = P0[0, j];
                Mur_X1[1, j] = P0[1, j];
                Mur_X1[2, j] = P0[Nx - 2, j];
                Mur_X1[3, j] = P0[Nx - 1, j];
            }

            #endregion


            #region source

            double freq = 1;

            double Cos = Math.Cos(2 * Math.PI * freq * tCurrent);
            double Sin = Math.Sin(2 * Math.PI * freq * tCurrent);
            double Amp = (1 - Cos) / 2 * Sin;

            if(tCurrent < ((double)1 / freq) / dt)
                P0[Nx / 2, Ny / 2] += Amp;

            #endregion
        }

        private void FillBoundaries_1ndOrder()
        {
            #region Mur's 1st Order Absorption

            for(int i = 1; i < Nx - 1; i++)
            {
                P0[i, 0] = Mur_Y1[i, 1] + (velocity * dt - dx) / (velocity * dt + dx) * (P0[i, 1] - Mur_Y1[i, 0]);
                P0[i, Ny - 1] = Mur_Y1[i, 2] + (velocity * dt - dx) / (velocity * dt + dx) * (P0[i, Ny - 2] - Mur_Y1[i, 3]);
            }

            for(int j = 1; j < Ny - 1; j++)
            {
                P0[0, j] = Mur_X1[1, j] + (velocity * dt - dx) / (velocity * dt + dx) * (P0[1, j] - Mur_X1[0, j]);
                P0[Nx - 1, j] = Mur_X1[2, j] + (velocity * dt - dx) / (velocity * dt + dx) * (P0[Nx - 2, j] - Mur_X1[3, j]);
            }

            #endregion

            #region Copy Previous Values

            for(int i = 0; i < Nx; i++)
            {
                Mur_Y1[i, 0] = P0[i, 0];
                Mur_Y1[i, 1] = P0[i, 1];
                Mur_Y1[i, 2] = P0[i, Ny - 2];
                Mur_Y1[i, 3] = P0[i, Ny - 1];
            }
            for(int j = 0; j < Ny; j++)
            {
                Mur_X1[0, j] = P0[0, j];
                Mur_X1[1, j] = P0[1, j];
                Mur_X1[2, j] = P0[Nx - 2, j];
                Mur_X1[3, j] = P0[Nx - 1, j];
            }

            #endregion


            #region source

            double freq = 1;

            double Cos = Math.Cos(2 * Math.PI * freq * tCurrent);
            double Sin = Math.Sin(2 * Math.PI * freq * tCurrent);
            double Amp = (1 - Cos) / 2 * Sin;

            if(tCurrent < ((double)1 / freq) / dt)
                P0[Nx / 2, Ny / 2] += Amp;

            #endregion
        }

        #endregion
        //

        private void UpdateV()
        {
            for(int i = 1; i < Nx; i++)
                for(int j = 0; j < Ny; j++)
                    Vx[i, j] += -dt_dx_ro * (P0[i, j] - P0[i - 1, j]);

            for(int i = 0; i < Nx; i++)
                for(int j = 1; j < Ny; j++)
                    Vy[i, j] += -dt_dy_ro * (P0[i, j] - P0[i, j - 1]);
        }
        private void UpdateP()
        {
            for(int i = 0; i < Nx; i++)
                for(int j = 0; j < Ny; j++)
                    P0[i, j] +=
                        -dt_dx_vro * (Vx[i + 1, j] - Vx[i, j]) +
                        -dt_dy_vro * (Vy[i, j + 1] - Vy[i, j]) - dt * this.F(i, j, tCurrent);
        }


        #region save

        public void SaveCurrentValuesP0(String fileName)
        {
            using StreamWriter outWriter = new StreamWriter(fileName, false, System.Text.Encoding.Default);

            for(int i = 0; i < Nx; i++)
                for(int j = 0; j < Ny; j++)
                    outWriter.WriteLine((dx * i + " " + dy * j + " " + P0[i, j]).Replace(',', '.'));
        }
        public void SaveCurrentValuesVx(String fileName)
        {
            using StreamWriter outWriter = new StreamWriter(fileName, false, System.Text.Encoding.Default);

            for (int i = 0; i < Nx; i++)
                for (int j = 0; j < Ny; j++)
                    outWriter.WriteLine((dx * i + " " + dy * j + " " + Vx[i, j]).Replace(',', '.'));
        }
        public void SaveCurrentValuesVy(String fileName)
        {
            using StreamWriter outWriter = new StreamWriter(fileName, false, System.Text.Encoding.Default);

            for (int i = 0; i < Nx; i++)
                for (int j = 0; j < Ny; j++)
                    outWriter.WriteLine((dx * i + " " + dy * j + " " + Vy[i, j]).Replace(',', '.'));
        }

        public double GettCurrent()
        {
            return tCurrent;
        }

        #endregion
    }
}
