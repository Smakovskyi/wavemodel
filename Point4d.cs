using System;
using System.Collections.Generic;
using System.Text;

namespace wavemodel
{
    public class Point4d
    {

        public double x;
        public double y;
        public double z;
        public double w;

        public Point4d()
        {
        }

        public Point4d(double x, double y, double z, double w)
        {
            this.x = x;
            this.y = y;
            this.z = z;
            this.w = w;
        }
    }
}
