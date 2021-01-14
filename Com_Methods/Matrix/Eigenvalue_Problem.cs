using System;
using System.Collections.Generic;

namespace Com_Methods
{
    class Eigenvalue_Problem
    {
        public static double Rayleigh_Shift(Matrix A)
        {
            return A.Elem[A.N - 1][A.N - 1];
        }
        public static double Wilkinson_Shift(Matrix A)
        {
            int N = A.N;
            double a = A.Elem[N - 2][N - 2],
                   b = A.Elem[N - 1][N - 1],
                   c = A.Elem[N - 1][N - 2],
                   d = A.Elem[N - 2][N - 1];
            double D = (a + b) * (a + b) - 4 * (a * b - c * d);
            if (D < 0) throw new Exception("The matrix has complex eigen value...");
            return ((a + b) + Math.Sqrt(D)) * 0.5;
        }

        public static void Shift(Matrix A, double Shift)
        {
            for (int i = 0; i < A.M; i++) A.Elem[i][i] += Shift;
        }

    }
}
