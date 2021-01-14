using System;

namespace Com_Methods
{
    /// <summary>
    /// преобразования Хаусхолдера
    /// </summary>
    class Householder_Transformation
    {
        /// <summary>
        /// реализация процедуры ортогонализации по методу отражений Хаусхолдера
        /// </summary>
        /// <param name="A - исходная матрица"></param>
        /// <param name="Q - ортогональная матрица преобразований"></param>
        /// <param name="R - результат"></param>
        public static void Householder_Orthogonalization(Matrix A, Matrix Q, Matrix R)
        {
            //инициализация вектора отражения
            Vector p = new Vector(A.M);

            //вспомогательные переменные
            double s, beta, mu;

            //запись матрицы А в R
            R.Copy(A);

            //алгоритм отражений Хаусхолдера
            for (int i = 0; i < R.N - 1; i++)
            {
                //находим квадрат нормы столбца для обнуления
                s = 0;
                for (int I = i; I < R.M; I++) s += Math.Pow(R.Elem[I][i], 2);

                //если есть ненулевые элементы под диагональю, тогда 
                //квадрат нормы вектора для обнуления не совпадает с квадратом диагонального элемента
                if (Math.Sqrt(Math.Abs(s - R.Elem[i][i] * R.Elem[i][i])) > CONST.EPS)
                {
                    //выбор знака слагаемого beta = sign(-x1)
                    if (R.Elem[i][i] < 0) beta = Math.Sqrt(s);
                    else beta = -Math.Sqrt(s);

                    //вычисляем множитель в м.Хаусхолдера mu = 2 / ||p||^2
                    mu = 1.0 / beta / (beta - R.Elem[i][i]);

                    //формируем вектор p
                    for (int I = 0; I < R.M; I++) { p.Elem[I] = 0; if (I >= i) p.Elem[I] = R.Elem[I][i]; }

                    //изменяем диагональный элемент
                    p.Elem[i] -= beta;

                    //вычисляем новые компоненты матрицы A = Hm * H(m-1) * ... * A
                    for (int m = i; m < R.N; m++)
                    {
                        //произведение S = At * p
                        s = 0;
                        for (int n = i; n < R.M; n++) { s += R.Elem[n][m] * p.Elem[n]; }
                        s *= mu;
                        //A = A - 2 * p * (At * p)^t / ||p||^2
                        for (int n = i; n < R.M; n++) { R.Elem[n][m] -= s * p.Elem[n]; }
                    }

                    //вычисляем новые компоненты матрицы Q = Q * H1 * H2 * ...
                    for (int m = 0; m < R.M; m++)
                    {
                        //произведение Q * p
                        s = 0;
                        for (int n = i; n < R.M; n++) { s += Q.Elem[m][n] * p.Elem[n]; }
                        s *= mu;
                        //Q = Q - p * (Q * p)^t
                        for (int n = i; n < R.M; n++) { Q.Elem[m][n] -= s * p.Elem[n]; }
                    }
                }
            }
        }

        /// <summary>
        /// метод зануления i-го столбца с j-ой позиции
        /// <param name="A - трансформирующаяся матрица"></param>
        /// <param name="U - ортогональная матрица преобразования (инициализирована единицами на диагонали)"></param>
        /// <param name="i - номер столбца"></param>
        /// <param name="j - позиция в столбце, с которой будет выполнено зануление"></param>
        /// </summary>
        public static void Column_Transformation(Matrix A, Matrix U, int i, int j)
        {
            //вектор отражения
            Vector p = new Vector(A.M);

            //вспомогательные переменные
            double s, beta, mu;

            //находим квадрат нормы столбца для обнуления
            s = 0;
            for (int I = j; I < A.M; I++) s += Math.Pow(A.Elem[I][i], 2);

            //если ненулевые элементы под диагональю есть: 
            //квадрат нормы вектора для обнуления не совпадает с квадратом зануляемого элемента
            if (Math.Sqrt(Math.Abs(s - A.Elem[j][i] * A.Elem[j][i])) > CONST.EPS)
            {
                //выбор знака слагаемого beta = sign(-x1)
                if (A.Elem[j][i] < 0) beta = Math.Sqrt(s);
                else beta = -Math.Sqrt(s);

                //вычисляем множитель в м.Хаусхолдера mu = 2 / ||p||^2
                mu = 1.0 / beta / (beta - A.Elem[j][i]);

                //формируем вектор p
                for (int I = 0; I < A.M; I++) { p.Elem[I] = 0; if (I >= j) p.Elem[I] = A.Elem[I][i]; }

                //изменяем элемент, с которого начнётся обнуление
                p.Elem[j] -= beta;

                //вычисляем новые компоненты матрицы A = ... * U2 * U1 * A
                for (int m = 0; m < A.N; m++)
                {
                    //произведение S = St * p
                    s = 0;
                    for (int n = j; n < A.M; n++) { s += A.Elem[n][m] * p.Elem[n]; }
                    s *= mu;
                    //S = S - 2 * p * (St * p)^t / ||p||^2
                    for (int n = j; n < A.M; n++) { A.Elem[n][m] -= s * p.Elem[n]; }
                }

                //вычисляем новые компоненты матрицы U = ... * H2 * H1 * U
                for (int m = 0; m < A.M; m++)
                {
                    //произведение S = Ut * p
                    s = 0;
                    for (int n = j; n < A.M; n++) { s += U.Elem[m][n] * p.Elem[n]; }
                    s *= mu;
                    //U = U - 2 * p * (Ut * p)^t / ||p||^2
                    for (int n = j; n < A.M; n++) { U.Elem[m][n] -= s * p.Elem[n]; }
                }
            }
        }

        //-----------------------------------------------------------------------------

        /// <summary>
        /// метод зануления i-ой строки с j-ой позиции
        /// <param name="A - трансформирующаяся матрица"></param>
        /// <param name="V - ортогональная матрица преобразования (инициализирована единицами на диагонали)"></param>
        /// <param name="i - номер строки"></param>
        /// <param name="j - позиция в строке, с которой будет выполнено зануление"></param>
        /// </summary>
        public static void Row_Transformation(Matrix A, Matrix V, int i, int j)
        {
            //вектор отражения
            Vector p = new Vector(A.N);

            //вспомогательные переменные
            double s, beta, mu;

            //находим квадрат нормы строки для обнуления
            s = 0;
            for (int I = j; I < A.N; I++) s += Math.Pow(A.Elem[i][I], 2);

            //если ненулевые элементы под диагональю есть: 
            //квадрат нормы вектора для обнуления не совпадает с квадратом зануляемого элемента
            if (Math.Sqrt(Math.Abs(s - A.Elem[i][j] * A.Elem[i][j])) > CONST.EPS)
            {
                //выбор знака слагаемого beta = sign(-x1)
                if (A.Elem[i][j] < 0) beta = Math.Sqrt(s);
                else beta = -Math.Sqrt(s);

                //вычисляем множитель в м.Хаусхолдера mu = 2 / ||p||^2
                mu = 1.0 / beta / (beta - A.Elem[i][j]);

                //формируем вектор p
                for (int I = 0; I < A.N; I++) { p.Elem[I] = 0; if (I >= j) p.Elem[I] = A.Elem[i][I]; }

                //изменяем диагональный элемент
                p.Elem[j] -= beta;

                //вычисляем новые компоненты матрицы A = A * H1 * H2 ...
                for (int m = 0; m < A.M; m++)
                {
                    //произведение A * p
                    s = 0;
                    for (int n = j; n < A.N; n++) { s += A.Elem[m][n] * p.Elem[n]; }
                    s *= mu;
                    //A = A - p * (A * p)^t
                    for (int n = j; n < A.N; n++) { A.Elem[m][n] -= s * p.Elem[n]; }
                }

                //вычисляем новые компоненты матрицы V = V * H1 * H2 * ...
                for (int m = 0; m < A.N; m++)
                {
                    //произведение V * p
                    s = 0;
                    for (int n = j; n < A.N; n++) { s += V.Elem[m][n] * p.Elem[n]; }
                    s *= mu;
                    //V = V - p * (V * p)^t
                    for (int n = j; n < A.N; n++) { V.Elem[m][n] -= s * p.Elem[n]; }
                }
            }
        }
    }
}
