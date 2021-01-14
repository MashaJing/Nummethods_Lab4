using System;
using System.Collections.Generic;
using System.Threading;
using System.IO;

namespace Com_Methods
{
    //интерфейс матрицы
    public interface IMatrix
    {
        //размер матрицы
        int M { set; get; } //строки
        int N { set; get; } //столбцы
    }

    //класс матрицы
    //атрибут Serializable разрешает сериализацию экзепляров данного класса
    [Serializable]
    public class Matrix : IMatrix
    {
        //размер матрицы
        public int M { set; get; }
        public int N { set; get; }
        //элементы матрицы
        public double[][] Elem { set; get; }

        //конструктор по умолчанию
        public Matrix(){}

        //конструктор нуль-матрицы m X n
        public Matrix(int m, int n)
        {
            N = n; M = m;
            Elem = new double[m][];
            for (int i = 0; i < m; i++) Elem[i] = new double[n];
        }

        //редукция размера матрицы
        public void Size_Reduction(int New_M, int New_N)
        {
            M = New_M;
            N = New_N;
            var Rows = Elem;
            Array.Resize<double[]>(ref Rows, M);
            for (int Row = 0; Row < M; Row++)
            {
                Array.Resize<double>(ref Rows[Row], N);
            }
            Elem = Rows;
        }

        //запись вектора в i-ый столбец матрицы
        public void Vector_To_Column(Vector V, int i)
        {
            if (V.N != M)
                throw new Exception("Error in Vector_To_Column: invalid size...");
            if (i < 0 || i>=N)
                throw new Exception("Error in Vector_To_Column: invalid indices...");

            for (int row = 0; row < M; row++)
            {
                Elem[row][i] = V.Elem[row];
            }
        }


        //запись вектора в i-ую строку матрицы
        public void Vector_To_Row(Vector V, int i)
        {
            if (V.N != N)
                throw new Exception("Error in Vector_To_Row: size is not correct...");
            if (i < 0 || i>=M)
                throw new Exception("Error in Vector_To_Row: indices are not correct...");

            for (int col = 0; col < N; col++)
            {
                Elem[i][col] = V.Elem[col];
            }

        }

        //запись i-ого столбца матрицы в отдельный вектор
        public Vector Column_To_Vector(int i)
        {
            if (i < 0 || i >= N)
                throw new Exception("Error in Column_To_Vector: indices are not correct...");

            Vector V = new Vector(M);
            for (int row=0; row<N; row++) V.Elem[row] =Elem[row][i];
            return V;    
        }
        
        //запись i-ой строки матрицы в отдельный вектор
        public Vector Row_To_Vector(int i)
        {
            if (i < 0 || i >= M)
                throw new Exception("Error in Row_To_Vector: indices are not correct...");

            Vector V = new Vector(N);
            for (int col=0; col<N; col++) V.Elem[col] =Elem[i][col];
            return V;    
        }

        //перестановка столбцов с номерами i и j
        public void Column_Transposition(int i, int j)
        {
            if (i < 0 || i > N || j < 0 || j > N) 
                throw new Exception("Error in Column_Transposition: indices are not correct...");
            for (int row = 0; row < M; row++)
            {
                var elem = Elem[row][i];
                Elem[row][i] = Elem[row][j];
                Elem[row][j] = elem;
            }
        }

        //умножение на скаляр с выделением памяти под новую матрицу
        public static Matrix operator *(Matrix T, double Scal)
        {
            Matrix RES = new Matrix(T.M, T.N);

            for (int i = 0; i < T.M; i++)
            {
                for (int j = 0; j < T.N; j++)
                {
                    RES.Elem[i][j] = T.Elem[i][j] * Scal;
                }
            }
            return RES;
        }

        //умножение на скаляр, результат запишется в исходную матрицу
        public void Dot_Scal (double Scal)
        {
            for (int i = 0; i < M; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    Elem[i][j] *= Scal;
                }
            }
        }

        //умножение матрицы на вектор
        public static Vector operator *(Matrix T, Vector V)
        {
            if (T.N != V.N) throw new Exception("M * V: dim(Matrix) != dim(Vector)...");

            Vector RES = new Vector(T.M);

            for (int i = 0; i < T.M; i++)
            {
                for (int j = 0; j < T.N; j++)
                {
                    RES.Elem[i] += T.Elem[i][j] * V.Elem[j];
                }
            }
            return RES;
        }

        //умножение транспонированной матрицы на вектор (умножение вектор-строки на матрицу)
        public Vector Multiplication_Trans_Matrix_Vector (Vector V)
        {
            if (M != V.N) throw new Exception("Mt * V: dim(Matrix) != dim(Vector)...");

            Vector RES = new Vector(N);

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    RES.Elem[i] += Elem[j][i] * V.Elem[j];
                }
            }
            return RES;
        }

        //умножение матрицы на матрицу
        public static Matrix operator *(Matrix T1, Matrix T2)
        {
            if (T1.N != T2.M) throw new Exception("M * M: dim(Matrix1) != dim(Matrix2)...");

            Matrix RES = new Matrix(T1.M, T2.N);

            for (int i = 0; i < T1.M; i++)
            {
                for (int j = 0; j < T2.N; j++)
                {
                    for (int k = 0; k < T1.N; k++)
                    {
                        RES.Elem[i][j] += T1.Elem[i][k] * T2.Elem[k][j];
                    }
                }
            }
            return RES;
        }

        //сложение матриц с выделением памяти под новую матрицу
        public static Matrix operator +(Matrix T1, Matrix T2)
        {
            if (T1.M != T2.M || T1.N != T2.N) throw new Exception("dim(Matrix1) != dim(Matrix2)...");

            Matrix RES = new Matrix(T1.M, T2.N);

            for (int i = 0; i < T1.M; i++)
            {
                for (int j = 0; j < T2.N; j++)
                {
                    RES.Elem[i][j] = T1.Elem[i][j] + T2.Elem[i][j];
                }
            }
            return RES;
        }


        //сложение матриц без выделения памяти под новую матрицу (добавление в ту же матрицу)
        public void Add(Matrix T2)
        {
            if (M != T2.M || N != T2.N) throw new Exception("dim(Matrix1) != dim(Matrix2)...");

            for (int i = 0; i < M; i++)
            {
                for (int j = 0; j < T2.N; j++)
                {
                    Elem[i][j] += T2.Elem[i][j];
                }
            }
        }

        //транспонирование матрицы
        public Matrix Transpose_Matrix()
        {
            var RES = new Matrix(N, M);

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    RES.Elem[i][j] = Elem[j][i];
                }
            }
            return RES;
        }

        //копирование матрицы A
        public void Copy(Matrix A)
        {
            if (N != A.N || M != A.M) throw new Exception("Copy: dim(matrix1) != dim(matrix2)...");
            for (int i = 0; i < M; i++)
                for (int j = 0; j < N; j++)
                    Elem[i][j] = A.Elem[i][j];
        }

        public void Hessenberg_Matrix()
        {
            Vector w = new Vector(M);

            double s, beta, mu;

            //алгоритм отражений Хаусхолдера
            for (int i=0; i < N - 2; i++)
            {
                //находим квадрат нормы столбца для обнуления
                s = 0;
                for (int I = i + 1; I < M; I++) s += Math.Pow(Elem[I][i], 2);

                //если ненулевые элементы под диагональю есть:
                //норма вектора для обнуления не совпала с квадратом диагонального элемента
                if (Math.Sqrt(Math.Abs(s - Elem[i+1][i]* Elem[i + 1][i]))>CONST.EPS)
                {
                    //выбор знака слагаемого beta
                    if (Elem[i + 1][i] < 0) beta = Math.Sqrt(s);
                    else beta = -Math.Sqrt(s);

                    //вычисляем множитель в м. Хаусхолдера mu = 2/ ||w||^2
                    mu = 1.0 / beta / (beta - Elem[i + 1][i]);

                    //формируем вектор w
                    for (int I = 0; I < M; I++) { w.Elem[I] = 0; if (I >= i + 1) w.Elem[I] = Elem[I][i]; }

                    //изменяем диагональный эелемент
                    w.Elem[i + 1] -= beta;

                    //вычисляем новые компоненты матрицы A = Hm * H(m-1) ... * A
                    for (int m = i; m < N; m++)
                    {
                        //произведение S = At * w
                        s = 0;
                        for (int n = i; n < M; n++) { s += Elem[n][m] * w.Elem[n]; }
                        s *= mu;

                        //A = A - 2 * w * (At * w)^t / ||w||^2
                        for (int n = i; n < M; n++) { Elem[n][m] -= s * w.Elem[n]; }
                    }

                    //правое произведение матрицы A = A * H1 * H2 * ...
                    for (int m =0; m < M; m++)
                    {
                        //произведение A * w
                        s = 0;
                        for (int n = 0; n < M; n++) { s += Elem[m][n] * w.Elem[n]; }
                        s *= mu;

                        //A = A - 2 * w * (At * w)^t
                        for (int n = i; n < M; n++) { Elem[m][n] -= s * w.Elem[n]; }
                    }
                }
            }
        }

        public List<double> Eigenvalues(QR_Decomposition.QR_Algorithm Method)
        {
            if (N != M) throw new Exception("Eigenvalues: invalid size of matrix...");

            //создадим копию данной матрицы, чтобы не изменять её
            Matrix T = new Matrix(N, N);
            T.Copy(this);
            List<double> RES = new List<double>();


            //приведём матрицу к верхней форме Хессенберга
             T.Hessenberg_Matrix();
            int t = 1;
            //QR-итерации
            while (T.M != 1)
            {
                Console.WriteLine("Iteration {0} ", t);
                T.Console_Write_Matrix();
                Console.WriteLine(" ");

                for (int i = T.M - 1; i > 0; i--)
                {
                    //если элемент А[i][i-1] == 0, то A[i][i] - собственное значение
                    if (Math.Abs(T.Elem[i][i - 1]) < 1e-6)
                    {
                        RES.Add(T.Elem[i][i]);
                        //исключаем i-ые строку и столбец
                        T.Delete_Row_Column(i);
                        i = T.M;
                    }
                }

                if (T.M == 1) break;

                //прямой сдвиг
                double shift = Eigenvalue_Problem.Rayleigh_Shift(T);
                Eigenvalue_Problem.Shift(T, -shift);

                var QR = new QR_Decomposition(T,  Method);

                T = QR.R * QR.Q;

                Eigenvalue_Problem.Shift(T, shift);
                
                //аварийное завершение цикла
                t++;
                if (t == 100) break;
            }

            //дополняем список последним оставшимся собственным значением
            RES.Add(T.Elem[0][0]);

            //сортируем по возрастанию и формируем результат
            RES.Sort((double X, double Y) => { if (X > Y) return 0; return 1; });

            return RES;
        }

        public void Delete_Row_Column(int i)
        {
            //исключаем i-ые строку и столбец
            for (int I=i; I< M-1; I++)
            {
                Elem[I] = Elem[I + 1];
                //перемещаем столбцы в i-ой строке
                for (int J = i; J < N - 1; J++)
                {
                    Elem[I][J] = Elem[I][J + 1];
                }
            }
            Size_Reduction(M - 1, N - 1);
        }

        //вывод матрицы на консоль
        public void Console_Write_Matrix()
        {
            for (int i = 0; i < M; i++)
            {
                for (int j = 0; j < N; j++)
                    Console.Write(String.Format("{0, -22}", Elem[i][j].ToString("E5")));
                Console.WriteLine();
            }
        }

    }
}