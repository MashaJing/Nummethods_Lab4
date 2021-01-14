using System;
using System.Collections.Generic;

namespace Com_Methods
{  
    //SVD на базе ортогональных преобразований
    class SVD
    {
        //левые сингулярные векторы
        public Matrix U { set; get; }
        //правые сингулярные векторы
        public Matrix V { set; get; }
        //сингулярные числа
        public Matrix Sigma { set; get; }
        //результат перемножения U*Sigma*Vt
        public Matrix RES { set; get; }
        
        //перечисление методов диагонализации на 2 этапе
        public enum Algorithm
        {
            Householder = 1,
            Givens
        }
        //---------------------------------------------------------------------------------------------

        /// <summary>
        /// конструктор
        /// </summary>
        /// <param name="A - матрица для SVD"></param>
        public SVD(Matrix A)
        {
            //вызов метода, который выполнит построение SVD
            Start_SVD(A);
        }

        public SVD(Matrix A, double Eps)
        {
            //вызов метода, который выполнит построение SVD алгоритмом
            //исчерпывания с заданной точностью Eps
            Start_SVD_exhaust(A, Eps);
        }
        
        //---------------------------------------------------------------------------------------------

        /// <summary>
        /// проверка на положительность сингулярных чисел
        /// </summary>
        private void Check_Singular_Values()
        {
            //наименьшее измерение
            int Min_Size = Math.Min(Sigma.M, Sigma.N);

            //проверка сингулярных чисел на положительность
            for (int i = 0; i < Min_Size; i++)
            {   
                if (Sigma.Elem[i][i] < 0)
                {
                    Sigma.Elem[i][i] = -Sigma.Elem[i][i];

                    for (int j = 0; j < U.M; j++)
                        U.Elem[j][i] = -U.Elem[j][i];
                }
            }
        }

        //---------------------------------------------------------------------------------------------

        /// <summary>
        /// сортировка сингулярных чисел
        /// </summary>
        private void Sort_Singular_Values()
        {
            //наименьшее измерение
            int Min_Size = Math.Min(Sigma.M, Sigma.N);

            //сортировка сингулярных чисел
            for (int I = 0; I < Min_Size; I++)
            {
                var Max_Elem = Sigma.Elem[I][I];
                int Index = I;
                for (int i = I + 1; i < Min_Size; i++)
                {
                    if (Sigma.Elem[i][i] > Max_Elem)
                    {
                        Max_Elem = Sigma.Elem[i][i];
                        Index = i;
                    }
                }
                //найден наибольший элемент
                if (I != Index)
                {
                    Sigma.Elem[Index][Index] = Sigma.Elem[I][I];
                    Sigma.Elem[I][I] = Max_Elem;
                    U.Column_Transposition(I, Index);
                    V.Column_Transposition(I, Index);
                }
            }
        }

        //---------------------------------------------------------------------------------------------

        /// <summary>
        /// редукция SVD 
        /// </summary>
        /// <param name="Reduction - порог отбрасывания сингулярных чисел"></param>
        public void Reduction_SVD(double Reduction)
        {
            //наименьшее измерение
            int Min_Size = Math.Min(Sigma.M, Sigma.N);

            //проверка на возможность редукции по сингулярным числам
            for (int i = 0; i < Min_Size; i++)
            {
                if (Math.Abs(Sigma.Elem[i][i]) < Reduction)
                {
                    Min_Size = i;
                    break;
                }
            }
            //редукция размерности матриц
            Sigma.Size_Reduction(Min_Size, Min_Size);
            U.Size_Reduction(U.M, Min_Size);
            V.Size_Reduction(V.M, Min_Size);
            

            RES = new Matrix(U.M, V.M);
            RES = U * Sigma * V.Transpose_Matrix();
        }
        
        //---------------------------------------------------------------------------------------------
        /// SVD-алгоритм на базе алгоритма с исчерпыванием
        /// 
        public void Start_SVD_exhaust(Matrix A, double Eps)
        {
            //инициализация матрицы левых сингулярных векторов 
            U = new Matrix(A.M, A.M);
            
            //матрица сингулярных чисел
            Sigma = new Matrix(A.M, A.N);

            //инициализация матрицы правых сингулярных векторов 
            V = new Matrix(A.N, A.N);

            //инициализация матриц для SVD
            for (int i = 0; i < A.M; i++)
            {
                U.Elem[i][i] = 0.0;
                for (int j = 0; j < A.N; j++) Sigma.Elem[i][j] = 0.0;
            }
            for (int i = 0; i < A.N; i++) V.Elem[i][i] = 0.0;
            
            //инициализация вектор-столбца u на k-ой итерации
            var Uk = new Vector(U.M);

            //инициализация вектор-строки v на k-ой итерации
            var Vk = new Vector(V.N);
            
            int k = 0;
            double S = 0;
            
            //**************** Алгоритм с исчерпыванием *************************
            
            do
            {

            double max_norm = 0;
            int Index = 0;
            
            //выбираем строку с наибольшей Евклидовой нормой
            for (int i = 0; i < A.M; i++)
            {
                S = 0;
                for (int j=0; j<A.N; j++) S += Math.Pow(A.Elem[i][j], 2);
                S = Math.Sqrt(S);

                if (max_norm < S)
                {
                  Index = i;
                  max_norm = S;
                }
            }
            
            for (int i = 0; i < Vk.N; i++) Vk.Elem[i] = A.Elem[Index][i];
            if (S < Eps) return;
            
            Vector Vk_new = new Vector(Vk.N);
            Vector Uk_new = new Vector(Uk.N);
            Vk_new = Vk;
            Uk_new = Uk;

            do
            {
            Vk = Vk_new;
            Uk = Uk_new;
            //найдём вектор-столбец u для k-ой итерации
            double VecPow2 = 1/(Vk*Vk);
            Uk_new = A*Vk*VecPow2;

            //найдём вектор-строку v для k-ой итерации
            VecPow2 = 1/(Uk_new*Uk_new);
            
            Vk_new = A.Multiplication_Trans_Matrix_Vector(Uk_new)*VecPow2;

            //если вектор не меняется (модуль разности норм 
            //векторов не превосходит заданной точности), резултат
            //сошёлся к доминирующему собственному вектору
            }
            while (Math.Abs(Vk.Norm() - Vk_new.Norm()) > Eps || (Math.Abs(Uk.Norm() - Uk_new.Norm()) > Eps));

            //найдём k-ое сингулярное число
            Sigma.Elem[k][k] = Uk.Norm()*Vk.Norm();
            
            //Вставляем в матрицы V и U нормированные векторы
            Uk.Normalizing();
            Vk.Normalizing();
            
            //добавляем на соответствующие позиции в матрицы U и V
            U.Vector_To_Column(Uk, k);
            V.Vector_To_Column(Vk, k); //Чтобы получить V, а не VT
            
            A = A + (Uk.Multiply(Vk))*(-Sigma.Elem[k][k]);
            
            k++;
            }
            while (k < Rank());

            //убираем отрицательные сингулярные числа
            Check_Singular_Values();

            //сортируем по возрастанию сингулярные числа
           // Sort_Singular_Values(); //???

            RES = new Matrix(U.M, V.M);
            RES = U * Sigma * V.Transpose_Matrix();

            
        }

        
        //---------------------------------------------------------------------------------------------

        /// <summary>
        /// двухэтапный SVD-алгоритм
        /// </summary>
        /// <param name="A - матрица для SVD"></param>
        public void Start_SVD(Matrix A)
        {
            //наименьшее измерение
            int Min_Size = Math.Min(A.M, A.N);

            //размеры нижней и верхней внешних диагоналей
            int Up_Size = Min_Size - 1, Down_Size = Min_Size - 1;

            //инициализация матрицы левых сингулярных векторов 
            U = new Matrix(A.M, A.M);
            
            //матрица сингулярных чисел
            Sigma = new Matrix(A.M, A.N);

            //инициализация матрицы правых сингулярных векторов 
            V = new Matrix(A.N, A.N);

            //инициализация матриц для SVD
            for (int i = 0; i < A.M; i++)
            {
                U.Elem[i][i] = 1.0;
                for (int j = 0; j < A.N; j++) Sigma.Elem[i][j] = A.Elem[i][j];
            }
            for (int i = 0; i < A.N; i++) V.Elem[i][i] = 1.0;

            //**************** этап I: бидиагонализация *************************
            
            for (int i = 0; i < Min_Size - 1; i++)
            {
                Householder_Transformation.Column_Transformation(Sigma, U, i, i);
                Householder_Transformation.Row_Transformation(Sigma, V, i, i + 1);
            }

            //ситуация M > N - строк больше => дополнительное умножение слева 
            if (A.M > A.N)
            {
                Householder_Transformation.Column_Transformation(Sigma, U, A.N - 1, A.N - 1);
                //нижняя побочная диагональ длиннее на 1
                Down_Size += 1;
            }

            //ситуация M < N - столбцов больше => дополнительное умножение справа
            if (A.M < A.N)
            {
                Householder_Transformation.Row_Transformation(Sigma, V, A.M - 1, A.M);
                //верхняя побочная диагональ длиннее на 1
                Up_Size += 1;
            }

            Sigma.Console_Write_Matrix();

            //**************** этап II: преследование ************
            //********* приведение к диагональному виду **********

            //для хранения изменяющихся элементов верхней диагонали
            var Up = new double[Up_Size];
            var Down = new double[Down_Size];
            //число неизменившихся элементов над главной диагональю
            int CountUpElements;

            //процедура преследования
            do
            {
                CountUpElements = 0;

                //обнуление верхней диагонали
                for (int i = 0; i < Up_Size; i++)
                {
                    if (Math.Abs(Up[i] - Sigma.Elem[i][i + 1]) > CONST.EPS)
                    {
                        Up[i] = Sigma.Elem[i][i + 1];
                        Givens_Transformation.Row_Transformation(Sigma, V, i, i);
                    }
                    else
                        CountUpElements++;
                }

                //обнуление нижней диагонали
                for (int i = 0; i < Down_Size; i++)
                {
                    if (Math.Abs(Down[i] - Sigma.Elem[i + 1][i]) > CONST.EPS)
                    {
                        Down[i] = Sigma.Elem[i + 1][i];
                        Givens_Transformation.Column_Transformation(Sigma, U, i, i);
                    }
                }
            }
            while (CountUpElements != Up_Size);

            //убираем отрицательные сингулярные числа
            Check_Singular_Values();
            //сортируем по возрастанию сингулярные числа
            Sort_Singular_Values();

            RES = new Matrix(U.M, V.M);
            RES = U * Sigma * V.Transpose_Matrix();

        }

        //---------------------------------------------------------------------------------------------

        //ранг матрицы
        public int Rank()
        {
            return Sigma.M;
        }

        //---------------------------------------------------------------------------------------------

        //модуль определителя матрицы
        public double Abs_Det()
        {
            //ранг матрицы
            int Size = Rank();

            if (Size == 0) throw new Exception("Error in SVD.Rank: SVD is not built ...");
            double det = 1;
            for (int i = 0; i < Size; i++)
                det *= Sigma.Elem[i][i];
            return det;
        }

        //---------------------------------------------------------------------------------------------

        //число обусловленности матрицы
        public double Cond()
        {
            //ранг матрицы
            int Size = Rank();

            if (Size == 0) throw new Exception("Error in SVD.Rank: SVD is not built ...");
            return Sigma.Elem[0][0] / Sigma.Elem[Size - 1][Size - 1];
        }

        //---------------------------------------------------------------------------------------------

        //нормальное псевдорешение СЛАУ Ax = F => x = V * S^(-1) * Ut * F
        public Vector Start_Solver(Vector F)
        {
            //ранг матрицы
            int Size = Rank();

            //int Min_Size = Math.Min(U.M, U.N);

            if (Size == 0) throw new Exception("Error in SVD.Rank: SVD is not built ...");

          //  if (Size != F.N) throw new Exception("Error in SVD.Rank: invalid size of F ...");
            
            //UtF = Ut * F
            var UtF = U.Multiplication_Trans_Matrix_Vector(F);

            //UtF = S^(-1) * Ut * F
            for (int i = 0; i < UtF.N; i++) UtF.Elem[i] /= Sigma.Elem[i][i];

            //Res = V * S^(-1) * Ut * F
            var Res = V * UtF;
            return Res;
        }

        public double Quality(Matrix A)
        {
            double q = 0;
            double q_max = 0;

            for (int i=0; i<A.M; i++)
                for (int j=0; j<A.N; j++)
                {
                    if  (A.Elem[i][j] != 0) 
                    {q = Math.Abs(A.Elem[i][j] - RES.Elem[i][j])/A.Elem[i][j];
                    if (q_max < q) q_max = q;}
                }
            return q_max;
        }

    }
}
