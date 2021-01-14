using System;
using System.IO;

namespace Com_Methods
{
    //интерфейс вектора
    public interface IVector
    {
        //размер
        int N { set; get; }
    }

    //класс вектора: 
    //атрибут Serializable разрешает сериализацию экзепляров данного класса
    //сериализация будет нужна при тестировании ваших программ
    [Serializable]
    public class Vector : IVector
    {
        //размер вектора
        public int N { set; get; }
        public double[] Elem { set; get; }

        //конструктор по умолчанию
        public Vector()
        {
        }

        //конструктор нуль-вектора по размеру n
        public Vector(int n)
        {
            N = n;
            Elem = new double[n];
        }

        //конструктор вектора по бинарному файлу
        public Vector(String PATH)
        {
            //чтение размера системы
            using (var Reader = new BinaryReader(File.Open(PATH + "Size.bin", FileMode.Open)))
            {
                N = Reader.ReadInt32();
            }

            Elem = new double[N];

            //чтение вектора правой части
            using (var Reader = new BinaryReader(File.Open(PATH + "F.bin", FileMode.Open)))
            {
                // считываем данные из файла
                try
                {
                    for (int j = 0; j < N; j++)
                    {
                        Elem[j] = Reader.ReadDouble();
                    }
                }
                catch { throw new Exception("Vector: data file is not correct..."); }
            }
        }

        //редукция размера вектора
        public void Size_Reduction(int New_N)
        {
            N = New_N;
            var Rows = Elem;
            Array.Resize<double>(ref Rows, N);
        }

        //умножение на скаляр с выделением памяти под новый вектор
        public static Vector operator *(Vector T, double Scal)
        {
            Vector RES = new Vector(T.N);

            for (int i = 0; i < T.N; i++)
            {
                RES.Elem[i] = T.Elem[i] * Scal;
            }
            return RES;
        }

        //умножение на скаляр, результат записывается в тот же вектор
        public void Dot_Scal (double Scal)
        {
            for (int i = 0; i < N; i++)
            {
                Elem[i] = Elem[i] * Scal;
            }
        }

        //произведение векторов == матрица
        public Matrix Multiply(Vector V2)
        {
            var RES = new Matrix(N, V2.N);

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < V2.N; j++)
                {
                    RES.Elem[i][j] = Elem[i] * V2.Elem[j];
                }
            }
            return RES;

        }

        //скалярное произведение векторов
        public static double operator *(Vector V1, Vector V2)
        {
            if (V1.N != V2.N) throw new Exception("V1 * V2: dim(vector1) != dim(vector2)...");

            Double RES = 0.0;

            for (int i = 0; i < V1.N; i++)
            {
                RES += V1.Elem[i] * V2.Elem[i];
            }
            return RES;
        }

        //сумма векторов с выделением памяти под новый вектор
        public static Vector operator +(Vector V1, Vector V2)
        {
            if (V1.N != V2.N) throw new Exception("V1 + V2: dim(vector1) != dim(vector2)...");
            Vector RES = new Vector(V1.N);

            for (int i = 0; i < V1.N; i++)
            {
                RES.Elem[i] = V1.Elem[i] + V2.Elem[i];
            }
            return RES;
        }

        //сумма векторов без выделения памяти под новый вектор
        public void Add (Vector V2)
        {
            if (N != V2.N) throw new Exception("V1 + V2: dim(vector1) != dim(vector2)...");

            for (int i = 0; i < N; i++)
            {
                Elem[i] += V2.Elem[i];
            }
        }

        //копирование вектора V2
        public void Copy (Vector V2)
        {
            if (N != V2.N) throw new Exception("Copy: dim(vector1) != dim(vector2)...");
            for (int i = 0; i < N; i++) Elem[i] = V2.Elem[i];
        }

        //копирование вектора-столбца матрицы с номером NUM
        public void Copy_Column (Matrix T, int NUM)
        {
            if (N != T.M) throw new Exception("Copy: dim(vector) != dim(matrix)...");
            for (int i = 0; i < N; i++) Elem[i] = T.Elem[i][NUM];
        }

        //норма вектора: ключевое слово this ссылается на текущий экземпляр класса
        public double Norm()
        {
            return Math.Sqrt(this * this);
        }

        //нормирование веткора
        public void Normalizing()
        {
            double norma = Norm();
            if (norma < CONST.EPS) throw new Exception("Error in Vector.Normalizing: ||V|| = 0...");
            for (int i = 0; i < N; i++) Elem[i] /= norma;
        }

        //вывод вектора на консоль
        public void Console_Write_Vector ()
        {
            for (int i = 0; i < N; i++) Console.WriteLine(Elem[i]);
        }
    }
}