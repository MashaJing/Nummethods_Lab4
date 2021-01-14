using System;
using System.IO;
using ExcelDataReader;
using Excel = Microsoft.Office.Interop.Excel;

namespace Com_Methods
{
    /*
        установка пакета для чтения данных:
        1. пройдите в "Средства - Диспетчер пакетов NuGet - Управление пакетами NuGet для решения"
        2. на вкладке "Обзор" введите в поле поиска ExcelDataReader и нажмите "Поиск"
        3. в выпавшем меню найдите ExcelDataReader (должен быть самым первым) и нажмите по продукту левой кнопкой мыши
        4. в меню справа поставьте галочку напротив вашего проекта, нажмите "установить"
        
        установка пакета для записи данных:
        1. пройдите в "Средства - Диспетчер пакетов NuGet - Управление пакетами NuGet для решения"
        2. на вкладке "Обзор" введите в поле поиска Microsoft.Office.Interop.Excel и нажмите "Поиск"
        3. в выпавшем меню найдите Microsoft.Office.Interop.Excel и нажмите по продукту левой кнопкой мыши
        4. в меню справа поставьте галочку напротив вашего проекта, нажмите "установить"
    
        ***************** ОСОБО ОБРАТИТЬ ВНИМАНИЕ! ***************************

        для работы с файлами требуется указать рабочий каталог:
        1. с помощью правой кнопки мыши перейдите на вкладку "Свойства проекта"
        2. на вкладке "Сборка" найдите поле "Выходные данные"
        3. в поле "Выходной путь" с помощью кнопки "Обзор" задайте выходной путь для чтения и записи файлов для всех конфигураций
    */

    /// <summary>
    /// класс, который позволяет считать/записать таблицу Excel
    /// </summary>
    class Excel_Reader_Writer
    {
        //чтение матрицы из файла Excel: PATH - путь к файлу, например, "Data\\Matrix.xls" (каталог Data предполагается созданным)
        public Matrix Read_Matrix (string PATH)
        {
            //результат чтения
            Matrix RES = new Matrix();

            //чтение данных из excel-файла: потоку указываем путь, метод открытия файла, уровень доступа  
            using (var stream = File.Open(Path.Combine(PATH), FileMode.Open, FileAccess.Read))
            {
                //читатель для файлов с расширением *.xlsx
                var excelReader = ExcelReaderFactory.CreateReader(stream);

                //число строк и столбцов в матрице
                RES.M = excelReader.RowCount;
                RES.N = excelReader.FieldCount;

                //инициализируем хранилище элементов
                RES.Elem = new double[RES.M][];
                for (int i = 0; i < RES.M; i++) 
                {
                    RES.Elem[i] = new double[RES.N];

                    //считываем очередную строку из книги Excel
                    excelReader.Read();

                    //цикл по столбцам книги
                    for (int j = 0; j < RES.N; j++) { RES.Elem[i][j] = excelReader.GetDouble(j); }
                }

                //после завершения чтения освобождаем ресурсы
                excelReader.Close();
            }
            return RES;
        }

        //запись данных в файл Excel: M - матрица
        public int Write_Matrix (Matrix T, string Name = "Matrix")
        {
            //приложение, которое откроет excel-файл по завершению алгоритма программы
            Excel.Application Excel_App = new Excel.Application();

            //создаём экземпляр рабочей книги Excel
            Excel.Workbook Work_Book = Excel_App.Workbooks.Add();

            //создаём экземпляр листа Excel: 1 - номер листа из списка листов (если лист будет не один в книге)
            Excel.Worksheet Work_Sheet = (Excel.Worksheet)Work_Book.Worksheets.get_Item(1);

            //указываем имя листа книги
            Work_Sheet.Name = Name;
           
            //заполняем строки таблицы числами
            for (int i = 1; i <= T.M; i++)
            {
                for (int j = 1; j <= T.N; j++)
                {
                    Work_Sheet.Cells[j][i] = T.Elem[i - 1][j - 1];
                }
            }

            //открываем созданный excel-файл
            //Excel_App.Visible = true;
            //Excel_App.UserControl = true;

            return 0;
        }
    }
}
