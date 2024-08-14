#include <iostream>
#include <locale.h>
#include <cmath>
#include "omp.h"

using namespace std; // пока можно так

double modul_in3(const double* initial, const int num, const int i, const int j)
{
    double xi = initial[i * num + 1];
    double xj = initial[j * num + 1];
    double yi = initial[i * num + 2];
    double yj = initial[j * num + 2];
    double zi = initial[i * num + 3];
    double zj = initial[j * num + 3];
    double c = (xi - xj) * (xi - xj) + (yi - yj) * (yi - yj) + (zi - zj) * (zi - zj);
    return (pow(c, 1.5));
}

// расчет ускорения тела
double accelpart_xyz(const int type, const double* initial, const int num, const int i, const int j, const double e)
{
    double G = 6.67e-11;
    double typei = initial[i * num + type];
    double typej = initial[j * num + type];
    double c = -(G * initial[j * num + 0] * (typei - typej)) / (max(modul_in3(initial, num, i, j), e * e * e));
    return (c);
}

double coord_next_type(const int type, const double* initial, const int num, const int i, const double dt)
{
    double c = initial[i * num + type] + initial[i * num + type + 3] * dt;
    return (c);
}

double velo_next_type(const int type, const double* initial, const int num, const int i, const double dt, const double accel)
{
    return (initial[i * num + type + 3] + accel * dt);
}

int main()
{
    double e = 1.0e-10; // очень малое число, нужно для замены нуля на него в операции деления
    setlocale(LC_ALL, "Russian");
    int num = 7; // Структура массива для каждого тела [mass, x, y, z, Vx, Vy, Vz]
    double ax, ay, az;//компоненты ускорения
    int num_bodies = 10000; //количество небесных тел
    double* initial = new double[num * num_bodies]; //создание двумерного массива (развернутого в одномерный) исходных данных
    double dt = 1e-5;//шаг интегрирования
    double tk = 5e-5; //конечное время интегрирования
    
    // параллельное задание рандомных начальных условий
    #pragma omp parallel for 
    for (int i = 0; i < num_bodies; ++i)//Ввод данных из файла
    {
        for (int j = 0; j < 7; ++j)
        {
            initial[i * num + j] = rand() * 1.0 / 10e+5;
        }
    }

    double t1 = omp_get_wtime(); // отсечка начального момента времени
    for (int t = 1; (t * dt <= tk); ++t)// Цикл по времени. Распараллелить по нему нельзя
    {
        {
            // параллельно для каждого тела интегрируется уравнеие движения
            #pragma omp parallel for private(ax,ay,az)
            for (int i = 0; i < num_bodies; ++i)
            {
                ax = 0.0;
                ay = 0.0;
                az = 0.0;
                for (int j = 0; j < num_bodies; ++j)
                {
                    ax += accelpart_xyz(1, initial, num, i, j, e);
                    ay += accelpart_xyz(2, initial, num, i, j, e);
                    az += accelpart_xyz(3, initial, num, i, j, e);
                }
                for (int j = 1; j <= 3; ++j)
                {
                    initial[i * num + j] = coord_next_type(j, initial, num, i, dt);
                }
                initial[i * num + 3 + 1] = velo_next_type(1, initial, num, i, dt, ax);
                initial[i * num + 3 + 2] = velo_next_type(2, initial, num, i, dt, ay);
                initial[i * num + 3 + 3] = velo_next_type(3, initial, num, i, dt, az);
            }
        }
    }
    double t2 = omp_get_wtime(); // отсечка конечного времени расчета
    cout << "Время затраченное на один цикл по времени - " << (t2 - t1) / (tk / dt);
    delete[] initial;
}
