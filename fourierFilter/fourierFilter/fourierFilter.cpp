// fourierFilter.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"
#include <iostream>
#include "fft.h"
#include <ctime>
#include <chrono>
using namespace std;
void run1() {

	int n = 256;
	std::cout << "Hello World!\n";
	double* tem = new double[2 * n];
	double* res = new double[2 * n];
	double* rev = new double[2 * n];
	for (int i = 0; i < n; i++) {
		tem[2*i] = sin(PI*double(i) / double(n / 2));
		tem[2 * i + 1] = 0.0;
	}
	/*for (int i = 0; i < n; i++)
		cout << i << tem[i] << endl;
	cout << " result" << endl;*/
	SerialFFT1D(tem, res, n, 1, 0);
	SerialInverseFFT1D(res, rev, n, 1, 0);
	double re, im;
	/*for (int i = 0; i < n; i++) {
		re = (abs(res[i].real()) < 0.000001) ? 0.0 : res[i].real();
		im = (abs(res[i].imag()) < 0.000001) ? 0.0 : res[i].imag();
		cout << i << "  " << re << "   " << im << endl;
	}*/
	cout << " reverse" << endl;
	for (int i = 0; i < n; i++)
		cout << i << "  " << tem[i] << "  " << rev[i] << endl;
}

void run2() {
	auto begin = std::chrono::steady_clock::now();
	int n = 1024;
	std::cout << "Hello World!\n";
	double* tem = new double[2 * n* 2 * n];
	double* res = new double[2 * n* 2 * n];
	double* rev = new double[2 * n* 2 * n];

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			tem[2 * (i*n + j)] = sin(PI*double(i) / double(n / 2))*sin(PI*double(j) / double(n / 2));
			tem[2 * (i*n + j)+1] = 0;
		}
	}
	/*for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cout << tem[i*n + j] << " ";
		cout << endl;
	}*/
	SerialFFT2D(tem, res, n, n);
	SerialInverseFFT2D(res, rev,n, n);
	
	auto end = std::chrono::steady_clock::now();
	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	auto elapsed_s = std::chrono::duration_cast<std::chrono::seconds>(end - begin);
	cout << "work time: " << elapsed_s.count() <<endl;
	cout << " reverse" << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cout << tem[2*(i*n + j)] <<","<<tem[2*(i*n + j)+1] << "|" << rev[2*(i*n + j)] <<","<<rev[2*(i*n + j)+1] << " " ;
		cout << endl;
	}
	delete[] tem;
	delete[] res;
	delete[] rev;
}


int main()
{
	run2();

}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
