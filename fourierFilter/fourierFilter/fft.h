#pragma once
#include <complex>
#include <vector>
#include <thread>

using namespace std;

#define PI (3.14159265358979323846)

void SerialFFT2D(double* inputSignal, double* outputSignal, int w, int h);

void SerialInverseFFT2D(double* inputSignal, double* outputSignal, int w, int h);
void SerialFFT1D(double* inputSignal, double* outputSignal, int size, int step, int offset);
void SerialInverseFFT1D(double* inputSignal, double* outputSignal, int size, int step, int offset);
