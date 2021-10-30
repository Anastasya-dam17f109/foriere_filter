#include "pch.h"
#include "fft.h"

void BitReversing(double* inputSignal, double* outputSignal, int size, int step, int offset)
{
	int bitsCount = 0;

	//ќпределение количества бит дл€ бит-реверсировани€. ѕолучаетс€, что bitsCount = log2(size).
	for (int tmp_size = size; tmp_size > 1; tmp_size /= 2, bitsCount++);

	//¬ыполнение бит-реверсировани€
	//ind - индекс элемента в  массиве input
	//revInd - соответствующие индексу ind индекс (бит-реверсивный) в массиве output
	for (int ind = 0; ind < size; ind++)
	{
		int mask = 1 << (bitsCount - 1);
		int revInd = 0;

		for (int i = 0; i < bitsCount; i++)
		{
			bool val = ind & mask;
			revInd |= val << i;
			mask = mask >> 1;
		}

		outputSignal[2*revInd*step + 2 * offset] = inputSignal[2 * ind*step + 2 * offset];
		outputSignal[2 * revInd*step + 2 * offset+1] = inputSignal[2 * ind*step + 2 * offset+1];
	}
}


void tread_inverse(double* inputSignal, double* tem, int w, int h, int coeff, int beg, int end, int flag) {
	for (int i = beg; i < end; i++) {
		if (flag == true)
			SerialInverseFFT1D(inputSignal, tem, h, w, coeff*i);
		else
			SerialFFT1D(inputSignal, tem, h, w, coeff*i);
	}

}

void Butterfly(double* signal, double uR, double uI, int offset, int butterflySize, int step, int off, int j)
{
	//complex<double> tem = signal[off + step * (offset + butterflySize)] * u;
	
	double  tem_real = signal[2 * off + 2*step * (j+offset + butterflySize)] * uR - signal[2 * off + 1 + 2 * step * (j + offset + butterflySize)] * uI;
	double  tem_img = signal[2 * off + 2*step * (j+offset + butterflySize)] * uI + signal[2 * off + 1 + 2 * step * (j + offset + butterflySize)] * uR;
	signal[2 * off + 2 * step * (j + offset + butterflySize)] = signal[2 * off + 2 * step * (j + offset)] - tem_real;
	signal[2 * off +1+ 2 * step * (j + offset + butterflySize)] = signal[2 * off +1+ 2 * step * (j+offset)] - tem_img;
	signal[2 * off + 2 * step * (j + offset)] += tem_real;
	signal[2 * off +1+ 2 * step * (j + offset)] += tem_img;
	/*signal[off + step * (offset + butterflySize)] = signal[off + step * (offset)] - tem;
	signal[off + step * (offset)] += tem;*/
}

//void SerialFFTCalculation(double* signal, int first, int size, int step, int offset, bool forward = true)
//{
//	if (size == 1)
//		return;
//
//	double const coeff = 2.0*PI / size;
//
//	SerialFFTCalculation(signal, first, size / 2, step, offset, forward);
//	SerialFFTCalculation(signal, first + size / 2, size / 2, step, offset, forward);
//
//	for (int j = first; j < first + size / 2; j++)
//		if (forward)
//			Butterfly(signal, complex<double>(cos(-j * coeff), sin(-j * coeff)), j, size / 2, step, offset);
//		else
//			Butterfly(signal, complex<double>(cos(j*coeff), sin(j*coeff)), j, size / 2, step, offset);
//}
void SerialFFTCalculation(double* signal, int first, int size, int step, int off, bool forward = true)
{
	int m = 0;
	for (int tmp_size = size; tmp_size > 1; tmp_size /= 2, m++);//size = 2^m

	double*_mas = signal + 2* off;

	for (int p = 0; p < m; p++)
	{
		int butterflyOffset = 1 << (p + 1);
		int butterflySize = butterflyOffset >> 1;

		double alpha = -2.0*PI / butterflyOffset;
		double wR = cos(alpha), wI = sin(alpha);

		for (int i = 0; i < size / butterflyOffset; i++)
		{
			double uR = 1.0, uI = 0.0;
			double uRtmp;

			int offset = i * butterflyOffset;

			for (int j = 0; j < butterflySize; j++)
			{
				
				Butterfly(_mas, uR, uI, offset, butterflySize, step, 0,j);
				uRtmp = uR;
				uR = uR * wR - uI * wI;
				uI = uI * wR + uRtmp * wI;
			}
		}
	}
}
	

	

void SerialFFT1D(double* inputSignal, double*outputSignal, int size, int step, int offset)
{
	BitReversing(inputSignal, outputSignal, size, step, offset);
	SerialFFTCalculation(outputSignal, 0, size, step, offset);
}

void SerialFFT2D(double* inputSignal, double* outputSignal, int w, int h)
{
	double* tem = new double[2*w*h];
	std::thread thr1(tread_inverse, inputSignal, tem, w, h, 1, 0, w / 4, false);
	std::thread thr2(tread_inverse, inputSignal, tem, w, h, 1, w / 4, w / 2, false);
	std::thread thr3(tread_inverse, inputSignal, tem, w, h, 1, w / 2, 3 * w / 4, false);
	std::thread thr4(tread_inverse, inputSignal, tem, w, h, 1, 3 * w / 4, w, false);
	thr1.join();
	thr2.join();
	thr3.join();
	thr4.join();
	/*for (int i = 0; i < w; i++)
		SerialInverseFFT1D(inputSignal, tem, h, w, i);*/
	std::thread thr5(tread_inverse, tem, outputSignal, 1, w, h, 0, h / 4, false);
	std::thread thr6(tread_inverse, tem, outputSignal, 1, w, h, h / 4, h / 2, false);
	std::thread thr7(tread_inverse, tem, outputSignal, 1, w, h, h / 2, 3 * h / 4, false);
	std::thread thr8(tread_inverse, tem, outputSignal, 1, w, h, 3 * h / 4, h, false);
	thr5.join();
	thr6.join();
	thr7.join();
	thr8.join();
	/*for (int i = 0; i < w; i++)
		SerialFFT1D(inputSignal, tem, h, w, i);


	for (int j = 0; j < h; j++)
		SerialFFT1D(tem, outputSignal, w, 1, h*j);*/
	delete[] tem;
}
void SerialInverseFFTCalculation(double*mas, int size, int step, int off)
{
	int m = 0;
	for (int tmp_size = size; tmp_size > 1; tmp_size /= 2, m++);//size = 2^m

	double*_mas = mas + 2*off;

	for (int p = 0; p < m; p++)
	{
		int butterflyOffset = 1 << (p + 1);
		int butterflySize = butterflyOffset >> 1;

		double alpha = 2.0*PI / butterflyOffset;
		double wR = cos(alpha), wI = sin(alpha);

		for (int i = 0; i < size / butterflyOffset; i++)
		{
			double uR = 1.0, uI = 0.0;
			double uRtmp;

			int offset = i * butterflyOffset;

			for (int j = 0; j < butterflySize; j++)
			{
				Butterfly(_mas, uR, uI, offset, butterflySize, step, 0, j);
				uRtmp = uR;
				uR = uR * wR - uI * wI;
				uI = uI * wR + uRtmp * wI;
			}
		}
	}
}



void SerialInverseFFT1D(double* inputSignal, double* outputSignal, int size, int step, int offset)
{
	BitReversing(inputSignal, outputSignal, size, step, offset);
	SerialInverseFFTCalculation(outputSignal, size, step, offset);

	for (int j = 0; j < size; j++) {
		outputSignal[2 * (offset + step * j)] /= size;
		outputSignal[2 * (offset + step * j)+1] /= size;
	}
		
}


void SerialInverseFFT2D(double* inputSignal, double* outputSignal, int w, int h)
{
	double* tem = new double[2 * w*h];
	std::thread thr1(tread_inverse, inputSignal, tem, w, h, 1, 0, w/4, true);
	std::thread thr2(tread_inverse, inputSignal, tem, w, h, 1, w / 4, w / 2, true);
	std::thread thr3(tread_inverse, inputSignal, tem, w, h, 1, w / 2, 3*w / 4, true);
	std::thread thr4(tread_inverse, inputSignal, tem, w, h, 1, 3 * w / 4, w, true);
	thr1.join();
	thr2.join();
	thr3.join();
	thr4.join();
	/*for (int i = 0; i < w; i++)
		SerialInverseFFT1D(inputSignal, tem, h, w, i);*/
	std::thread thr5(tread_inverse, tem, outputSignal, 1, w, h, 0, h / 4, true);
	std::thread thr6(tread_inverse, tem, outputSignal, 1, w, h, h / 4, h / 2, true);
	std::thread thr7(tread_inverse, tem, outputSignal, 1, w, h, h / 2, 3 * h / 4, true);
	std::thread thr8(tread_inverse, tem, outputSignal, 1, w, h, 3 * h / 4, h, true);
	thr5.join();
	thr6.join();
	thr7.join();
	thr8.join();
	/*for (int j = 0; j < h; j++)
		SerialInverseFFT1D(tem, outputSignal, w, 1, h*j);*/
	delete[] tem;
}