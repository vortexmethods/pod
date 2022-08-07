/*--------------------------------*- POD -*------------------*---------------*\
|    #####    ####   #####      |                            | Version 1.0    |
|    ##  ##  ##  ##  ##  ##     |  POD: Proper Orthogonal    | 2022/08/01     |
|    #####   ##  ##  ##  ##     |  Decomposition method      *----------------*
|    ##      ##  ##  ##  ##     |  Open Source Code                           |
|    ##       ####   #####      |  https://www.github.com/vortexmethods/pod   |
|                                                                             |
| Copyright (C) 2020-2022 Ilia Marchevsky, Irina Soldatova, Kseniia Sokol     |
*-----------------------------------------------------------------------------*
| File name: CMakeLists.txt                                                   |
| Info: Source code of POD                                                    |
|                                                                             |
| This file is part of POD.                                                   |
| POD is free software: you can redistribute it and/or modify it              |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| POD is distributed in the hope that it will be useful, but WITHOUT          |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with POD.  If not, see <http://www.gnu.org/licenses/>.                |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief POD method for vtk-files parsing from VM2D
\author Марчевский Илья Константинович
\author Солдатова Ирина Александровна
\author Сокол Ксения Сергеевна
\version 1.0
\date 01 августа 2022 г.
*/


#include <bitset>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>


#include "Logo.h"
#include "Params.h"

#include "omp.h"

//#define EIGEN_USE_MKL_ALL
#include "../include/Eigen/Core"
#include "../include/Eigen/Dense"
using namespace Eigen;

typedef double doubleType;
typedef VectorXd eigDoubleType;
typedef MatrixXd matrEigDoubleType;
const char doubleTypeName[] = "double";

typedef float doubleTypeToRead;
typedef VectorXf eigDoubleTypeToRead;
typedef MatrixXf matrEigDoubleTypeToRead;

typedef int intType;

const char eolnBIN[] = "\n";



using namespace POD;

int main()
{
	PrintLogoToStream(std::cout);

	//Принудительное задание числа нитей в OpenMP (при необходимости)
	//int np = 1;
	//std::cout << "np = " << np << std::endl;
	//omp_set_num_threads(np);
	
	//Проверка, как хранятся числа в памяти
	uint16_t x = 0x0001;
	littleEndian = (*((uint8_t*)&x));

	PrintProperties();

	//Считывание параметров задачи
	std::ifstream file("./info.txt");
	if (!file)
	{
		std::cout << "ERROR: File \"info.txt\" not open\n";
		return -1;
	}		
	file >> TSIZE >> NSIZE >> startStep >> dStep >> dirData >> fileNameLength >> threshold;
	PrintConfiguration();
	

	//Считывание эталонного файла и координат узлов сетки
	std::ifstream eFile(dirData + fileNameStep("VelPres", fileNameLength, startStep, "vtk"), std::ios::in | std::ios::binary);
	
	//Считывание числа узлов в сетке в MSIZE
	std::string buff;
	char buff2[255];
	std::stringstream strStream;

	for (int q = 0; q < 5; ++q)
		getline(eFile, buff);
	strStream << buff;
	strStream >> buff2;
	strStream >> MSIZE;

	std::cout << "GridSize = " << MSIZE << std::endl << std::endl;

	//Заготовка памяти для исходных данных, их же в другом типе данных и матрицы ковариаций
	matrEigDoubleTypeToRead pDataRead = matrEigDoubleTypeToRead::Zero(MSIZE, TSIZE);
	matrEigDoubleType pData = matrEigDoubleType::Zero(MSIZE, TSIZE);
	matrEigDoubleType pCov = matrEigDoubleType::Zero(TSIZE, TSIZE);

	matrEigDoubleTypeToRead vDataRead = matrEigDoubleTypeToRead::Zero(MSIZE * 3, TSIZE);
	matrEigDoubleType vData = matrEigDoubleType::Zero(MSIZE * 3, TSIZE);
	matrEigDoubleType vCov = matrEigDoubleType::Zero(TSIZE, TSIZE);

	//Заготовка памяти под координаты узлов и их же в другом типе данных
	std::vector<doubleTypeToRead> rDataRead(3 * MSIZE, 0);
	std::vector<doubleType> rData(3 * MSIZE, 0);

	//Коды ячеек
	std::vector<intType> cells(2 * MSIZE);
	std::vector<intType> cellsTypes(MSIZE);

	//Считываем координаты узлов сетки
	eFile.read(reinterpret_cast<char*>(rDataRead.data()), MSIZE * 3 * sizeof(doubleTypeToRead));

	//Превращаем их в little-endian
	if (littleEndian)
		for (int i = 0; i < 3 * MSIZE; ++i)
			SwapEnd(rDataRead[i]);
	
	//Приводим к новому типу
	for (size_t q=0; q < rDataRead.size(); ++q)
		rData[q] = rDataRead[q];

	//Превращаем обратно в big-endian
	if (littleEndian)
		for (int i = 0; i < 3 * MSIZE; ++i)
		{
			SwapEnd(rData[i]);
			SwapEnd(rDataRead[i]);
		}
	

	//Считали конец строки и стрчку со словом CELLS
	eFile.get();
	getline(eFile, buff);

	//Номера ячеек, предваренные единицами
	eFile.read(reinterpret_cast<char*>(cells.data()), 2 * MSIZE * sizeof(int));

	//Считали конец строки и стрчку со словом CELLS_TYPES
	eFile.get();
	getline(eFile, buff);

	//Типы ячеек
	eFile.read(reinterpret_cast<char*>(cellsTypes.data()), MSIZE * sizeof(int));

	eFile.close();


	double startReadFiles = omp_get_wtime();
	// Чтение файлов с данными
	std::cout << "Reading input files... " << std::flush;
	std::string timeFile;

#pragma  omp parallel for private(timeFile) shared(vData,pData,vDataRead,pDataRead)
	for (int k = 0; k < TSIZE; k++)
	{
		timeFile = dirData + fileNameStep("VelPres", fileNameLength, startStep + dStep * k, "vtk");

		std::ifstream file(timeFile, std::ios::in | std::ios::binary);

		if (!file.is_open())
		{
			std::cout << "ERROR: Unable to open file " << timeFile << std::endl;
			exit(-1);
		}

		std::string buff;

		for (int q = 0; q < 5; ++q)
			getline(file, buff);

		//Смещаем указатель чтения, пропуская сетку и сведения о CELLS
		file.seekg(MSIZE * 3 * sizeof(doubleTypeToRead), std::ios::cur);
		file.get();
		getline(file, buff);

		file.seekg(2 * MSIZE * sizeof(int), std::ios::cur);
		file.get();
		getline(file, buff);

		file.seekg(MSIZE * sizeof(int), std::ios::cur);
		file.get();
		getline(file, buff);
		getline(file, buff);

		//Загружаем поле скоростей в узлах
		file.read(reinterpret_cast<char*>(vDataRead.col(k).data()), MSIZE * 3 * sizeof(doubleTypeToRead));
		file.get();

		if (littleEndian)
			for (int i = 0; i < MSIZE * 3; ++i)
				SwapEnd(vDataRead(i, k));

		//меняем тип	
		for (int i = 0; i < vDataRead.col(k).size(); ++i)
			vData.col(k)(i) = vDataRead.col(k)(i);

		getline(file, buff);
		getline(file, buff);

		//Загружаем поле давления в узлах
		file.read(reinterpret_cast<char*>(pDataRead.col(k).data()), MSIZE * sizeof(doubleTypeToRead));
		file.get();

		if (littleEndian)
			for (int i = 0; i < MSIZE; ++i)
				SwapEnd(pDataRead(i, k));

		//меняем тип
		for (int i = 0; i < pDataRead.col(k).size(); ++i)
			pData.col(k)(i) = pDataRead.col(k)(i);
		
		file.close();
	}//for k


	//Установка полей V и P в ноль правее прямой x = threshold	

	if (littleEndian)
		for (int i = 0; i < 3 * MSIZE; ++i)
			SwapEnd(rData[i]);

	for (int i = 0; i < MSIZE; i++)
	{
//		if (rData[3 * i] * rData[3 * i] + rData[3 * i + 1] * rData[3 * i + 1] < 0.25)
		if (rData[3 * i + 0] > threshold)
			for (int t = 0; t < TSIZE; t++)
			{
				for (int j = 0; j < 3; j++)
					vData(3*i + j, t) = 0.0;

				pData(i, t) = 0.0;
			}
	}

	if (littleEndian)
		for (int i = 0; i < 3 * MSIZE; ++i)
			SwapEnd(rData[i]);

	double endReadFiles = omp_get_wtime();
	double timeReadFiles = (endReadFiles - startReadFiles);
	std::cout << "Done" << std::endl;
	std::cout << "Time reading input files = " << timeReadFiles << std::endl << std::endl;

		
	//Расчет ковариационной матрицы

	std::cout << "Computing covariation matrices... " << std::flush;
	double startCalcMatrix = omp_get_wtime();
	
#pragma omp parallel for
	for (int i = 0; i < TSIZE; ++i)
		for (int j = 0; j < TSIZE; ++j)
		{
			pCov(i, j) = (1.0 / TSIZE) * (pData.col(i).dot(pData.col(j).transpose()));
			vCov(i, j) = (1.0 / TSIZE) * (vData.col(i).dot(vData.col(j).transpose()));
		}
	

	std::cout << "Done" << std::endl;
	double endCalcMatrix = omp_get_wtime();
	double timeCalcMatrix = (endCalcMatrix - startCalcMatrix);
	std::cout << "Time computing covariation matrices = " << timeCalcMatrix << std::endl << std::endl;


	//Собственные числа и собственные векторы

	double startEigen = omp_get_wtime();
	std::cout << "Computing eigenvalues/vectors... " << std::flush;
	SelfAdjointEigenSolver<matrEigDoubleType>
		eigensolverP(pCov);

	SelfAdjointEigenSolver<matrEigDoubleType>
		eigensolverV(vCov);

	if ((eigensolverP.info() != Success) || (eigensolverV.info() != Success))
	{
		std::cout << "ERROR: calculating eigenvalues" << std::endl;
		exit(-1);
	}

	eigDoubleType eigvalP = eigensolverP.eigenvalues().reverse();
	matrEigDoubleType eigvecP = eigensolverP.eigenvectors().rowwise().reverse();

	eigDoubleType eigvalV = eigensolverV.eigenvalues().reverse();
	matrEigDoubleType eigvecV = eigensolverV.eigenvectors().rowwise().reverse();

	std::cout << "Done" << std::endl;
	double endEigen = omp_get_wtime();
	double timeEigen = (endEigen - startEigen);
	std::cout << "Time computing eigenvectors = " << timeEigen << std::endl << std::endl;

	
	// Вычисление POD-мод
	double startModes = omp_get_wtime();
	std::cout << "Computing POD modes... " << std::flush;
	matrEigDoubleType podP = matrEigDoubleType::Zero(MSIZE, TSIZE);
	matrEigDoubleType podV = matrEigDoubleType::Zero(MSIZE * 3, TSIZE);
	
	double cft = 0.0;
#pragma omp parallel for private (cft)
	for (int i = 0; i < TSIZE; ++i)
	{
		cft = 1.0 / sqrt(eigvalP(i) * TSIZE);
		for (int j = 0; j < TSIZE; ++j)
			podP.col(i) += cft * eigvecP(j, i) * pData.block(0, j, MSIZE, 1);
	}

#pragma  omp parallel for private (cft)
	for (int i = 0; i < TSIZE; ++i)
	{
		cft = 1.0 / sqrt(eigvalV(i) * TSIZE);
		for (int j = 0; j < TSIZE; j++)
			podV.col(i) += cft * eigvecV(j, i) * vData.block(0, j, MSIZE * 3, 1);
	}

	std::cout << "Done" << std::endl;
	double endModes = omp_get_wtime();
	double timeModes = (endModes - startModes);
	std::cout << "Time computing POD modes = " << timeModes << std::endl << std::endl;


	CreateDirectory(dirData, "output");

	//Сохранение собственных значений
	std::cout << "Writing eigenvalues... " << std::flush;
	std::ofstream writeEigvalP(dirData + "/output/eigenP.txt");
	if (writeEigvalP.is_open())
	{
		writeEigvalP << std::scientific << std::setprecision(16) << eigvalP;
		writeEigvalP.close();
	}

	std::ofstream writeEigvalV(dirData + "/output/eigenV.txt");
	if (writeEigvalV.is_open())
	{
		writeEigvalV << std::scientific << std::setprecision(16) << eigvalV;
		writeEigvalV.close();
	}

	std::cout << "Done" << std::endl << std::endl;

	//Запись временных коэффициентов
	std::cout << "Writing time coefficients... " << std::flush;
	//double startTimeCoeff = omp_get_wtime();
	std::vector<std::vector<double>> timeCftP(NSIZE), timeCftV(NSIZE);
		
#pragma omp parallel for
	for (int i = 0; i < NSIZE; ++i)
	{
		timeCftP[i].resize(TSIZE);
		timeCftV[i].resize(TSIZE);
			
		std::ofstream writeChronos(fileNameStep(dirData + "/output/timecft", 3, i, "dat"));
		for (int j = 0; j < TSIZE; j++)
		{
			timeCftP[i][j] = pData.col(j).dot(podP.col(i));
			timeCftV[i][j] = vData.col(j).dot(podV.col(i));
			writeChronos << std::scientific << std::setprecision(16) << timeCftP[i][j] << " " << timeCftV[i][j] << '\n';
		}
		writeChronos.close();
	}
	std::cout << "Done" << std::endl << std::endl;


	//Запись мод
	std::cout << "Writing POD modes... " << std::flush;
	//double startWritingModes = omp_get_wtime();

#pragma omp parallel for
	for (int k = 0; k < NSIZE; k++)
	{
		std::ofstream outfile(fileNameStep(dirData + "/output/mode", 3, k, "vtk"), std::ios::out | std::ios::binary);	

		outfile << "# vtk DataFile Version 3.0" << eolnBIN << "POD result: "  << eolnBIN;
		outfile << "BINARY" << eolnBIN;
		outfile << "DATASET UNSTRUCTURED_GRID" << eolnBIN << "POINTS " << MSIZE << " " << doubleTypeName << eolnBIN;
		
		outfile.write(reinterpret_cast<char*>(rData.data()), MSIZE * 3 * sizeof(doubleType));

		// CELLS
		outfile << eolnBIN << "CELLS " << MSIZE << " " << MSIZE * 2 << eolnBIN;
		outfile.write(reinterpret_cast<char*>(cells.data()), MSIZE * 2 * sizeof(int));
		outfile << eolnBIN << "CELL_TYPES " << MSIZE << eolnBIN;
		outfile.write(reinterpret_cast<char*>(cellsTypes.data()), MSIZE * sizeof(int));

		//VECTORS V
		outfile << eolnBIN << "POINT_DATA " << MSIZE << eolnBIN;
		outfile << "VECTORS V " << doubleTypeName << eolnBIN;
		if (littleEndian)
			for (int i = 0; i < MSIZE * 3; ++i)
				SwapEnd(podV(i, k));
		
		outfile.write(reinterpret_cast<char*>(podV.col(k).data()), MSIZE * 3 * sizeof(doubleType));
		
		if (littleEndian)
			for (int i = 0; i < MSIZE * 3; ++i)
				SwapEnd(podV(i, k));

		//SCALARS P	
		outfile << eolnBIN << "SCALARS P " << doubleTypeName << " 1" << eolnBIN;
		outfile << "LOOKUP_TABLE default" << eolnBIN;
		if (littleEndian)
			for (int i = 0; i < MSIZE; ++i)
				SwapEnd(podP(i, k));
		
		outfile.write(reinterpret_cast<char*>(podP.col(k).data()), MSIZE * sizeof(doubleType));
		
		if (littleEndian)
			for (int i = 0; i < MSIZE; ++i)
				SwapEnd(podP(i, k));

		outfile << eolnBIN;

		outfile.close();
	}
	std::cout << "Done" << std::endl << std::endl;



	std::cout << "Writing approxSol... " << std::flush;
	double startApproxSol = omp_get_wtime();
		
#pragma omp parallel for
	for (int t = 0; t < TSIZE; ++t)
	{
		std::vector<doubleType> TmpSolV(3 * MSIZE, 0.);
		std::vector<doubleType> TmpSolP(MSIZE, 0.);

		double tmpx = 0.0, tmpy = 0.0, tmpp = 0.0;
		for (int q = 0; q < MSIZE; ++q)
		{
			tmpx = tmpy = tmpp = 0.0;
			for (int md = 0; md < NSIZE; ++md)
			{
				tmpx += timeCftV[md][t] * podV.col(md)(3 * q + 0);
				tmpy += timeCftV[md][t] * podV.col(md)(3 * q + 1);
				tmpp += timeCftP[md][t] * podP.col(md)(q);
			}

			TmpSolV[3 * q + 0] = tmpx;
			TmpSolV[3 * q + 1] = tmpy;
			TmpSolP[q] = tmpp;
		}

		if (littleEndian)
			for (int i = 0; i < MSIZE * 3; ++i)
				SwapEnd(TmpSolV[i]);

		if (littleEndian)
			for (int i = 0; i < MSIZE; ++i)
				SwapEnd(TmpSolP[i]);
			   
		std::ofstream outfile(fileNameStep(dirData + "/output/Sol", 5, startStep + dStep * t, "vtk"), std::ios::out | std::ios::binary);

		outfile << "# vtk DataFile Version 3.0" << eolnBIN << "POD result: " << eolnBIN;
		outfile << "BINARY" << eolnBIN;
		outfile << "DATASET UNSTRUCTURED_GRID" << eolnBIN << "POINTS " << MSIZE << " " << doubleTypeName << eolnBIN;

		outfile.write(reinterpret_cast<char*>(rData.data()), MSIZE * 3 * sizeof(doubleType));

		// CELLS
		outfile << eolnBIN << "CELLS " << MSIZE << " " << MSIZE * 2 << eolnBIN;
		outfile.write(reinterpret_cast<char*>(cells.data()), MSIZE * 2 * sizeof(int));
		outfile << eolnBIN << "CELL_TYPES " << MSIZE << eolnBIN;
		outfile.write(reinterpret_cast<char*>(cellsTypes.data()), MSIZE * sizeof(int));

		//VECTORS V
		outfile << eolnBIN << "POINT_DATA " << MSIZE << eolnBIN;
		outfile << "VECTORS V " << doubleTypeName << eolnBIN;
		
		outfile.write(reinterpret_cast<char*>(TmpSolV.data()), MSIZE * 3 * sizeof(doubleType));
	
		//SCALARS P	
		outfile << eolnBIN << "SCALARS P " << doubleTypeName << " 1" << eolnBIN;
		outfile << "LOOKUP_TABLE default" << eolnBIN;
		
		outfile.write(reinterpret_cast<char*>(TmpSolP.data()), MSIZE * sizeof(doubleType));

		outfile << eolnBIN;

		outfile.close();

		if (littleEndian)
			for (int i = 0; i < MSIZE * 3; ++i)
				SwapEnd(TmpSolV[i]);

		if (littleEndian)
			for (int i = 0; i < MSIZE; ++i)
				SwapEnd(TmpSolP[i]);
	}


	double endApproxSol = omp_get_wtime();
	double timeApproxSol = endApproxSol - startApproxSol;
	std::cout << " Done" << std::endl;
	std::cout << "Time Approximate Solution = " << timeApproxSol << std::endl << std::endl;
	
	/*	
	"Для разминки"
	1. Подсчитать разность между исходными полями и полученными с использованием удерживаемого количества мод
	2. Сравнить полученный результат с оценкой через собственные числа
	3. Проанализировать формы 1 и 2 (считая, что они нумеруются с 0-й!) для установившегося следа за цилиндром, взяв их комбинации с периодическими множителями cos(t) и sin(t) - можно в WM, можно в С++
	4. То же с формами 3 и 4
	5. Понять, не получится ли результат в целом (по погрешности) лучше, если считать Vx и Vy отдельными полями (модифицировав код)

	"По-серьезному"
	1. Познакомиться с понятием "информационная энтропия" https://ru.wikipedia.org/wiki/Информационная_энтропия 
	2. Для последовательных кадров velPres посчитать разность полей в узлах сетки (из следующего каждый раз вычитаем предыдущий).
	3. Полученную разность отнормировать, разделив на максимальную (по всей задаче) вариацию скорости и давления соответственно, так что результат будет выражен числом из диапазона (-1, 1)
	4. Эти дробные числа заменить рациональными, считая что рациональное число кодируется целым числом битов: т.е. числами вида x/16, x/32, x/64... 
	   Разумеется, хранить надо только числители этих чисел, т.е. целые числа.
	   Здравый смысл подсказывает, что знаменатель должен быть довольно большим, не менее 2^11 = 2048.
	   Последнее будет означать, что погрешность будет составлять максимум 0.05 %, а для хранения числа будет достаточно 11 битов (для сравнения - float = 32 бита, double = 64 бита).
	   Скорее всего, знаменатель должен быть даже заметно больше по следующей причине: 
	   добавки на каждом шаге будут как правило малыми, поэтому большие числа почти не будут встречаться, а малые нужно различать как следует.	   
	5. Оценить величину ошибки после обратного восстановления кадров, отсюда найти разумное число бит для кодирования рациональных чисел.
	6. Производим расчет информационной энтропии каждого такого "прибавочно-отнормированного кадра".
	7. Количество информации (в битах) в каждом кадре будет равно полученной энтропии, умноженной на количество чисел в кадре (т.е. число узлов сетки для скалярного поля)
	8. Просуммировать количество информации по всем кадрам и оценить "информационную содержательность" результата расчета
	*/
}