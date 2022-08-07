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
\brief Параметры решаемой задачи
\author Марчевский Илья Константинович
\author Солдатова Ирина Александровна
\author Сокол Ксения Сергеевна
\version 1.0
\date 01 августа 2022 г.
*/

#ifndef PARAMS_H_
#define PARAMS_H_

#if defined(_WIN32)
#include <direct.h>
#endif

#include <iostream>

namespace POD
{
	// Способ хранения многобайтных чисел
	static bool littleEndian;
	
	// Размерность матрицы ковариаций (общее число точек в сетке, на которой вычислены поля)
	static int MSIZE;

	// Число слоев по времени (количество обрабатываемых кадров)
	static int TSIZE;

	// Количество сохраняемых мод
	static int NSIZE;

	// Номер кадра, с которого начинается обработка
	static int startStep;

	// Шаг номеров кадров
	static int dStep;
	
	// Отсечка правой границы (все поля справа от линии x = threshold обнуляются)
	static double threshold;
	
	// Папка, где лежат данные в виде vtk-файлов
	static std::string dirData;

	// Количество разрядов в имени файла
	static int fileNameLength;

	// Файл с исходными параметрами
	static const std::string infofileName = "./info.txt";




	/// \brief Формирование имени файла
	///
	/// \param[in] name константная ссылка на префикс имени файла
	/// \param[in] length количество знаков под номер
	/// \param[in] number номер для имени файла
	/// \param[in] ext константная ссылка на раширение (дописывается, если непустое)
	/// \return строку --- имя текстового файла
	inline std::string fileNameStep(const std::string& name, int length, size_t number, const std::string& ext)
	{
		std::string fname(name);

		size_t dec = 1;

		for (int i = 1; i < length; ++i)
		{
			dec *= 10;
			if (number < dec)
				fname += "0";
		}

		std::ostringstream ss;
		ss << number;
		fname += ss.str();

		if (ext.size() > 0)
		{
			fname += ".";
			fname += ext;
		}

		return fname;
	}


	template<typename T>
	void SwapEnd(T& var)
	{
		char* varArray = reinterpret_cast<char*>(&var);
		for (long i = 0; i < static_cast<long>(sizeof(var) / 2); ++i)
			std::swap(varArray[sizeof(var) - 1 - i], varArray[i]);
	}
	

	/// \brief Создание каталога
	///
	/// \param[in] dir константная ссылка на имя текущего каталога (со слешем на конце)
	/// \param[in] name константная ссылка на имя создаваемого подкаталога
	inline void CreateDirectory(const std::string& dir, const std::string& name)
	{
#if defined(_WIN32)
		_mkdir((dir + name).c_str());
#else
		mkdir((dir + name).c_str(), S_IRWXU | S_IRGRP | S_IROTH);
#endif
	}


}//namespace POD

#endif