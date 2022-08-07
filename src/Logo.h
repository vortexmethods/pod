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
\brief Логотип
\author Марчевский Илья Константинович
\author Солдатова Ирина Александровна
\author Сокол Ксения Сергеевна
\version 1.0
\date 01 августа 2022 г.
*/

#ifndef LOGO_H_
#define LOGO_H_

#include <fstream>

#include "omp.h"

#include "Params.h"

namespace POD
{

	void PrintLogoToStream(std::ostream& str)
	{
		str <<
			"/*--------------------------------*- POD -*------------------*---------------*\\" << '\n' << \
			"|    #####    ####   #####      |                            | Version 1.0    |" << '\n' << \
			"|    ##  ##  ##  ##  ##  ##     |  POD: Proper Orthogonal    | 2022/08/01     |" << '\n' << \
			"|    #####   ##  ##  ##  ##     |  Decomposition method      *----------------*" << '\n' << \
			"|    ##      ##  ##  ##  ##     |  Open Source Code                           |" << '\n' << \
			"|    ##       ####   #####      |  https://www.github.com/vortexmethods/pod   |" << '\n' << \
			"|                                                                             |" << '\n' << \
			"| Copyright (C) 2020-2022 Ilia Marchevsky, Irina Soldatova, Kseniia Sokol     |" << '\n' << \
			"\\*---------------------------------------------------------------------------*/" << '\n';
	}

	void PrintProperties()
	{
		printf("\n");
		printf("                             CPU Properties                             \n");
		printf("------------------------------------------------------------------------\n");
		printf("Multi-bytes numbers storing:  %s\n", (littleEndian ? "little-endian" : "big-endian   "));
		printf("OpenMP max threads:                 %d\n", omp_get_max_threads());
		//printf("OpenMP nested parallelism:        %s\n", (omp_get_nested() ? "on" : "off"));		
		//printf("Nested OpenMP up to tree level:     %d\n", maxLevelOmp);
		printf("------------------------------------------------------------------------\n");
	}

	void PrintConfiguration()
	{
		std::cout << std::endl;
		std::cout << "                     Configuration for task                             " << std::endl;
		std::cout << "------------------------------------------------------------------------" << std::endl;
		//std::cout << "Folder:      " << dirData.c_str() << std::endl;
		std::cout << "FirstFile: " << dirData + fileNameStep("VelPres", fileNameLength, startStep, "vtk") << std::endl;
		std::cout << "LastFile:  " << dirData + fileNameStep("VelPres", fileNameLength, startStep + (TSIZE-1) * dStep, "vtk") << std::endl;
		std::cout << "           (number of snapshots = " << TSIZE << ", step = " << dStep << ")" << std::endl;
		//std::cout << "Start snapshot = " << startStep << ", step = " << dStep << ", number of snapshots = " << TSIZE << std::endl;
		//std::cout << "Covariation matrix size = " << MSIZE << ",  right threshold = " << threshold << std::endl;
		std::cout << "Number of POD modes to be stored = " << NSIZE << std::endl;
		std::cout << "------------------------------------------------------------------------" << std::endl;
	}

/*
	void PrintStatistics(int run, int runs,
		const double* timing, const double* mintiming, const double* avtiming,
		double runtime, double minruntime, double avruntime)
	{
		if (run == 0)
		{
			std::cout << std::endl;
			std::cout << "                         Time statistics                                " << std::endl;
			std::cout << "------------------------------------------------------------------------" << std::endl;
			std::cout << "  run/runs  total        TBK     SKK     FCK      ker.time              " << std::endl;
		}


		printf("# %3d/%-3d: %6.3lf s  (", run + 1, runs, runtime);

		printf(" %6.3f ", timing[1]);
		printf(" %6.3f ", timing[2]);
		printf(" %6.3f ", timing[3]);
		printf(") = %6.3f s\n", timing[4]);
		
		if (run == runs - 1)
		{
			std::cout << "------------------------------------------------------------------------" << std::endl;

			printf("     min : %6.3lf s  (", minruntime);

			printf(" %6.3f ", mintiming[1]);
			printf(" %6.3f ", mintiming[2]);
			printf(" %6.3f ", mintiming[3]);
			printf(") | %6.3f s\n", mintiming[4]);


			printf("    aver : %6.3lf s  (", avruntime / runs);

			printf(" %6.3f ", avtiming[1] / runs);
			printf(" %6.3f ", avtiming[2] / runs);
			printf(" %6.3f ", avtiming[3] / runs);
			printf(") | %6.3f s\n", avtiming[4] / runs);

#ifdef calcOp
			std::cout << "------------------------------------------------------------------------" << std::endl;
			std::cout << "Operations:   " << op << " operations                                   " << std::endl;
#endif

		}
	}


	void PrintAccuracyHead()
	{
		std::cout << std::endl;
		std::cout << "                         Accuracy control                               " << std::endl;
		std::cout << "------------------------------------------------------------------------" << std::endl;
	}

	void PrintAccuracyError(double val)
	{
		std::cout << "Relative error:     " << val << std::endl;

		std::cout << "------------------------------------------------------------------------" << std::endl;
	}
*/

}
#endif