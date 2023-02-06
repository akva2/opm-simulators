/*
  Copyright 2020 Equinor.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>

#include <ebos/hdf5serializer.hh>

#include <opm/simulators/timestepping/SimulatorTimer.hpp>

#include <boost/date_time.hpp>

#include <fmt/format.h>

#include <array>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Need one parameter, the .SAVE file to inspect\n";
        return 1;
    }

    Opm::HDF5Serializer ser(argv[1], Opm::HDF5File::OpenMode::READ);

    std::array<std::string,4> strings;
    try {
        ser.read(strings, "/", "simulator_info");
    } catch(...) {
        std::cerr << "Error reading data from file, is it really a .SAVE file?\n";
        return 2;
    }

    std::cout << "Info for " << argv[1] <<":\n";
    std::cout << fmt::format("\tSimulator name: {}\n"
                             "\tSimulator version: {}\n"
                             "\tCompile time stamp: {}\n"
                             "\tCase name: {}\n",
                             strings[0], strings[1], strings[2], strings[3]);

    int lastGroup = ser.lastGroup();
    std::cout << fmt::format("\tNumber of report steps: {}\n", lastGroup);
    for (int i = 1; i <= lastGroup; ++i) {
        Opm::SimulatorTimer timer;
        try {
            ser.read(timer, "/" + std::to_string(i), "simulator_timer");
        } catch (...) {
            std::cerr << "*** Failed to read timer info for level " << i << std::endl;
        }
        std::cout << "\t\tReport step id " << i << ": Time "
                  << timer.currentDateTime() << std::endl;
    }

    return 0;
}
