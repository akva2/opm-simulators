/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#ifndef HDF5_FILE_HPP
#define HDF5_FILE_HPP

#include <hdf5.h>

#include <string>
#include <vector>

namespace Opm {

class HDF5File {
public:
    enum class OpenMode {
        APPEND,
        READ,
        WRITE
    };

    HDF5File(const std::string& fileName, OpenMode mode);

    ~HDF5File();

    bool write(const std::string& group,
               const std::string& dset,
               const std::vector<char>& buffer);

    bool read(const std::string& group,
              const std::string& dst,
              std::vector<char>& buffer);

    int lastGroup() const;

private:
    bool groupExists(hid_t parent, const std::string& path);

    hid_t m_file;
};

}

#endif
