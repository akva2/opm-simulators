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
#ifndef ECL_HDF5_SERIALIZER_HH
#define ECL_HDF5_SERIALIZER_HH

#include <opm/common/utility/Serializer.hpp>

#include <opm/simulators/utils/HDF5File.hpp>
#include <opm/simulators/utils/HDF5Packer.hpp>

namespace Opm {

template<>
template<class B, class A>
struct Serializer<Serialization::Packer>::is_vector<Dune::BlockVector<B,A>> {
    constexpr static bool value = true;
};

//! \brief Class for (de-)serializing using HDF5.
class HDF5Serializer : public Serializer<Serialization::Packer> {
public:
    HDF5Serializer(const std::string& fileName, HDF5File::OpenMode mode)
        : Serializer<Serialization::Packer>(m_packer)
        , m_packer{}
        , m_h5file(fileName, mode)
    {}

    //! \brief Serialize and write data to restart file.
    //!
    //! \tparam T Type of class to write
    //! \param data Simulator to write restart data for
    template<class T>
    void write(T& data,
               const std::string& group,
               const std::string& dset)
    {
        try {
            this->pack(data);
        } catch (...) {
            m_packSize = std::numeric_limits<size_t>::max();
            throw;
        }

        m_h5file.write(group, dset, m_buffer);
    }

    //! \brief read data and deserialize from restart file.
    //!
    //! \tparam T Type of class to read
    //! \param data Simulator to read restart data for
    template<class T>
    void read(T& data,
              const std::string& group,
              const std::string& dset)
    {
        m_h5file.read(group, dset, m_buffer);
        this->unpack(data);
    }

protected:
private:
    const Serialization::Packer m_packer; //!< Packer instance
    HDF5File m_h5file;
};

}

#endif
