#pragma once

#include <iostream>
#include <memory>
#include <map>
#include <deque>
#include <queue>
#include <vector>
#include <unordered_map>
#include <functional>

#include <boost/program_options.hpp>
#include <boost/format.hpp>


#include <irrlicht/irrlicht.h>
#include <irrlicht/quaternion.h>
#include <irrlicht/vector3d.h>

#include <boost/algorithm/string.hpp>
#include <boost/cstdint.hpp>

namespace po = boost::program_options;
namespace ba = boost::algorithm;

typedef int8_t          int8;
typedef uint8_t         uint8;
typedef int16_t         int16;
typedef uint16_t        uint16;
typedef int32_t         int32;
typedef uint32_t        uint32;
typedef int64_t         int64;
typedef uint64_t        uint64;

using std::string;
using std::deque;
using std::unique_ptr;
using std::make_shared;
using std::cout;
using std::cerr;
using std::endl;

using std::ostream;
using std::istream;




typedef irr::core::quaternion quaternion;
typedef irr::core::vector3df vector3d;
typedef irr::core::vector3di vector3i;

class mas_exception : public std::exception
{
  public:
    mas_exception(const std::string& msg);
    mas_exception(const boost::format& msg);

    virtual const char*  what() const throw()    { return m_msg; }

  private:
    char m_msg[1024];
};
