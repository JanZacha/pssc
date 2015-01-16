#include "standardlibs.h"

mas_exception::mas_exception(const string& msg)
{
        snprintf(m_msg, sizeof(m_msg), "%s", msg.c_str());
}

mas_exception::mas_exception(const boost::format& msg)
{
        snprintf(m_msg, sizeof(m_msg), "%s", msg.str().c_str());
}