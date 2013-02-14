/* sdsl - succinct data structures library
    Copyright (C) 2008 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
#include "sdsl/testutils.hpp"
#include <cxxabi.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <unistd.h> // for get_file_size, also contains clock_gettime
#include <sys/types.h> // for get_file_size
#include <sys/stat.h>  // for get_file_size

namespace sdsl
{

timeval stop_watch::m_first_t = {0,0};
rusage stop_watch::m_first_r = {{0,0},{0,0}};

void stop_watch::start()
{
    gettimeofday(&m_timeOfDay1, 0);
    getrusage(RUSAGE_SELF, &m_ruse1);
}

void stop_watch::stop()
{
    getrusage(RUSAGE_SELF, &m_ruse2);
    gettimeofday(&m_timeOfDay2, 0);
}

double stop_watch::get_user_time()
{
    timeval t1, t2;
    t1 = m_ruse1.ru_utime;
    t2 = m_ruse2.ru_utime;
    return ((double)(t2.tv_sec*1000000 + t2.tv_usec - (t1.tv_sec*1000000 + t1.tv_usec)))/1000.0;
}

double stop_watch::get_sys_time()
{
    timeval t1, t2;
    t1 = m_ruse1.ru_stime;
    t2 = m_ruse2.ru_stime;
    return ((double)(t2.tv_sec*1000000 + t2.tv_usec - (t1.tv_sec*1000000 + t1.tv_usec)))/1000.0;
}

double stop_watch::get_real_time()
{
    double result = ((double)((m_timeOfDay2.tv_sec*1000000 + m_timeOfDay2.tv_usec)-(m_timeOfDay1.tv_sec*1000000 + m_timeOfDay1.tv_usec)))/1000.0;
    if (result < get_sys_time() + get_user_time())
        return get_sys_time()+get_user_time();
    return result;
}

uint64_t stop_watch::get_abs_real_time()
{
    uint64_t result = (((m_timeOfDay2.tv_sec*1000000 + m_timeOfDay2.tv_usec - (m_first_t.tv_sec*1000000 + m_first_t.tv_usec))))/1000;
    return result;
}

uint64_t stop_watch::get_abs_user_time()
{
    timeval t1, t2;
    t1 = m_first_r.ru_utime;
    t2 = m_ruse2.ru_utime;
    return (t2.tv_sec*1000000 + t2.tv_usec - (t1.tv_sec*1000000 + t1.tv_usec))/1000;
}


uint64_t stop_watch::get_abs_sys_time()
{
    timeval t1, t2;
    t1 = m_first_r.ru_stime;
    t2 = m_ruse2.ru_stime;
    return (t2.tv_sec*1000000 + t2.tv_usec - (t1.tv_sec*1000000 + t1.tv_usec))/1000;
}

uint64_t stop_watch::get_abs_page_faults()
{
    return m_ruse2.ru_majflt - m_first_r.ru_majflt; // does not work on my platform
}

std::string clock::get_time_string()
{
    time_t rawtime;
    struct tm* timeinfo;
    char buffer[1024];
    time(&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer, 1024, "%Y-%m-%d-%H%M%S", timeinfo);
    return buffer;
}


pattern_file::size_type
pattern_file::generate(const char* text_file_name, size_type pattern_cnt, size_type pattern_len)
{
    char* text = NULL;
    size_t n = util::file::read_text(text_file_name,text);

    m_pattern_cnt = pattern_cnt;
    m_pattern_len = pattern_len;

    char* buf = new char[pattern_len+1];
    size_t j=0;
    m_patterns.resize(m_pattern_cnt);
    for (size_type i=0; j<pattern_cnt; i++) {
        size_type pos = rand() % (n-pattern_len);
        strncpy(buf,(const char*)(text+pos),pattern_len);
        buf[pattern_len-1] = 0;

        if (strlen(buf) == pattern_len) {
            m_patterns[i] = buf;
            j++;
        }
    }
    delete [] buf;
    delete [] text;
    return m_pattern_cnt;
}


pattern_file::size_type
pattern_file::serialize(std::ostream& out, structure_tree_node*, std::string) const
{
    out << "num_patterns=" << m_pattern_cnt << std::endl;
    out << "pattern_len=" << m_pattern_len << std::endl;

    for (size_t i=0; i<m_pattern_cnt; i++) {
        out << util::encode_base64(m_patterns[i]) << std::endl;
    }
    return m_pattern_cnt*m_pattern_len;
}

void
pattern_file::load(std::istream& in)
{
    std::string buf;
    /* parse args argument */
    std::getline(in,buf);
    sscanf(buf.c_str(),"num_patterns=%ld",&m_pattern_cnt);
    std::getline(in,buf);
    sscanf(buf.c_str(),"pattern_len=%ld",&m_pattern_len);

    m_patterns.resize(m_pattern_cnt);
    for (size_t i=0; i<m_pattern_cnt; i++) {
        std::getline(in,buf);
        m_patterns[i] = util::decode_base64(buf);
    }
}


} // end of namespace sdsl
