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
/*! \file testutils.hpp
 *  \brief testutils.hpp contains a "stopwatch" class for performance meassurement of program pieces.
 *  \author Simon Gog
 */
#ifndef INCLUDE_SDSL_TESTUTILS
#define INCLUDE_SDSL_TESTUTILS

#include "util.hpp"
#include "uintx_t.hpp"
#include <sys/time.h> // for struct timeval
#include <sys/resource.h> // for struct rusage
#include <iomanip>
#include <iostream>
#include <string>

#include "int_vector.hpp"

using namespace std;
using namespace sdsl;

namespace sdsl
{

//template<uint8_t, class size_type_class>
//class int_vector;    // forward declaration

//! A helper class to meassure the time consumption of program pieces.
/*! stop_watch is a stopwatch based on the commands getrusage and
 *  gettimeofday. Where getrusage is used to determine the user and system time
 *  and gettimeofday to determine the elapsed real time.
 */
class stop_watch
{
    private:
        rusage m_ruse1, m_ruse2;
        timeval m_timeOfDay1, m_timeOfDay2;
        static timeval m_first_t;
        static rusage m_first_r;
    public:

        stop_watch() : m_ruse1(), m_ruse2(), m_timeOfDay1(), m_timeOfDay2() {
            timeval t;
            t.tv_sec = 0; t.tv_usec = 0;
            m_ruse1.ru_utime = t; m_ruse1.ru_stime = t; // init m_ruse1
            m_ruse2.ru_utime = t; m_ruse2.ru_stime = t; // init m_ruse2
            m_timeOfDay1 = t; m_timeOfDay2 = t;
            if (m_first_t.tv_sec == 0) {
                gettimeofday(&m_first_t, 0);
            }
            if (m_first_r.ru_utime.tv_sec == 0 and m_first_r.ru_utime.tv_usec ==0) {
                getrusage(RUSAGE_SELF, &m_first_r);
            }
        }
        //! Start the stopwatch.
        /*! \sa stop
         */
        void start();

        //! Stop the stopwatch.
        /*! \sa start
         */
        void stop();

        //! Get the elapsed user time in milliseconds between start and stop.
        /*! \sa start, stop, get_real_time, get_sys_time
         */
        double get_user_time();

        //! Get the elapsed system time in milliseconds between start and stop.
        /*! \sa start, stop, get_real_time, get_user_time
         */
        double get_sys_time();

        //! Get the elapsed real time in milliseconds between start and stop.
        /*! \sa start, stop, get_sys_time, get_user_time
         */
        double get_real_time();

        //! Get the elapsed user time in milliseconds since the first construction of a stop_watch in the current process.
        /*! \sa get_user_time
         */
        uint64_t get_abs_user_time();

        //! Get the elapsed system time in milliseconds since the first construction of a stop_watch in the current process.
        /*! \sa get_sys_time
         */
        uint64_t get_abs_sys_time();

        //! Get the elapsed real time in milliseconds since the first construction of a stop_watch in the current process.
        /*! \sa get_real_time
         */
        uint64_t get_abs_real_time();

        uint64_t get_abs_page_faults();
};

//! Write stopwatch output in readable format
inline void write_R_output(std::string data_structure, std::string action,
                           std::string state="begin", uint64_t times=1, uint64_t check=0)
{
    if (util::verbose) {
        stop_watch _sw;
        _sw.stop();
        std::cout << data_structure << "\t" << action << "\t" << state << "\t"
                  << std::setw(9)<< times << "\t" << std::setw(9) << check << "\t"
                  << std::setw(9) << _sw.get_abs_real_time() << "\t "
                  << std::setw(9) << _sw.get_abs_user_time() << "\t"
                  << std::setw(9) << _sw.get_abs_sys_time() << std::endl;
    }
}

//! A helper class to get time information.
class clock
{
    public:
        static std::string get_time_string();
};


class pattern_file
{
    public:
        typedef bit_vector::size_type size_type;
    private:
        size_type m_pattern_cnt;
        size_type m_pattern_len;
        vector<std::string> m_patterns;
    public:
        pattern_file();
        size_type generate(const char* text_file_name, size_type pattern_cnt, size_type pattern_len);
        template<class tCst>
        void generate_restricted(const tCst& cst, const char* text_file_name, size_type pattern_cnt,
                                 size_type pattern_len, size_type min_occ, size_type max_occ);
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="") const;
        void load(std::istream& in);
        size_type num_patterns() {
            return m_pattern_cnt;
        };
        size_type pattern_length() {
            return m_pattern_len;
        };
        const char* operator[](size_type i) const;
};


template<class tCst>
void pattern_file::generate_restricted(const tCst& cst, const char* text_file_name, size_type pattern_cnt,
                                       size_type pattern_len, size_type min_occ, size_type max_occ)
{
    typedef typename tCst::node_type node_type;
    m_pattern_cnt = pattern_cnt;
    m_pattern_len = pattern_len;

    vector<size_type> candidates;
    for (typename tCst::const_iterator it = cst.begin(), end=cst.end(); it!=end; ++it) {
        if (it.visit() == 1) {
            node_type v = *it;
            if (cst.depth(v) >= pattern_len) {  // if the depth >= pattern_len we determine
                // if min_occ <= occ <= max_occ
                size_type occ = cst.leaves_in_the_subtree(v);
                if (min_occ <= occ and occ <= max_occ) {
                    candidates.push_back(cst.lb(v));
                }
                it.skip_subtree();
            } else { // if the depth < pattern_len we process the subtree
                // if the subtree size is to small we don't process the subtree
                if (min_occ > cst.leaves_in_the_subtree(v)) {
                    it.skip_subtree();
                }
            }
        }
    }
    sdsl::bit_vector already_used(candidates.size(), 0);

    if (candidates.size() == 0) {
        m_pattern_cnt = 0;
        std::cerr << "Found no candidates for pattern_len = "<< pattern_len << std::endl;
        std::cerr << "min_occ = "<< min_occ << " and max_occ = " << max_occ << std::endl;
    } else if (candidates.size() < 100) {
        std::cerr << "Less then 100 candidates for pattern_len = " << pattern_len << std::endl;
        std::cerr << "min_occ = "<< min_occ << " and max_occ = " << max_occ << std::endl;
    } else {
        std::cerr << candidates.size() << " candidates for pattern_len = " << pattern_len << std::endl;
        std::cerr << "min_occ = "<< min_occ << " and max_occ = " << max_occ << std::endl;
    }


    m_patterns.resize(m_pattern_cnt);
    if (candidates.size() > 0) {
        char* pattern_buf = new char[pattern_len+1];
        pattern_buf[pattern_len] = '\0';
        ifstream text_stream(text_file_name);
        if (text_stream) {
            for (size_type i=0, j, k; i < pattern_cnt; ++i) {
                j=rand()%candidates.size();
                if (already_used[j]) {
                    // find next free location
                    size_type j1=j+1;
                    while (j1 < already_used.size() and already_used[j1]) {
                        ++j1;
                    }
                    if (j1 >= already_used.size()) {
                        j1 = j;
                    }
                    j = j1;
                }
                k = cst.csa[ candidates[j] ];
                // extract the pattern
                text_stream.seekg(k ,std::ios::beg);
                text_stream.read(pattern_buf, pattern_len);

                m_patterns[i] = std::string(pattern_buf);
                already_used[j] = 1;
            }
            text_stream.close();
        } else {
            std::cerr << "ERROR: Could not open text file " << text_file_name << endl;
        }
        delete [] pattern_buf;
    }

    return m_pattern_cnt;
}

} // end of namespace sdsl

#endif
