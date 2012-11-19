/* sdsl - succinct data structures library
    Copyright (C) 2012 Matthias Petri

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
/*! \file select_support_ggmn.hpp
    \brief select_support_ggmn.hpp contains classes that support a sdsl::bit_vector with log(n) select operations.
	\author Matthias Petri
*/
#ifndef INCLUDED_SDSL_SELECT_SUPPORT_GGMN
#define INCLUDED_SDSL_SELECT_SUPPORT_GGMN

#include "select_support.hpp"
#include "int_vector.hpp"

#include <queue>

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! A class supporting select queries in log(n) time.
/*! The implementation is a version of the data structure
 *  proposed by Rodrigo Gonz\'alez and
             Szymon Grabowski and
             Veli M{\"akinen} and
             Gonzalo Navarro. in 2005.
 *
 * @ingroup select_support_group
 */
class select_support_ggmn : public select_support
{
    public:
        typedef select_support::size_type   size_type;
        typedef bit_vector bit_vector_type;
    private:
        size_type m_logn;
        int_vector<0> m_superblockrank;
        int_vector<0> m_blockrank;
        int_vector<64> m_rank_samples;   /* space for additional rank samples */
    public:
        explicit select_support_ggmn(const int_vector<1>* v=NULL);
        select_support_ggmn(const select_support_ggmn& rs);
        ~select_support_ggmn();
        void init(const int_vector<1>* v=NULL);
        inline const size_type select(size_type idx) const;
        inline const size_type operator()(size_type idx)const;
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const;
        void load(std::istream& in, const int_vector<1>* v=NULL);
        void set_vector(const int_vector<1>* v);
        //! Assign Operator
        /*! Required for the Assignable Concept of the STL.
         */
        select_support_ggmn& operator=(const select_support_ggmn& ss);
        //! swap Operator
        /*! Swap two select_support_ggmn in constant time.
            Required for the Container Concept of the STL.
         */
        void swap(select_support_ggmn& ss);

        // precondition: m_rank_samples.size() <= m_superblocks
        void init_rank_samples() {
            size_type idx = 0;
            std::queue<size_type> lbs, rbs;
            lbs.push(0); rbs.push(m_superblockrank.size()-1);
            while (!lbs.empty()) {
                size_type lb = lbs.front(); lbs.pop();
                size_type rb = rbs.front(); rbs.pop();
                if (/*lb < rb and*/ idx < m_rank_samples.size()) {
                    size_type mid = lb + (rb-lb)/2; // select mid \in [lb..rb)
                    m_rank_samples[ idx++ ] = m_superblockrank[mid];
                    lbs.push(lb); rbs.push(mid);
                    lbs.push(mid+1); rbs.push(rb);
                }
            }
        }

};

inline select_support_ggmn::select_support_ggmn(const int_vector<1>* v)
{
    m_logn = 0;
    init(v);
}

inline select_support_ggmn::select_support_ggmn(const select_support_ggmn& ss) : select_support()
{
    set_vector(ss.m_v);
    m_superblockrank = ss.m_superblockrank;
    m_blockrank		 = ss.m_blockrank;
    m_rank_samples   = ss.m_rank_samples;
}

inline select_support_ggmn& select_support_ggmn::operator=(const select_support_ggmn& ss)
{
    if (this != &ss) {
        set_vector(ss.m_v);
        m_superblockrank = ss.m_superblockrank;
        m_blockrank		= ss.m_blockrank;
        m_rank_samples   = ss.m_rank_samples;
    }
    return *this;
}

inline void select_support_ggmn::swap(select_support_ggmn& ss)
{
    if (this != &ss) { // if ss and _this_ are not the same object
        std::swap(m_logn, ss.m_logn);
        m_superblockrank.swap(ss.m_superblockrank);
        m_blockrank.swap(ss.m_blockrank);
        m_rank_samples.swap(ss.m_rank_samples);
    }
}

inline void select_support_ggmn::init(const int_vector<1>* v)
{
    set_vector(v);
    if (m_v == NULL) return;
    if (m_v->empty()) {
        m_blockrank.set_int_width(1); m_superblockrank.set_int_width(1);
        m_blockrank.resize(1);		m_superblockrank.resize(1);
        m_blockrank[0] = 0; m_superblockrank[0] = 0;
        return;
    }
    m_blockrank.set_int_width(12);
    m_blockrank.resize((m_v->capacity()>>6) + (0==(m_v->size()&0x3F)));     // n/64 + 2*loglog 64
    m_superblockrank.set_int_width(m_logn);
    m_superblockrank.resize((m_blockrank.size()+63)>>6);

    m_blockrank[0]=0;
    m_superblockrank[0]=0;
    size_type cnt = 0, blockcnt = 0, wcnt = 0;
    const uint64_t* data = m_v->data();
    size_type i;
    for (i = 1; i < (m_v->capacity()>>6) ; ++i) {
        wcnt = bit_magic::b1Cnt(*data);
        ++data;
        blockcnt += wcnt;
        cnt 	 += wcnt;
        if ((i & 0x3F) == 0) {
            m_superblockrank[i>>6] = cnt;
            blockcnt = 0;
        }
        m_blockrank[i] = blockcnt;
    }
    if (0 == (m_v->size()&0x3F)) {
        wcnt = bit_magic::b1Cnt(*data);
        blockcnt += wcnt;
        cnt		 += wcnt;
        if ((i & 0x3F) == 0) {
            m_superblockrank[i>>6] = cnt;
            blockcnt = 0;
        }
        m_blockrank[i] = blockcnt;
    }

    if (m_superblockrank.size() > 512) {
        // we store at most m_superblocks+1 rank_samples:
        // we do a cache efficient binary search for the select on X=1024
        // or X=the smallest power of two smaller than m_superblock
        m_rank_samples.resize(std::min(1024ULL, 1ULL << bit_magic::l1BP(m_superblockrank.size())));
    }
    init_rank_samples();
}


inline const select_support_ggmn::size_type select_support_ggmn::select(size_type idx) const
{
    size_type lb = 0, rb = m_superblockrank.size()-1; // search interval [lb..rb)
    size_type res = 0;

    /* binary search over super blocks */
    while (lb < rb) {
        size_type mid = (lb+rb)/2; // select mid \in [lb..rb)
        if (m_superblockrank[mid] >= idx)
            rb = mid;
        else
            lb = mid + 1;
    }

    const uint64_t* w = m_v->data();
    if (rb) {
        res = (rb-1) << 12; // *4096
        w = w + ((rb-1)<<6);
        idx -= m_superblockrank[rb-1];  // subtract the cumulative sum before the superblock
    }
    size_type ones = bit_magic::b1Cnt(*w);
    while (ones < idx) {
        idx -= ones; ++w;
        ones = bit_magic::b1Cnt(*w);
        res += 64;
    }
    /* handle last word */
    res += bit_magic::i1BP(*w, idx);

    return res;
}

inline const select_support_ggmn::size_type select_support_ggmn::operator()(size_type idx) const
{
    return select(idx);
}

inline select_support_ggmn::size_type select_support_ggmn::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    size_type written_bytes = 0;
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    written_bytes += m_blockrank.serialize(out, child, "blockrank");
    written_bytes += m_superblockrank.serialize(out, child, "superblockrank");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

inline void select_support_ggmn::load(std::istream& in, const int_vector<1>* v)
{
    set_vector(v);
    assert(m_v != NULL); // supported bit vector should be known
    m_blockrank.load(in);
    m_superblockrank.load(in);
}

inline select_support_ggmn::~select_support_ggmn()
{
}

inline void select_support_ggmn::set_vector(const int_vector<1>* v)
{
    if (v != NULL) {
        m_v = v;
        m_logn = bit_magic::l1BP(m_v->capacity())+1;
    }
}


}// end namespace sds

#endif // end file
