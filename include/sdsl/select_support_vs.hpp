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
/*! \file select_support_vs.hpp
    \brief select_support_vs.hpp contains classes that support a sdsl::bit_vector with constant time select information.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_SELECT_SUPPORT_VS
#define INCLUDED_SDSL_SELECT_SUPPORT_VS

#include "int_vector.hpp"
#include "select_support.hpp"

//#define SDSL_DEBUG_SELECT_SUPPORT_JMC

#ifdef SDSL_DEBUG_SELECT_SUPPORT_VS
#include "testutils.hpp"
#endif

#include <algorithm>

#define VS_MAX_ONES_PER_INVENTORY (8192)
#define VS_MAX_LOG2_LONGWORDS_PER_SUBINVENTORY (3)

#define VS_ONES_STEP_4 ( 0x1111111111111111ULL )
#define VS_VS_ONES_STEP_8 ( 0x0101010101010101ULL )
#define VS_ONES_STEP_9 ( 1ULL << 0 | 1ULL << 9 | 1ULL << 18 | 1ULL << 27 | 1ULL << 36 | 1ULL << 45 | 1ULL << 54 )
#define VS_ONES_STEP_16 ( 1ULL << 0 | 1ULL << 16 | 1ULL << 32 | 1ULL << 48 )
#define VS_MSBS_STEP_4 ( 0x8ULL * VS_ONES_STEP_4 )
#define VS_VS_MSBS_STEP_8 ( 0x80ULL * VS_VS_ONES_STEP_8 )
#define VS_VS_MSBS_STEP_9 ( 0x100ULL * VS_ONES_STEP_9 )
#define VS_VS_MSBS_STEP_16 ( 0x8000ULL * VS_ONES_STEP_16 )
#define VS_INCR_STEP_8 ( 0x80ULL << 56 | 0x40ULL << 48 | 0x20ULL << 40 | 0x10ULL << 32 | 0x8ULL << 24 | 0x4ULL << 16 | 0x2ULL << 8 | 0x1 )

#define VS_ONES_STEP_32 ( 0x0000000100000001ULL )
#define VS_MSBS_STEP_32 ( 0x8000000080000000ULL )

#define VS_COMPARE_STEP_8(x,y) ( ( ( ( ( (x) | VS_MSBS_STEP_8 ) - ( (y) & ~VS_MSBS_STEP_8 ) ) ^ (x) ^ ~(y) ) & VS_MSBS_STEP_8 ) >> 7 )
#define VS_LEQ_STEP_8(x,y) ( ( ( ( ( (y) | VS_MSBS_STEP_8 ) - ( (x) & ~VS_MSBS_STEP_8 ) ) ^ (x) ^ (y) ) & VS_MSBS_STEP_8 ) >> 7 )

#define VS_UCOMPARE_STEP_9(x,y) ( ( ( ( ( ( (x) | VS_MSBS_STEP_9 ) - ( (y) & ~VS_MSBS_STEP_9 ) ) | ( x ^ y ) ) ^ ( x | ~y ) ) & VS_MSBS_STEP_9 ) >> 8 )
#define VS_UCOMPARE_STEP_16(x,y) ( ( ( ( ( ( (x) | VS_MSBS_STEP_16 ) - ( (y) & ~VS_MSBS_STEP_16 ) ) | ( x ^ y ) ) ^ ( x | ~y ) ) & VS_MSBS_STEP_16 ) >> 15 )
#define VS_ULEQ_STEP_9(x,y) ( ( ( ( ( ( (y) | VS_MSBS_STEP_9 ) - ( (x) & ~VS_MSBS_STEP_9 ) ) | ( x ^ y ) ) ^ ( x & ~y ) ) & VS_MSBS_STEP_9 ) >> 8 )
#define VS_ULEQ_STEP_16(x,y) ( ( ( ( ( ( (y) | VS_MSBS_STEP_16 ) - ( (x) & ~VS_MSBS_STEP_16 ) ) | ( x ^ y ) ) ^ ( x & ~y ) ) & VS_MSBS_STEP_16 ) >> 15 )
#define VS_ZCOMPARE_STEP_8(x) ( ( ( x | ( ( x | VS_MSBS_STEP_8 ) - VS_ONES_STEP_8 ) ) & VS_MSBS_STEP_8 ) >> 7 )

#define VS_EASY_LEQ_STEP_8(x,y) ( ( ( ( ( (y) | VS_MSBS_STEP_8 ) - ( x ) ) ) & VS_MSBS_STEP_8 ) >> 7 )

//! Namespace for the succinct data structure library.
namespace sdsl
{


template<uint8_t bit_pattern, uint8_t pattern_len>
struct select_support_vs_trait {
    typedef select_support::size_type	size_type;

    /* Count the number of arguments for the specific select support */
    static size_type arg_cnt(const bit_vector&) {
        return 0;
    }

    static size_type args_in_the_first_word(uint64_t, uint8_t, uint64_t) {
        return 0;
    }

    static size_type ith_arg_pos_in_the_first_word(uint64_t, size_type, uint8_t, uint64_t) {
        return 0;
    }

    static size_type args_in_the_word(uint64_t, uint64_t&) {
        return 0;
    }

    static size_type ith_arg_pos_in_the_word(uint64_t, size_type, uint64_t) {
        return 0;
    }

    static bool found_arg(size_type, const bit_vector&) {
        return 0;
    }

    static uint64_t init_carry(const uint64_t*, size_type) {
        return 0;
    }
    static uint64_t get_carry(uint64_t) {
        return 0;
    }
};

template<>
struct select_support_vs_trait<0,1> {
    typedef select_support::size_type	size_type;

    static size_type arg_cnt(const bit_vector& v) {
        return v.bit_size()-util::get_one_bits(v);
    }

    static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t) {
        return bit_magic::b1Cnt((~w)& bit_magic::Li0Mask[offset]);
    }

    static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t) {
        return bit_magic::i1BP(~w & bit_magic::Li0Mask[offset], i);
    }
    static size_type args_in_the_word(uint64_t w, uint64_t&) {
        return bit_magic::b1Cnt(~w);
    }
    static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t) {
        return bit_magic::i1BP(~w, i);
    }
    static bool found_arg(size_type i, const bit_vector& v) {
        return !v[i];
    }
    static uint64_t init_carry(const uint64_t*, size_type) {
        return 0;
    }
    static uint64_t get_carry(uint64_t) {
        return 0;
    }
};

template<>
struct select_support_vs_trait<1,1> {
    typedef select_support::size_type	size_type;

    /* Count the number of arguments for the specific select support */
    static size_type arg_cnt(const bit_vector& v) {
        return util::get_one_bits(v);
    }

    static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t) {
        return bit_magic::b1Cnt(w & bit_magic::Li0Mask[offset]);
    }

    static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t) {
        return bit_magic::i1BP(w & bit_magic::Li0Mask[offset], i);
    }

    static size_type args_in_the_word(uint64_t w, uint64_t&) {
        return bit_magic::b1Cnt(w);
    }

    static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t) {
        return bit_magic::i1BP(w, i);
    }

    static bool found_arg(size_type i, const bit_vector& v) {
        return v[i];
    }

    static uint64_t init_carry(const uint64_t*, size_type) {
        return 0;
    }
    static uint64_t get_carry(uint64_t) {
        return 0;
    }
};

template<>
struct select_support_vs_trait<10,2> {
    typedef select_support::size_type	size_type;

    /* Count the number of arguments for the specific select support */
    static size_type arg_cnt(const bit_vector& v) {
        return util::get_onezero_bits(v);
    }

    static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t carry) {
        return bit_magic::b1Cnt(bit_magic::b10Map(w, carry) & bit_magic::Li0Mask[offset]);
    }

    static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t carry) {
        return bit_magic::i1BP(bit_magic::b10Map(w, carry) & bit_magic::Li0Mask[offset], i);
    }

    static size_type args_in_the_word(uint64_t w, uint64_t& carry) {
        return bit_magic::b10Cnt(w, carry);
    }

    static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t carry) {
        return bit_magic::i1BP(bit_magic::b10Map(w, carry), i);
    }

    static bool found_arg(size_type i, const bit_vector& v) {
        if (i > 0 and v[i-1] and !v[i])
            return true;
        return false;
    }

    static uint64_t init_carry(const uint64_t* data, size_type word_pos) {
        return word_pos ? (*(data-1)>>63) : 0;
    }
    static uint64_t get_carry(uint64_t w) {
        return w>>63;
    }
};

template<>
struct select_support_vs_trait<01,2> {
    typedef select_support::size_type	size_type;

    /* Count the number of arguments for the specific select support */
    static size_type arg_cnt(const bit_vector& v) {
        return util::get_zeroone_bits(v);
    }

    static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t carry) {
        return bit_magic::b1Cnt(bit_magic::b01Map(w, carry) & bit_magic::Li0Mask[offset]);
    }

    static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t carry) {
        return bit_magic::i1BP(bit_magic::b01Map(w, carry) & bit_magic::Li0Mask[offset], i);
    }

    static size_type args_in_the_word(uint64_t w, uint64_t& carry) {
        return bit_magic::b01Cnt(w, carry);
    }

    static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t carry) {
        return bit_magic::i1BP(bit_magic::b01Map(w, carry), i);
    }

    static bool found_arg(size_type i, const bit_vector& v) {
        if (i > 0 and !v[i-1] and v[i])
            return true;
        return false;
    }

    static uint64_t init_carry(const uint64_t* data, size_type word_pos) {
        return word_pos ? (*(data-1)>>63) : 1;
    }
    static uint64_t get_carry(uint64_t w) {
        return w>>63;
    }
};


//! A class supporting constant time select queries (proposed by Munro/Clark, 1996) enhanced by broadword computing tricks.
/*!
 * \par Space usage
 *      The space usage of the data structure depends on the number of \f$ m \f$ of ones in the
 *      original bitvector $b$. We store the position of every $4096$th set bit
 *      (called L1-sampled bits) of $b$.
 *      This takes in the worst case \f$\frac{m}{4096} \log{n} \leq \frac{64}{n}\f$ bits.
 *      Next,
 *      (1) if the distance of two adjacent L1-sampled bits $b[i]$ and $b[j]$
 *      is greater or equal than $\log^4 n$, then
 *      we store each of the 4096 positions of the set $b$ in [i..j-1] with
 *      $\log{n}$ bits. This results in at most
 *      \$ \frac{4096\cdot \log n}{\log^4 n}=\frac{4096}{\log^3 n}\$ bits per bit.
 *      For a bitvector of 1GB, i.e. \f$ \log n = 35 \f$ we get about 0.01 bits per bit.
 *      If the $j-i+1 < \log^4 n$ then
 *      (2) we store the relative position of every $64$th set bit (called L2-sampled bits)
 *      in b[i..j-1] in at most $4\log\log n$ bits per L2-sampled bits.
 *      An pessimistic upper bound for the space would be
 *      \f$ \frac{4\log\log n}{64} \leq \frac{24}{64} = 0.375\f$ bit per
 *      bit (since $\log\log n\leq 6$. It is very pessimistic, since we store
 *      the relative position in $\log\log(j-i+1)\leq \log\log n$ bits.
 * @ingroup select_support_group
 */
template<uint8_t b=1, uint8_t pattern_len=1>
class select_support_vs : public select_support
{
    public:
        typedef bit_vector bit_vector_type;
    public:
        int_vector<64> m_inventory;
        int_vector<64> m_subinventory;
        int_vector<64> m_exact_spill;
        int log2_ones_per_inventory, log2_ones_per_sub16, log2_ones_per_sub64, log2_longwords_per_subinventory,
            ones_per_inventory, ones_per_sub16, ones_per_sub64, longwords_per_subinventory,
            ones_per_inventory_mask, ones_per_sub16_mask, ones_per_sub64_mask;

        uint64_t num_words, inventory_size, subinventory_size, exact_spill_size, num_ones;

        void copy(const select_support_vs<b, pattern_len>& ss);
    public:
        explicit select_support_vs(const int_vector<1>* v=NULL);
        select_support_vs(const select_support_vs<b,pattern_len>& ss);
        ~select_support_vs();
        void init(const int_vector<1>* v=NULL);
        //! Select function
        /*! \sa select_support.select
         */
        inline const size_type select(size_type i) const;
        //! Alias for select(i).
        inline const size_type operator()(size_type i)const;
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const;
        void load(std::istream& in, const int_vector<1>* v=NULL);
        void set_vector(const int_vector<1>* v=NULL);
        select_support_vs<b, pattern_len>& operator=(const select_support_vs& ss);

        //! Swap operator
        /*! This swap operator swaps two select_support_vss in constant time.
         */
        void swap(select_support_vs<b, pattern_len>& ss);

};


template<uint8_t b, uint8_t pattern_len>
select_support_vs<b,pattern_len>::select_support_vs(const int_vector<1>* f_v):select_support(f_v)
{
    init(f_v);
}

template<uint8_t b, uint8_t pattern_len>
select_support_vs<b,pattern_len>::select_support_vs(const select_support_vs& ss):select_support(ss.m_v)
{
    copy(ss);
}

template<uint8_t b, uint8_t pattern_len>
select_support_vs<b, pattern_len>& select_support_vs<b,pattern_len>::operator=(const select_support_vs& ss)
{
    if (this != &ss) {
        copy(ss);
    }
    return *this;
}

template<uint8_t b, uint8_t pattern_len>
void select_support_vs<b,pattern_len>::swap(select_support_vs& ss)
{
    std::swap(log2_ones_per_inventory, ss.log2_ones_per_inventory);
    std::swap(inventory_size, ss.inventory_size);
    std::swap(ones_per_inventory_mask, ss.ones_per_inventory_mask);
    std::swap(log2_longwords_per_subinventory, ss.log2_longwords_per_subinventory);
    std::swap(log2_ones_per_sub16, ss.log2_ones_per_sub16);
    std::swap(ones_per_sub64, ss.ones_per_sub64);
    std::swap(ones_per_sub16, ss.ones_per_sub16);
    std::swap(ones_per_sub64_mask, ss.ones_per_sub64_mask);
    std::swap(ones_per_sub16_mask, ss.ones_per_sub16_mask);
    std::swap(num_ones, ss.num_ones);
    std::swap(num_ones, ss.num_ones);

    m_inventory.swap(ss.m_inventory);
    m_subinventory.swap(ss.m_subinventory);
    m_exact_spill.swap(ss.m_exact_spill);
}

template<uint8_t b, uint8_t pattern_len>
void select_support_vs<b,pattern_len>::copy(const select_support_vs<b, pattern_len>& ss)
{
    m_v             = ss.m_v;           // copy pointer to the supported bit vector
    m_inventory = ss.m_inventory;
    m_subinventory = ss.m_subinventory;
    m_exact_spill = ss.m_exact_spill;

    log2_ones_per_inventory          = ss.log2_ones_per_inventory;
    inventory_size          = ss.inventory_size;
    ones_per_inventory_mask          = ss.ones_per_inventory_mask;
    log2_longwords_per_subinventory          = ss.log2_longwords_per_subinventory;
    log2_ones_per_sub16          = ss.log2_ones_per_sub16;
    ones_per_sub64          = ss.ones_per_sub64;
    ones_per_sub16          = ss.ones_per_sub16;
    ones_per_sub64_mask          = ss.ones_per_sub64_mask;
    ones_per_sub16_mask          = ss.ones_per_sub16_mask;
    num_ones          = ss.num_ones;
    num_words          = ss.num_words;

}

template<uint8_t b, uint8_t pattern_len>
select_support_vs<b,pattern_len>::~select_support_vs()
{
}

template<uint8_t b, uint8_t pattern_len>
void select_support_vs<b,pattern_len>::init(const int_vector<1>* v)
{
    //fprintf(stderr, "INIT_VS()\n");
    set_vector(v);

    if (m_v==NULL)
        return;

    /*    fprintf(stdout, "B[] = ");
        for(size_t i=0;i<v->size();i++) {
            int bit = (*v)[i];
            fprintf(stdout, "%d",bit);
        }
        fprintf(stdout, "\n");*/

    uint64_t num_bits = v->size();
    num_words = (num_bits + 63) / 64;

    num_ones = select_support_vs_trait<b,pattern_len>::arg_cnt(*v);

    ones_per_inventory = num_bits == 0 ? 0 : (num_ones * VS_MAX_ONES_PER_INVENTORY + num_bits - 1) / num_bits;

    // Make ones_per_inventory into a power of 2
    log2_ones_per_inventory = bit_magic::l1BP(ones_per_inventory);
    ones_per_inventory = 1ULL << log2_ones_per_inventory;
    ones_per_inventory_mask = ones_per_inventory - 1;
    inventory_size = (num_ones + ones_per_inventory - 1) / ones_per_inventory;

    //printf("Number of ones: %ld Number of ones per inventory item: %d\n", num_ones, ones_per_inventory );

    m_inventory.resize(inventory_size + 1);
    uint64_t* inventory = (uint64_t*) m_inventory.data();

    uint64_t d = 0;
    // First phase: we build an inventory for each one out of ones_per_inventory.
    const uint64_t* bits = m_v->data();
    for (uint64_t i = 0; i < num_words; i++) {
        for (int j = 0; j < 64; j++) {
            if (i * 64 + j >= num_bits) break;
            if (bits[ i ] & 1ULL << j) {
                if ((d & ones_per_inventory_mask) == 0) inventory[ d >> log2_ones_per_inventory ] = i * 64 + j;
                d++;
            }
        }
    }

    assert(num_ones == d);
    inventory[ inventory_size ] = num_bits;

    //printf("Inventory entries filled: %ld\n", inventory_size + 1 );

    log2_longwords_per_subinventory = std::min(VS_MAX_LOG2_LONGWORDS_PER_SUBINVENTORY, std::max(0, log2_ones_per_inventory - 2));
    longwords_per_subinventory = 1 << log2_longwords_per_subinventory;
    log2_ones_per_sub64 = std::max(0, log2_ones_per_inventory - log2_longwords_per_subinventory);
    log2_ones_per_sub16 = std::max(0, log2_ones_per_sub64 - 2);
    ones_per_sub64 = (1ULL << log2_ones_per_sub64);
    ones_per_sub16 = (1ULL << log2_ones_per_sub16);
    ones_per_sub64_mask = ones_per_sub64 - 1;
    ones_per_sub16_mask = ones_per_sub16 - 1;

    //printf("Longwords per subinventory: %d Ones per sub 64: %d sub 16: %d\n", longwords_per_subinventory, ones_per_sub64, ones_per_sub16 );

    if (ones_per_inventory > 1) {
        d = 0;
        int ones = 0;
        uint64_t spilled = 0, diff16 = 0, exact = 0, start = 0, span = 0, inventory_index;

        for (uint64_t i = 0; i < num_words; i++)
            // We estimate the subinventory and exact spill size
            for (int j = 0; j < 64; j++) {
                if (i * 64 + j >= num_bits) break;
                if (bits[ i ] & 1ULL << j) {
                    if ((d & ones_per_inventory_mask) == 0) {
                        inventory_index = d >> log2_ones_per_inventory;
                        start = inventory[ inventory_index ];
                        span = inventory[ inventory_index + 1 ] - start;
                        ones = std::min(num_ones - d, (uint64_t)ones_per_inventory);

                        // We must always count (possibly unused) diff16's. And we cannot store less then 4 diff16.
                        diff16 += std::max(4, (ones + ones_per_sub16 - 1) >> log2_ones_per_sub16);

                        // We accumulate space for exact pointers ONLY if necessary.
                        if (span >= (1<<16)) {
                            exact += ones;
                            if (ones_per_sub64 > 1) spilled += ones;
                        }

                    }
                    d++;
                }
            }

        //printf("Spilled entries: %lld exact: %lld diff16: %lld subinventory size: %lld\n", spilled, exact, diff16 - ( exact - spilled ) * 4, ( diff16 + 3 ) / 4 );

        subinventory_size = (diff16 + 3) / 4;
        exact_spill_size = spilled;
        m_subinventory.resize(subinventory_size);
        m_exact_spill.resize(exact_spill_size);

        uint16_t* p16 = NULL;
        uint64_t* p64 = NULL;
        int offset = 0;
        spilled = 0;
        d = 0;

        uint64_t* subinventory = (uint64_t*) m_subinventory.data();
        uint64_t* exact_spill = (uint64_t*) m_exact_spill.data();
        for (uint64_t i = 0; i < num_words; i++)
            for (int j = 0; j < 64; j++) {
                if (i * 64 + j >= num_bits) break;
                if (bits[ i ] & 1ULL << j) {
                    if ((d & ones_per_inventory_mask) == 0) {
                        inventory_index = d >> log2_ones_per_inventory;
                        start = inventory[ inventory_index ];
                        span = inventory[ inventory_index + 1 ] - start;
                        p16 = (uint16_t*)&subinventory[ inventory_index << log2_longwords_per_subinventory ];
                        p64 = &subinventory[ inventory_index << log2_longwords_per_subinventory ];
                        offset = 0;
                    }

                    if (span < (1<<16)) {
                        assert(i * 64 + j - start <= (1<<16));
                        if ((d & ones_per_sub16_mask) == 0) {
                            assert(offset < longwords_per_subinventory * 4);
                            assert(p16 + offset < (uint16_t*)(subinventory + subinventory_size));
                            p16[ offset++ ] = i * 64 + j - start;
                        }
                    } else {
                        if (ones_per_sub64 == 1) {
                            assert(p64 + offset < subinventory + subinventory_size);
                            p64[ offset++ ] = i * 64 + j;
                        } else {
                            assert(p64 < subinventory + subinventory_size);
                            if ((d & ones_per_inventory_mask) == 0) {
                                inventory[ inventory_index ] |= 1ULL << 63;
                                p64[ 0 ] = spilled;
                            }
                            assert(spilled < exact_spill_size);
                            exact_spill[ spilled++ ] = i * 64 + j;
                        }
                    }

                    d++;
                }
            }
    } else {
        exact_spill_size = subinventory_size = 0;
    }
}

template<uint8_t b, uint8_t pattern_len>
inline const typename select_support_vs<b,pattern_len>::size_type select_support_vs<b,pattern_len>::select(size_type rank)const
{
    rank--; // vigna starts with 0. we start with 1

//#ifdef DEBUG
    // printf( "Selecting %ld\n...", rank );
//#endif
    const uint64_t* subinventory = m_subinventory.data();
    const uint64_t* exact_spill = m_exact_spill.data();
    const uint64_t* inventory = m_inventory.data();
    const uint64_t* bits = m_v->data();

    const uint64_t inventory_index = rank >> log2_ones_per_inventory;
    assert(inventory_index < inventory_size);

    const int64_t inventory_rank = inventory[ inventory_index ];
    const int subrank = rank & ones_per_inventory_mask;
//#ifdef DEBUG
    // printf( "Rank: %ld inventory index: %ld inventory rank: %ld subrank: %d\n", rank, inventory_index, inventory_rank, subrank );
//#endif

//#ifdef DEBUG
    // if ( subrank == 0 ) puts( "Exact hit (no subrank); returning inventory" );
//#endif
    if (subrank == 0) return inventory_rank & ~(1ULL<<63);

    uint64_t start;
    int residual;

    if (inventory_rank >= 0) {
        start = inventory_rank + ((uint16_t*)(subinventory + (inventory_index << log2_longwords_per_subinventory)))[ subrank >> log2_ones_per_sub16 ];
        residual = subrank & ones_per_sub16_mask;
    } else {
        if (ones_per_sub64 == 1) return subinventory[(inventory_index << log2_longwords_per_subinventory) + subrank ];
        assert(subinventory[ inventory_index << log2_longwords_per_subinventory ] + subrank < exact_spill_size);
        return exact_spill[ subinventory[ inventory_index << log2_longwords_per_subinventory ] + subrank ];
    }

//#ifdef DEBUG
    // printf( "Differential; start: %ld residual: %d\n", start, residual );
    //if ( residual == 0 ) puts( "No residual; returning start" );
//#endif

    if (residual == 0) return start;

    uint64_t word_index = start / 64;
    uint64_t word = bits[ word_index ] & -1ULL << start;

    // fprintf(stdout, "word_index %ld\n",word_index);

    for (;;) {
        const int bit_count = bit_magic::b1Cnt(word);
        if (residual < bit_count) break;
        word = bits[ ++word_index ];
        residual -= bit_count;
    }

    return word_index * 64 + bit_magic::i1BP(word, residual + 1);
}

template<uint8_t b, uint8_t pattern_len>
inline const typename select_support_vs<b,pattern_len>::size_type select_support_vs<b,pattern_len>::operator()(size_type i)const
{
    return select(i);
}

template<uint8_t b, uint8_t pattern_len>
void select_support_vs<b,pattern_len>::set_vector(const int_vector<1>* v)
{
    m_v = v;
}

template<uint8_t b, uint8_t pattern_len>
typename select_support_vs<b,pattern_len>::size_type select_support_vs<b,pattern_len>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;

    written_bytes += m_inventory.serialize(out,child,"inventory");
    written_bytes += m_subinventory.serialize(out,child,"subinventory");
    written_bytes += m_exact_spill.serialize(out,child,"exact_spill");

    written_bytes += util::write_member(log2_ones_per_inventory, out, child, "log2_ones_per_inventory");
    written_bytes += util::write_member(inventory_size, out, child, "inventory_size");
    written_bytes += util::write_member(ones_per_inventory_mask, out, child, "ones_per_inventory_mask");
    written_bytes += util::write_member(log2_longwords_per_subinventory, out, child, "log2_longwords_per_subinventory");
    written_bytes += util::write_member(log2_ones_per_sub16, out, child, "log2_ones_per_sub16");
    written_bytes += util::write_member(ones_per_sub64, out, child, "ones_per_sub64");
    written_bytes += util::write_member(ones_per_sub16, out, child, "ones_per_sub16");
    written_bytes += util::write_member(ones_per_sub64_mask, out, child, "ones_per_sub64_mask");
    written_bytes += util::write_member(ones_per_sub16_mask, out, child, "ones_per_sub16_mask");
    written_bytes += util::write_member(num_ones, out, child, "num_ones");
    written_bytes += util::write_member(num_words, out, child, "num_words");

    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<uint8_t b, uint8_t pattern_len>
void select_support_vs<b,pattern_len>::load(std::istream& in, const int_vector<1>* v)
{
    set_vector(v);

    m_inventory.load(in);
    m_subinventory.load(in);
    m_exact_spill.load(in);

    util::read_member(log2_ones_per_inventory, in);
    util::read_member(inventory_size, in);
    util::read_member(ones_per_inventory_mask, in);
    util::read_member(log2_longwords_per_subinventory, in);
    util::read_member(log2_ones_per_sub16, in);
    util::read_member(ones_per_sub64, in);
    util::read_member(ones_per_sub16, in);
    util::read_member(ones_per_sub64_mask, in);
    util::read_member(ones_per_sub16_mask, in);
    util::read_member(num_ones, in);
    util::read_member(num_words, in);
}



}

#endif
