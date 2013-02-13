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
/*! \file select_support_v9.hpp
    \brief select_support_v9.hpp contains classes that support a sdsl::bit_vector with constant time select information.
	\author Simon Gog
*/
#ifndef INCLUDED_SDSL_select_support_v9
#define INCLUDED_SDSL_select_support_v9

#include "int_vector.hpp"
#include "select_support.hpp"

//#define SDSL_DEBUG_SELECT_SUPPORT_JMC

#ifdef SDSL_DEBUG_select_support_v9
#include "testutils.hpp"
#endif

#define ONES_PER_INVENTORY (512)
#define LOG2_ONES_PER_INVENTORY (9)
#define INVENTORY_MASK (ONES_PER_INVENTORY-1)

#define ONES_STEP_4 ( 0x1111111111111111ULL )
#define ONES_STEP_8 ( 0x0101010101010101ULL )
#define ONES_STEP_9 ( 1ULL << 0 | 1ULL << 9 | 1ULL << 18 | 1ULL << 27 | 1ULL << 36 | 1ULL << 45 | 1ULL << 54 )
#define ONES_STEP_16 ( 1ULL << 0 | 1ULL << 16 | 1ULL << 32 | 1ULL << 48 )
#define MSBS_STEP_4 ( 0x8ULL * ONES_STEP_4 )
#define MSBS_STEP_8 ( 0x80ULL * ONES_STEP_8 )
#define MSBS_STEP_9 ( 0x100ULL * ONES_STEP_9 )
#define MSBS_STEP_16 ( 0x8000ULL * ONES_STEP_16 )
#define INCR_STEP_8 ( 0x80ULL << 56 | 0x40ULL << 48 | 0x20ULL << 40 | 0x10ULL << 32 | 0x8ULL << 24 | 0x4ULL << 16 | 0x2ULL << 8 | 0x1 )

#define ONES_STEP_32 ( 0x0000000100000001ULL )
#define MSBS_STEP_32 ( 0x8000000080000000ULL )

#define COMPARE_STEP_8(x,y) ( ( ( ( ( (x) | MSBS_STEP_8 ) - ( (y) & ~MSBS_STEP_8 ) ) ^ (x) ^ ~(y) ) & MSBS_STEP_8 ) >> 7 )
#define LEQ_STEP_8(x,y) ( ( ( ( ( (y) | MSBS_STEP_8 ) - ( (x) & ~MSBS_STEP_8 ) ) ^ (x) ^ (y) ) & MSBS_STEP_8 ) >> 7 )

#define UCOMPARE_STEP_9(x,y) ( ( ( ( ( ( (x) | MSBS_STEP_9 ) - ( (y) & ~MSBS_STEP_9 ) ) | ( x ^ y ) ) ^ ( x | ~y ) ) & MSBS_STEP_9 ) >> 8 )
#define UCOMPARE_STEP_16(x,y) ( ( ( ( ( ( (x) | MSBS_STEP_16 ) - ( (y) & ~MSBS_STEP_16 ) ) | ( x ^ y ) ) ^ ( x | ~y ) ) & MSBS_STEP_16 ) >> 15 )
#define ULEQ_STEP_9(x,y) ( ( ( ( ( ( (y) | MSBS_STEP_9 ) - ( (x) & ~MSBS_STEP_9 ) ) | ( x ^ y ) ) ^ ( x & ~y ) ) & MSBS_STEP_9 ) >> 8 )
#define ULEQ_STEP_16(x,y) ( ( ( ( ( ( (y) | MSBS_STEP_16 ) - ( (x) & ~MSBS_STEP_16 ) ) | ( x ^ y ) ) ^ ( x & ~y ) ) & MSBS_STEP_16 ) >> 15 )
#define ZCOMPARE_STEP_8(x) ( ( ( x | ( ( x | MSBS_STEP_8 ) - ONES_STEP_8 ) ) & MSBS_STEP_8 ) >> 7 )

#define EASY_LEQ_STEP_8(x,y) ( ( ( ( ( (y) | MSBS_STEP_8 ) - ( x ) ) ) & MSBS_STEP_8 ) >> 7 )

//! Namespace for the succinct data structure library.
namespace sdsl
{


template<uint8_t bit_pattern, uint8_t pattern_len>
struct select_support_v9_trait {
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
struct select_support_v9_trait<0,1> {
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
struct select_support_v9_trait<1,1> {
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
struct select_support_v9_trait<10,2> {
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
struct select_support_v9_trait<01,2> {
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


const uint8_t select_in_byte[] = {
    0,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    7,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
    4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,

    0,0,0,1,0,2,2,1,0,3,3,1,3,2,2,1,
    0,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
    0,5,5,1,5,2,2,1,5,3,3,1,3,2,2,1,
    5,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
    0,6,6,1,6,2,2,1,6,3,3,1,3,2,2,1,
    6,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
    6,5,5,1,5,2,2,1,5,3,3,1,3,2,2,1,
    5,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
    0,7,7,1,7,2,2,1,7,3,3,1,3,2,2,1,
    7,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
    7,5,5,1,5,2,2,1,5,3,3,1,3,2,2,1,
    5,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
    7,6,6,1,6,2,2,1,6,3,3,1,3,2,2,1,
    6,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,
    6,5,5,1,5,2,2,1,5,3,3,1,3,2,2,1,
    5,4,4,1,4,2,2,1,4,3,3,1,3,2,2,1,

    0,0,0,0,0,0,0,2,0,0,0,3,0,3,3,2,
    0,0,0,4,0,4,4,2,0,4,4,3,4,3,3,2,
    0,0,0,5,0,5,5,2,0,5,5,3,5,3,3,2,
    0,5,5,4,5,4,4,2,5,4,4,3,4,3,3,2,
    0,0,0,6,0,6,6,2,0,6,6,3,6,3,3,2,
    0,6,6,4,6,4,4,2,6,4,4,3,4,3,3,2,
    0,6,6,5,6,5,5,2,6,5,5,3,5,3,3,2,
    6,5,5,4,5,4,4,2,5,4,4,3,4,3,3,2,
    0,0,0,7,0,7,7,2,0,7,7,3,7,3,3,2,
    0,7,7,4,7,4,4,2,7,4,4,3,4,3,3,2,
    0,7,7,5,7,5,5,2,7,5,5,3,5,3,3,2,
    7,5,5,4,5,4,4,2,5,4,4,3,4,3,3,2,
    0,7,7,6,7,6,6,2,7,6,6,3,6,3,3,2,
    7,6,6,4,6,4,4,2,6,4,4,3,4,3,3,2,
    7,6,6,5,6,5,5,2,6,5,5,3,5,3,3,2,
    6,5,5,4,5,4,4,2,5,4,4,3,4,3,3,2,

    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
    0,0,0,0,0,0,0,4,0,0,0,4,0,4,4,3,
    0,0,0,0,0,0,0,5,0,0,0,5,0,5,5,3,
    0,0,0,5,0,5,5,4,0,5,5,4,5,4,4,3,
    0,0,0,0,0,0,0,6,0,0,0,6,0,6,6,3,
    0,0,0,6,0,6,6,4,0,6,6,4,6,4,4,3,
    0,0,0,6,0,6,6,5,0,6,6,5,6,5,5,3,
    0,6,6,5,6,5,5,4,6,5,5,4,5,4,4,3,
    0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,3,
    0,0,0,7,0,7,7,4,0,7,7,4,7,4,4,3,
    0,0,0,7,0,7,7,5,0,7,7,5,7,5,5,3,
    0,7,7,5,7,5,5,4,7,5,5,4,5,4,4,3,
    0,0,0,7,0,7,7,6,0,7,7,6,7,6,6,3,
    0,7,7,6,7,6,6,4,7,6,6,4,6,4,4,3,
    0,7,7,6,7,6,6,5,7,6,6,5,6,5,5,3,
    7,6,6,5,6,5,5,4,6,5,5,4,5,4,4,3,

    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,
    0,0,0,0,0,0,0,5,0,0,0,5,0,5,5,4,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,
    0,0,0,0,0,0,0,6,0,0,0,6,0,6,6,4,
    0,0,0,0,0,0,0,6,0,0,0,6,0,6,6,5,
    0,0,0,6,0,6,6,5,0,6,6,5,6,5,5,4,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
    0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,4,
    0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,5,
    0,0,0,7,0,7,7,5,0,7,7,5,7,5,5,4,
    0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,6,
    0,0,0,7,0,7,7,6,0,7,7,6,7,6,6,4,
    0,0,0,7,0,7,7,6,0,7,7,6,7,6,6,5,
    0,7,7,6,7,6,6,5,7,6,6,5,6,5,5,4,

    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,
    0,0,0,0,0,0,0,6,0,0,0,6,0,6,6,5,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
    0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,5,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
    0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,6,
    0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,6,
    0,0,0,7,0,7,7,6,0,7,7,6,7,6,6,5,

    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,
    0,0,0,0,0,0,0,7,0,0,0,7,0,7,7,6,

    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7
};

__inline int select_in_word( const uint64_t x, const int rank ) {
    // Phase 1: sums by byte
    uint64_t byte_sums = x - ( ( x & 0xa * ONES_STEP_4 ) >> 1 );
    byte_sums = ( byte_sums & 3 * ONES_STEP_4 ) + ( ( byte_sums >> 2 ) & 3 * ONES_STEP_4 );
    byte_sums = ( byte_sums + ( byte_sums >> 4 ) ) & 0x0f * ONES_STEP_8;
    byte_sums *= ONES_STEP_8;
    // Phase 2: compare each byte sum with rank
    const uint64_t rank_plus_1_step_8 = ( rank + 0 ) * ONES_STEP_8;
    const int place = ( EASY_LEQ_STEP_8( byte_sums, rank_plus_1_step_8 ) * ONES_STEP_8 >> 53 ) & ~0x7;

    // Phase 3: Compute the rank in the relevant byte and use lookup table.
    const int byte_rank = ( rank + 0 ) - ( ( byte_sums << 8 ) >> place & 0xFF );
    return place + select_in_byte[ x >> place & 0xFF | byte_rank << 8 ];
}


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
class select_support_v9 : public select_support
{
    public:
        typedef bit_vector bit_vector_type;
    public:
        int_vector<64> m_counts;
        int_vector<64> m_inventory;
        int_vector<64> m_subinventory;
        uint64_t num_words, num_counts, inventory_size, ones_per_inventory, log2_ones_per_inventory, num_ones;

        void copy(const select_support_v9<b, pattern_len>& ss);
    public:
        explicit select_support_v9(const int_vector<1>* v=NULL);
        select_support_v9(const select_support_v9<b,pattern_len>& ss);
        ~select_support_v9();
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
        select_support_v9<b, pattern_len>& operator=(const select_support_v9& ss);

        //! Swap operator
        /*! This swap operator swaps two select_support_v9s in constant time.
         */
        void swap(select_support_v9<b, pattern_len>& ss);

};


template<uint8_t b, uint8_t pattern_len>
select_support_v9<b,pattern_len>::select_support_v9(const int_vector<1>* f_v):select_support(f_v)
{
    init(f_v);
}

template<uint8_t b, uint8_t pattern_len>
select_support_v9<b,pattern_len>::select_support_v9(const select_support_v9& ss):select_support(ss.m_v)
{
    copy(ss);
}

template<uint8_t b, uint8_t pattern_len>
select_support_v9<b, pattern_len>& select_support_v9<b,pattern_len>::operator=(const select_support_v9& ss)
{
    if (this != &ss) {
        copy(ss);
    }
    return *this;
}

template<uint8_t b, uint8_t pattern_len>
void select_support_v9<b,pattern_len>::swap(select_support_v9& ss)
{
    std::swap(num_words, ss.num_words);
    std::swap(num_counts, ss.num_counts);
    std::swap(ones_per_inventory, ss.ones_per_inventory);
    std::swap(log2_ones_per_inventory, ss.log2_ones_per_inventory);
    std::swap(num_ones, ss.num_ones);

    m_counts.swap(ss.m_inventory);
    m_inventory.swap(ss.m_inventory);
    m_subinventory.swap(ss.m_subinventory);
}

template<uint8_t b, uint8_t pattern_len>
void select_support_v9<b,pattern_len>::copy(const select_support_v9<b, pattern_len>& ss)
{
    m_v             = ss.m_v;           // copy pointer to the supported bit vector
    m_inventory = ss.m_inventory;
    m_subinventory = ss.m_subinventory;
    m_counts = ss.m_counts;

    ones_per_inventory = ss.ones_per_inventory;
    log2_ones_per_inventory = ss.log2_ones_per_inventory;
    num_counts          = ss.num_counts;
    num_ones          = ss.num_ones;
    num_words          = ss.num_words;
}

template<uint8_t b, uint8_t pattern_len>
select_support_v9<b,pattern_len>::~select_support_v9()
{
}

template<uint8_t b, uint8_t pattern_len>
void select_support_v9<b,pattern_len>::init(const int_vector<1>* v)
{
    //fprintf(stderr, "INIT_VS()\n");
    set_vector(v);

    if (m_v==NULL)
        return;

    size_t num_bits = v->size();
    uint64_t* bits = (uint64_t*) v->data();
    num_words = ( num_bits + 63 ) / 64;
    num_counts = ( ( num_bits + 64 * 8 - 1 ) / ( 64 * 8 ) ) * 2;

    // Init rank/select structure
    m_counts.resize(num_counts + 1);
    uint64_t* counts = (uint64_t*) m_counts.data();
    memset( counts, 0, ( num_counts + 1 ) * sizeof *counts );

    uint64_t c = 0;
    uint64_t pos = 0;
    for( uint64_t i = 0; i < num_words; i += 8, pos += 2 ) {
        counts[ pos ] = c;
        c += __builtin_popcountll( bits[ i ] );
        for( int j = 1;  j < 8; j++ ) {
            counts[ pos + 1 ] |= ( c - counts[ pos ] ) << 9 * ( j - 1 );
            if ( i + j < num_words ) c += __builtin_popcountll( bits[ i + j ] );
        }
    }

    counts[ num_counts ] = c;
#ifdef DEBUG
    printf("Number of ones: %lld\n", c );
#endif
    assert( c <= num_bits );

    inventory_size = ( c + ONES_PER_INVENTORY - 1 ) / ONES_PER_INVENTORY;

#ifdef DEBUG
    printf("Number of ones per inventory item: %d\n", ONES_PER_INVENTORY );
    assert( ONES_PER_INVENTORY <= 8 * 64 );
#endif

    m_inventory.resize(inventory_size + 1);
    uint64_t* inventory = (uint64_t*) m_inventory.data();
    memset( inventory, 0, ( inventory_size + 1 ) * sizeof *inventory );
    m_subinventory.resize(( num_words + 3 ) / 4);
    uint64_t* subinventory = (uint64_t*) m_subinventory.data();
    memset( subinventory, 0, ( ( num_words + 3 ) / 4 ) * sizeof *subinventory );

    uint64_t d = 0;
    for( uint64_t i = 0; i < num_words; i++ )
        for( int j = 0; j < 64; j++ )
            if ( bits[ i ] & 1ULL << j ) {
                if ( ( d & INVENTORY_MASK ) == 0 ) {
                    inventory[ d >> LOG2_ONES_PER_INVENTORY ] = i * 64 + j;
                    assert( counts[ ( i / 8 ) * 2 ] <= d );
                    assert( counts[ ( i / 8 ) * 2 + 2 ] > d );
                }

                d++;
            }

    assert( c == d );
    inventory[ inventory_size ] = ( ( num_words + 3 ) & ~3ULL ) * 64;

#ifdef DEBUG
    printf("Inventory entries filled: %lld\n", d / ONES_PER_INVENTORY + 1 );


    printf("First inventories: %lld %lld %lld %lld\n", inventory[ 0 ], inventory[ 1 ], inventory[ 2 ], inventory[ 3 ] );
#endif

    d = 0;
    int state;
    uint64_t *s, first_bit, index, span, block_span, block_left, counts_at_start;

    for( uint64_t i = 0; i < num_words; i++ ) {
        for( int j = 0; j < 64; j++ ) {
            if ( bits[ i ] & 1ULL << j ) {
                if ( ( d & INVENTORY_MASK ) == 0 ) {
                    first_bit = i * 64 + j;
                    index = d >> LOG2_ONES_PER_INVENTORY;
                    assert( inventory[ index ] == first_bit );
                    s = &subinventory[ ( inventory[ index ] / 64 ) / 4 ];
                    span = ( inventory[ index + 1 ] / 64 ) / 4 - ( inventory[ index ] / 64 ) / 4;
                    state = -1;
                    counts_at_start = counts[ ( ( inventory[ index ] / 64 ) / 8 ) * 2 ];
                    block_span = ( inventory[ index + 1 ] / 64 ) / 8 - ( inventory[ index ] / 64 ) / 8;
                    block_left = ( inventory[ index ] / 64 ) / 8;

                    if ( span >= 512 ) state = 0;
                    else if ( span >= 256 ) state = 1;
                    else if ( span >= 128 ) state = 2;
                    else if ( span >= 16 ) {
                        assert( ( block_span + 8 & -8LL ) + 8 <= span * 4 );

                        int k;
                        for( k = 0; k < block_span; k++ ) {
                            assert( ((uint16_t *)s)[ k + 8 ] == 0 );
                            ((uint16_t *)s)[ k + 8 ] = counts[ ( block_left + k + 1 ) * 2 ] - counts_at_start;
                        }

                        for( ; k < ( block_span + 8 & -8LL ); k++ ) {
                            assert( ((uint16_t *)s)[ k + 8 ] == 0 );
                            ((uint16_t *)s)[ k + 8 ] = 0xFFFFU;
                        }

                        assert( block_span / 8 <= 8 );

                        for( k = 0; k < block_span / 8; k++ ) {
                            assert( ((uint16_t *)s)[ k ] == 0 );
                            ((uint16_t *)s)[ k ] = counts[ ( block_left + ( k + 1 ) * 8 ) * 2 ] - counts_at_start;
                        }

                        for( ; k < 8; k++ ) {
                            assert( ((uint16_t *)s)[ k ] == 0 );
                            ((uint16_t *)s)[ k ] = 0xFFFFU;
                        }
                    }
                    else if ( span >= 2 ) {
                        assert( ( block_span + 8 & -8LL ) <= span * 4 );

                        int k;
                        for( k = 0; k < block_span; k++ ) {
                            assert( ((uint16_t *)s)[ k ] == 0 );
                            ((uint16_t *)s)[ k ] = counts[ ( block_left + k + 1 ) * 2 ] - counts_at_start;
                        }

                        for( ; k < ( block_span + 8 & -8LL ); k++ ) {
                            assert( ((uint16_t *)s)[ k ] == 0 );
                            ((uint16_t *)s)[ k ] = 0xFFFFU;
                        }
                    }
                }

                switch( state ) {
                    case 0:
                        assert( s[ d & INVENTORY_MASK ] == 0 );
                        s[ d & INVENTORY_MASK ] = i * 64 + j;
                        break;
                    case 1:
                        assert( ((uint32_t *)s)[ d & INVENTORY_MASK ] == 0 );
                        assert( i * 64 + j - first_bit < (1ULL << 32) );
                        ((uint32_t *)s)[ d & INVENTORY_MASK ] = i * 64 + j - first_bit;
                        break;
                    case 2:
                        assert( ((uint16_t *)s)[ d & INVENTORY_MASK ] == 0 );
                        assert( i * 64 + j - first_bit < (1 << 16) );
                        ((uint16_t *)s)[ d & INVENTORY_MASK ] = i * 64 + j - first_bit;
                        break;
                }

                d++;
            }
        }
    }
}

template<uint8_t b, uint8_t pattern_len>
inline const typename select_support_v9<b,pattern_len>::size_type select_support_v9<b,pattern_len>::select(size_type rank)const
{
    rank--; // vigna starts with 0. we start with 1

    const uint64_t* subinventory = m_subinventory.data();
    const uint64_t* counts = m_counts.data();
    const uint64_t* inventory = m_inventory.data();
    const uint64_t* bits = m_v->data();

    const uint64_t inventory_index_left = rank >> LOG2_ONES_PER_INVENTORY;
    assert( inventory_index_left < inventory_size );

    const uint64_t inventory_left = inventory[ inventory_index_left ];
    const uint64_t block_right = inventory[ inventory_index_left + 1 ] / 64;
    uint64_t block_left = inventory_left / 64;
    const uint64_t span = block_right / 4 - block_left / 4;
    const uint64_t * const s = &subinventory[ block_left / 4 ];
    uint64_t count_left, rank_in_block;

#ifdef DEBUG
    printf( "Initially, rank: %lld block_left: %lld block_right: %lld span: %lld\n", rank, block_left, block_right, span );
#endif

    if ( span < 2 ) {
        block_left &= ~7;
        count_left = block_left / 4 & ~1;
        assert( rank < counts[ count_left + 2 ] );
        rank_in_block = rank - counts[ count_left ];
#ifdef DEBUG
        printf( "Single span; rank_in_block: %lld block_left: %lld\n", rank_in_block, block_left );
#endif
#ifdef COUNTS
        single++;
#endif
    }
    else if ( span < 16 ) {
        block_left &= ~7;
        count_left = block_left / 4 & ~1;
        const uint64_t rank_in_superblock = rank - counts[ count_left ];
        const uint64_t rank_in_superblock_step_16 = rank_in_superblock * ONES_STEP_16;

        const uint64_t first = s[ 0 ], second = s[ 1 ];
        const int where = ( ULEQ_STEP_16( first, rank_in_superblock_step_16 ) + ULEQ_STEP_16( second, rank_in_superblock_step_16 ) ) * ONES_STEP_16 >> 47;
        assert( where >= 0 );
        assert( where <= 16 );
#ifdef DEBUG
        printf( "rank_in_superblock: %lld (%llx) %llx %llx\n", rank_in_superblock, rank_in_superblock, s[0], s[1] );
#endif

        block_left += where * 4;
        count_left += where;
        rank_in_block = rank - counts[ count_left ];
        assert( rank_in_block >= 0 );
        assert( rank_in_block < 512 );
#ifdef DEBUG
        printf( "Found where (1): %d rank_in_block: %lld block_left: %lld\n", where, rank_in_block, block_left );
        printf( "supercounts: %016llx %016llx\n", s[0], s[1] );
#endif
#ifdef COUNTS
        one_level++;
#endif
    }
    else if ( span < 128 ) {
        block_left &= ~7;
        count_left = block_left / 4 & ~1;
        const uint64_t rank_in_superblock = rank - counts[ count_left ];
        const uint64_t rank_in_superblock_step_16 = rank_in_superblock * ONES_STEP_16;

        const uint64_t first = s[ 0 ], second = s[ 1 ];
        const int where0 = ( ULEQ_STEP_16( first, rank_in_superblock_step_16 ) + ULEQ_STEP_16( second, rank_in_superblock_step_16 ) ) * ONES_STEP_16 >> 47;
        assert( where0 <= 16 );
        const uint64_t first_bis = s[ where0 + 2 ], second_bis = s[ where0 + 2 + 1 ];
        const int where1 = where0 * 8 + ( ( ULEQ_STEP_16( first_bis, rank_in_superblock_step_16 ) + ULEQ_STEP_16( second_bis, rank_in_superblock_step_16 ) ) * ONES_STEP_16 >> 47 );

        block_left += where1 * 4;
        count_left += where1;
        rank_in_block = rank - counts[ count_left ];
        assert( rank_in_block >= 0 );
        assert( rank_in_block < 512 );

#ifdef DEBUG
        printf( "Found where (2): %d rank_in_block: %lld block_left: %lld\n", where1, rank_in_block, block_left );
#endif
#ifdef COUNTS
        two_levels++;
#endif
    }
    else if ( span < 256 ) {
#ifdef COUNTS
        shorts++;
#endif
        return ((uint16_t*)s)[ rank % ONES_PER_INVENTORY ] + inventory_left;
    }
    else if ( span < 512 ) {
#ifdef COUNTS
        longs++;
#endif
        return ((uint32_t*)s)[ rank % ONES_PER_INVENTORY ] + inventory_left;
    }
    else {
#ifdef COUNTS
        longlongs++;
#endif
        return s[ rank % ONES_PER_INVENTORY ];
    }

    const uint64_t rank_in_block_step_9 = rank_in_block * ONES_STEP_9;
    const uint64_t subcounts = counts[ count_left + 1 ];
    const uint64_t offset_in_block = ( ULEQ_STEP_9( subcounts, rank_in_block_step_9 ) * ONES_STEP_9 >> 54 & 0x7 );

    const uint64_t word = block_left + offset_in_block;
    const uint64_t rank_in_word = rank_in_block - ( subcounts >> ( offset_in_block - 1 & 7 ) * 9 & 0x1FF );
#ifdef DEBUG
    printf( "rank_in_block: %lld offset_in_block: %lld rank_in_word: %lld compare: %016llx shift: %lld\n", rank_in_block, offset_in_block, rank_in_word, UCOMPARE_STEP_9( rank_in_block_step_9, subcounts ), subcounts >> ( offset_in_block - 1 ) * 9 & 0x1FF );
#endif
    assert( offset_in_block >= 0 );
    assert( offset_in_block <= 7 );

#ifdef DEBUG
    printf( "rank_in_block: %lld offset_in_block: %lld rank_in_word: %lld\n", rank_in_block, offset_in_block, rank_in_word );
    printf( "subcounts: " );
    for( int i = 0; i < 7; i++ ) printf( "%lld ", subcounts >> i * 9 & 0x1FF );
    printf( "\n" );
    fflush( stdout );
#endif

    assert( rank_in_word < 64 );
    assert( rank_in_word >= 0 );

#ifdef DEBUG
    printf( "rank_in_word: %ld" , rank_in_word ); fflush( stdout );
    printf( "Returning %ld %lld\n", rank_in_word, word * 64ULL + select_in_word( bits[ word ], rank_in_word ) );
#endif
    return word * 64ULL + select_in_word( bits[ word ], rank_in_word );
}

template<uint8_t b, uint8_t pattern_len>
inline const typename select_support_v9<b,pattern_len>::size_type select_support_v9<b,pattern_len>::operator()(size_type i)const
{
    return select(i);
}

template<uint8_t b, uint8_t pattern_len>
void select_support_v9<b,pattern_len>::set_vector(const int_vector<1>* v)
{
    m_v = v;
}

template<uint8_t b, uint8_t pattern_len>
typename select_support_v9<b,pattern_len>::size_type select_support_v9<b,pattern_len>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;

    written_bytes += m_inventory.serialize(out,child,"inventory");
    written_bytes += m_subinventory.serialize(out,child,"subinventory");
    written_bytes += m_counts.serialize(out,child,"counts");

    written_bytes += util::write_member(log2_ones_per_inventory, out, child, "log2_ones_per_inventory");
    written_bytes += util::write_member(inventory_size, out, child, "inventory_size");
    written_bytes += util::write_member(ones_per_inventory, out, child, "ones_per_inventory");
    written_bytes += util::write_member(num_counts, out, child, "num_counts");
    written_bytes += util::write_member(num_ones, out, child, "num_ones");
    written_bytes += util::write_member(num_words, out, child, "num_words");

    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<uint8_t b, uint8_t pattern_len>
void select_support_v9<b,pattern_len>::load(std::istream& in, const int_vector<1>* v)
{
    set_vector(v);

    m_inventory.load(in);
    m_subinventory.load(in);
    m_counts.load(in);

    util::read_member(log2_ones_per_inventory, in);
    util::read_member(inventory_size, in);
    util::read_member(ones_per_inventory, in);
    util::read_member(num_counts, in);
    util::read_member(num_ones, in);
    util::read_member(num_words, in);
}



}

#endif
