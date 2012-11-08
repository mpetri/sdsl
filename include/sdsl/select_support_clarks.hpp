/*! \file select_support_clarks.hpp
    \brief select_support_clarks.hpp contains classes that support a sdsl::bit_vector with constant time select information.
	\author Matthias Petri
*/
#ifndef INCLUDED_SDSL_select_support_clarks
#define INCLUDED_SDSL_select_support_clarks

#include "int_vector.hpp"
#include "select_support.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{


template<uint8_t bit_pattern, uint8_t pattern_len>
struct select_support_clarks_trait {
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
struct select_support_clarks_trait<0,1> {
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
struct select_support_clarks_trait<1,1> {
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
struct select_support_clarks_trait<10,2> {
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
struct select_support_clarks_trait<01,2> {
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


//! A class containing a faithful implementation of constant time
//  select queries proposed by clark in 1996.
//
//  The space usage is at most 3n/log(log(n)).
//
//   First the bitvector is partitioned into superblocks covering
//   log(n)*log(log(n)) elements each. Depending on the size r
//   of the superblock more positions are stored:
//
//    Case 1) r > log(n)^2* log(log(n))^2: all positions are stored
//              explicitly using log(n) bits each. this block
//              is considered a "long" block.
//
//    Case 2) r < log(n)^2* log(log(n))^2: subdivide each superblock
//              into subblocks covering log(r)*log(log(n)) positions
//              each. The size of each subblock is r'. We store the
//              first position of each block relative to the superblock
//              using log(r) bits.
//
//      Case 2a) If r' > log(r')*log(r)*log(log(n)) bits we additional
//              store the position of each individual position relative
//              to the position of each subblock at a cost of log(r') bits.
//              this is referred to as a "miniblock". A superblock
//              can have multiple subblocks. Each subblock can correspond
//              to one miniblock.
//      Case 2b) If r' < log(r')*log(r)*log(log(n)) bits we do not store
//              miniblocks as the range is too small. Instead, during query
//              time we scan a "constant" portion of the bitvector to
//              retrieve the final position.
//
//
//  The space usage is at most 3n/log(log(n)) in the case we store subblocks
//  and miniblocks for all superblocks.
//

// temp structure to process a superblock
typedef struct {
    size_t id;
    size_t r;
    size_t logr;
    int_vector<64>  elements;
    int_vector<0>   longblock;
    int_vector<0>   blocks;
    int_vector<64>  block_sizes;
    bit_vector      miniblocks;
    size_t num_miniblocks;
} superblock_t;

template<uint8_t b=1, uint8_t pattern_len=1>
class select_support_clarks : public select_support
{
    public:
        typedef bit_vector bit_vector_type;
    private:
        size_type m_count;                  // number of elements
        size_type m_numsuperblocks;         // number of superblocks
        uint32_t  m_elems_per_superblock;   // log(n)*log(log(n))
        uint32_t  m_longthreshold;          // log^2(n)*(log(log(n)))^2
        uint32_t  m_logn;                   // log(n)
        uint32_t  m_loglogn;                // log(log(n))
        int_vector<0>  m_superblocks;       // log(n) for each position
        int_vector<0>  m_longblocks;        // log(n) for each position
        int_vector<0>  m_block_pos;         // start positions of each block
        bit_vector     m_blocks_and_mini;   // store blocks and mini blocks
        size_type      m_size_miniblocks;
        size_type      m_size_blocks;
        size_type      m_num_long_elems;
    private:
        void copy(const select_support_clarks<b, pattern_len>& ss);
        void construct();
        void initData();

        void store_longblock(superblock_t& sb);
        void store_blocks(superblock_t& sb);
        void store_miniblock(superblock_t& sb,size_type block,int_vector<0>*& mini_blocks,size_t width);
        void process_superblock(superblock_t& sb,bool last,size_t& start_long,size_t& start_block);
        void compress_mini_blocks(superblock_t& sb,int_vector<0>*& miniblocks);
        void compress_superblock(superblock_t& sb,bool last,size_t& start_long,size_t& start_block);
        void check_compressed_superblock(superblock_t& sb,bool last);
    public:
        explicit select_support_clarks(const int_vector<1>* v=NULL);
        select_support_clarks(const select_support_clarks<b,pattern_len>& ss);
        void init(const int_vector<1>* v=NULL);

        void block_bytes(size_t& sb,size_t& lb,size_t& bb,size_t& mb) const;

        inline const size_type select(size_type i) const;
        //! Alias for select(i).
        inline const size_type operator()(size_type i)const;

        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const;
        void load(std::istream& in, const int_vector<1>* v=NULL);
        void set_vector(const int_vector<1>* v=NULL);
        select_support_clarks<b, pattern_len>& operator=(const select_support_clarks& ss);

        //! Swap operator
        /*! This swap operator swaps two select_support_clarkss in constant time.
         */
        void swap(select_support_clarks<b, pattern_len>& ss);

};


template<uint8_t b, uint8_t pattern_len>
void select_support_clarks<b,pattern_len>::construct()
{
    m_count	= 0;
    m_numsuperblocks = 0;
}

template<uint8_t b, uint8_t pattern_len>
select_support_clarks<b,pattern_len>::select_support_clarks(const int_vector<1>* f_v):select_support(f_v)
{
    construct();
    init(f_v);
}

template<uint8_t b, uint8_t pattern_len>
select_support_clarks<b,pattern_len>::select_support_clarks(const select_support_clarks& ss):select_support(ss.m_v)
{
    init();
    copy(ss);
}

template<uint8_t b, uint8_t pattern_len>
select_support_clarks<b, pattern_len>& select_support_clarks<b,pattern_len>::operator=(const select_support_clarks& ss)
{
    if (this != &ss) {
        copy(ss);
    }
    return *this;
}

template<uint8_t b, uint8_t pattern_len>
void select_support_clarks<b,pattern_len>::swap(select_support_clarks& ss)
{
    std::swap(m_logn, ss.m_logn);
    std::swap(m_loglogn, ss.m_loglogn);
    std::swap(m_longthreshold, ss.m_longthreshold);
    std::swap(m_elems_per_superblock, ss.m_elems_per_superblock);
    std::swap(m_count, ss.m_count);
    std::swap(m_numsuperblocks, ss.m_numsuperblocks);
    std::swap(m_size_blocks, ss.m_size_blocks);
    std::swap(m_size_miniblocks, ss.m_size_miniblocks);
    std::swap(m_num_long_elems, ss.m_num_long_elems);
    m_superblocks.swap(ss.m_superblocks);
    m_longblocks.swap(ss.m_longblocks);
    m_block_pos.swap(ss.m_block_pos);
    m_blocks_and_mini.swap(ss.m_blocks_and_mini);
}

template<uint8_t b, uint8_t pattern_len>
void select_support_clarks<b,pattern_len>::copy(const select_support_clarks<b, pattern_len>& ss)
{
    m_logn			 = ss.m_logn;
    m_loglogn        = ss.m_loglogn;
    m_longthreshold  = ss.m_longthreshold;
    m_elems_per_superblock = ss.m_elems_per_superblock;
    m_superblocks    = ss.m_superblocks;
    m_count          = ss.m_count;
    m_numsuperblocks = ss.m_numsuperblocks;
    m_superblocks    = ss.m_superblocks;
    m_longblocks    = ss.m_longblocks;
    m_block_pos    = ss.m_block_pos;
    m_blocks_and_mini    = ss.m_blocks_and_mini;
    m_size_blocks  = ss.m_size_blocks;
    m_size_miniblocks = ss.m_size_miniblocks;
    m_num_long_elems = ss.m_num_long_elems;
}

/**
 *  store all positions in a superblock explicitly using log(n) bits each.
 *  we create a "long" block.
 */
template<uint8_t b, uint8_t pattern_len>
void select_support_clarks<b,pattern_len>::store_longblock(superblock_t& sb)
{
    // we store all positions explicitly using log(n) bits
    sb.longblock.set_int_width(m_logn);
    sb.longblock.resize(sb.elements.size()-1);

    // copy the positions into the long block
    for (size_type j=0; j < sb.elements.size()-1; j++) {
        sb.longblock[j] = sb.elements[j+1];
    }
}

/**
 *   we subdivide a superblock and store position samples using log(r) bits each.
 */
template<uint8_t b, uint8_t pattern_len>
void select_support_clarks<b,pattern_len>::store_blocks(superblock_t& sb)
{
    // we don't count the first position as we store it explicitly
    size_type subblock_size = sb.logr * m_loglogn;
    size_type num_blocks = sb.elements.size()/subblock_size;

    //fprintf(stderr, "===========================================\n");
    //fprintf(stderr, "num blocks = %zu\n",num_blocks);
    //fprintf(stderr, "subblock_size = %zu\n",subblock_size);
    //fprintf(stderr, "sb.elements.size() = %zu\n",sb.elements.size());

    if (sb.elements.size() % subblock_size == 0) num_blocks--;


    if (num_blocks) {
        sb.blocks.set_int_width(sb.logr);
        sb.blocks.resize(num_blocks);
        sb.block_sizes.resize(num_blocks+1); // we have one more block
        // than we actually store

        // we sample every blocksize-th position and store the relative
        // position in log(r) bits
        size_type start = sb.elements[0];
        size_type prev = 0;
        for (size_type i = 0; i < num_blocks; ++i) {
            sb.blocks[i] = sb.elements[(i+1)*subblock_size]-start;
            sb.block_sizes[i] = sb.blocks[i]-prev;
            prev = sb.blocks[i];
        }
        // last block
        sb.block_sizes[num_blocks] = sb.r - prev;

        // for each block we might have to subdivide again
        int_vector<0>* mini_blocks = new int_vector<0>[sb.block_sizes.size()];
        sb.num_miniblocks = 0;
        for (size_type i = 0; i <  sb.block_sizes.size(); ++i) {
            size_type logrprime = bit_magic::l1BP(sb.block_sizes[i])+1; // log(r')
            size_type threshold = logrprime*sb.logr*m_loglogn*m_loglogn;
            //fprintf(stderr, "sb %zu block %zu block size = %zu threshold = %zu\n",sb.id,i,sb.block_sizes[i],threshold);
            if (sb.block_sizes[i] >= threshold) {
                sb.num_miniblocks++;
                store_miniblock(sb,i,mini_blocks,logrprime);
            }
        }

        // output sb
        /*
        std::cerr << "sb " << sb.id << std::endl;
        std::cerr << "r " << sb.r << std::endl;
        std::cerr << "logr " << sb.logr << std::endl;
        std::cerr << "subblock_size " << subblock_size << std::endl;
        std::cerr << "num_blocks " << num_blocks << std::endl;
        std::cerr << "num_miniblocks " << sb.num_miniblocks << std::endl;
        for (size_type i = 0; i < num_blocks; ++i) {
            std::cerr << "block[" << i << "] = " << sb.blocks[i] << std::endl;
        }
        for (size_type i = 0; i <  sb.block_sizes.size(); ++i) {
            size_type logrprime = bit_magic::l1BP(sb.block_sizes[i])+1; // log(r')
            size_type threshold = logrprime*sb.logr*m_loglogn*m_loglogn;
            //size_type threshold = m_loglogn*m_loglogn;
            std::cerr << "rprime " << sb.block_sizes[i] << std::endl;
            std::cerr << "logrprime " << logrprime << std::endl;
            std::cerr << "threshold " << threshold << std::endl;

            if( sb.block_sizes[i] >= threshold) {
                for (size_type j = 0; j < mini_blocks[i].size(); ++j) {
                    std::cerr << "block[" << i << "]-mini[" << j << "] = " << mini_blocks[i][j] << std::endl;
                }
            }
        }*/

        // if we created miniblocks compress them to save space
        if (sb.num_miniblocks) compress_mini_blocks(sb,mini_blocks);
        delete [] mini_blocks;
    }
}

/**
 *   store all positions in a given block explicitly relative to the subblock
 *   using log(r') bits each.
 */
template<uint8_t b, uint8_t pattern_len>
void
select_support_clarks<b,pattern_len>::store_miniblock(superblock_t& sb,size_type block,
        int_vector<0>*& mini_blocks,
        size_t width)
{
    // calc number of elements to store
    size_type num_elems = sb.logr*m_loglogn;
    if (block == sb.block_sizes.size()-1) { // last?
        num_elems = sb.elements.size() % (sb.logr*m_loglogn);
        if (num_elems == 0) num_elems = sb.logr*m_loglogn;
    }
    //fprintf(stderr, "STORE MINI: NUM ELEMS = %zu\n",num_elems);

    //fprintf(stderr, "num elements in miniblock of block %zu = %zu\n",block,num_elems);

    // store each relative position in log(r') bits
    mini_blocks[block] = int_vector<0>(num_elems-1 , 0 , width);

    // we store relative to the samples block pos
    size_type start = sb.elements[0];
    if (block!=0)
        start += sb.blocks[block-1];

    size_type skip = block*sb.logr*m_loglogn;
    for (size_type i = 1; i < num_elems; ++i) {
        mini_blocks[block][i-1] = sb.elements[skip+i] - start;
    }
}


/**
 *   process a superblock. first decide if we have a long block or if we
 *   divide into subblocks.
 */
template<uint8_t b, uint8_t pattern_len>
void select_support_clarks<b,pattern_len>::process_superblock(superblock_t& sb,bool last,
        size_t& start_long,size_t& start_block)
{
    // by convention we always store the last sb as a long block
    if (last || sb.r >= m_longthreshold) {
        store_longblock(sb);
    } else {
        // subdivide the super block
        store_blocks(sb);
    }
    compress_superblock(sb,last,start_long,start_block);
    check_compressed_superblock(sb,last);
}

/**
 *   process a superblock. first decide if we have a long block or if we
 *   divide into subblocks.
 */
template<uint8_t b, uint8_t pattern_len>
void select_support_clarks<b,pattern_len>::check_compressed_superblock(superblock_t& sb,bool last)
{
    // check that everything is stored in the right place
    if (last || sb.r >= m_longthreshold) {
        //////fprintf(stderr, "start longblocks = %zu\n",m_block_pos[sb.id]);
        for (size_type i=0; i<sb.longblock.size(); i++) {
            if (sb.longblock[i] != m_longblocks[ m_block_pos[sb.id] + i ]) {
                //fprintf(stderr, "[%zu] ERROR! long block value %zu (%u <-> %u)\n",sb.id,i,sb.longblock[i],m_longblocks[ m_block_pos[sb.id] + i ]);
                exit(EXIT_FAILURE);
            }
        }
    } else {

        size_type start = m_block_pos[sb.id];
        for (size_type i=0; i<sb.blocks.size(); i++) {
            uint64_t* block_ptr = (uint64_t*)(m_blocks_and_mini.data()+ (start>>6));
            uint8_t offset = start%64;
            uint64_t val = bit_magic::read_int(block_ptr,offset,sb.logr);
            if (sb.blocks[i] != val) {
                //fprintf(stderr, "[%zu] ERROR! block value %zu (%u <-> %u)\n",sb.id,i,sb.blocks[i],val);
                exit(EXIT_FAILURE);
            }
            start += sb.logr;
        }

        // check the miniblocks
        if (sb.num_miniblocks) {
            for (size_type i=0; i<sb.miniblocks.size(); i++) {
                if (sb.miniblocks[i] != m_blocks_and_mini[start+i]) {
                    //fprintf(stderr, "[%zu] ERROR! miniblock value %zu\n",sb.id,i);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
}

/**
 *   process a superblock. first decide if we have a long block or if we
 *   divide into subblocks.
 */
template<uint8_t b, uint8_t pattern_len>
void select_support_clarks<b,pattern_len>::compress_superblock(superblock_t& sb,bool last,
        size_t& start_long,size_t& start_block)
{
    // by convention we always store the last sb as a long block
    if (last || sb.r >= m_longthreshold) {
        // compress the long block
        m_block_pos[sb.id] = start_long;
        for (size_type i=0; i<sb.longblock.size(); i++) {
            if (m_longblocks.size() <= start_long+i) {
                m_longblocks.resize(m_longblocks.size()*2);
            }
            m_longblocks[start_long+i] = sb.longblock[i];
        }
        // move the start of the next long block we save
        start_long += sb.longblock.size();
        m_num_long_elems += sb.longblock.size();
    } else {
        // compress the block data

        // check if we have enough data
        if (m_blocks_and_mini.size() <= (start_block+sb.logr*sb.blocks.size()+sb.miniblocks.size())) {
            m_blocks_and_mini.resize(sb.logr*sb.blocks.size()+
                                     sb.miniblocks.size()+m_blocks_and_mini.size()*2);
        }

        // first write the blocks
        m_block_pos[sb.id] = start_block;
        uint64_t* block_ptr = (uint64_t*)(m_blocks_and_mini.data()+(start_block>>6));
        uint8_t offset = start_block%64;

        for (size_type i=0; i<sb.blocks.size(); i++) {
            bit_magic::write_int_and_move(block_ptr,sb.blocks[i],offset,sb.logr);
        }
        start_block += sb.logr*sb.blocks.size();
        m_size_blocks += sb.logr*sb.blocks.size();

        // write the miniblocks if there are any
        if (sb.num_miniblocks) {
            for (size_type i=0; i<sb.miniblocks.size(); i++) {
                m_blocks_and_mini[start_block+i] = sb.miniblocks[i];
            }
            start_block += sb.miniblocks.size();
            m_size_miniblocks += sb.miniblocks.size();
        }
    }
}

/**
 *   compress all miniblocks in a superblock into one bitvector.
 *   the layout is as follows:
 *
 *   [mask][offsetwidth][offsetpos][offsetpos][    data          ]
 *
 *   where mask is a 64bit integer marking the subblocks in the superblock
 *   we store explicitly.
 *
 *   offsetwidth is the size of the integers storing the starting positions
 *   of each stored miniblock in the bitvector
 *
 *   each offsetpos stored the starting position of a stored miniblock in
 *   the bitvector.
 */
template<uint8_t b, uint8_t pattern_len>
void
select_support_clarks<b,pattern_len>::compress_mini_blocks(superblock_t& sb,int_vector<0>*& miniblocks)
{

    uint64_t mask = 0;
    uint8_t stored = 0;
    // check if we have to do something first and count the size of the
    // resulting vector
    size_type bits = 0;
    for (size_type i = 0; i < sb.block_sizes.size(); ++i) {
        if (miniblocks[i].size() > 0) {
            bits += miniblocks[i].size()*miniblocks[i].get_int_width();
            mask |= (1ULL << i);
            stored++;
        }
    }

    if (stored) {
        size_type header_size = sizeof(mask) + sizeof(stored) + stored*32;
        size_type bits_per_pos = bit_magic::l1BP(header_size+bits)+1;

        // calc the total bitvector size which includes the header
        bits += sizeof(mask)*8;      // mask
        bits += sizeof(stored)*8;    // how many blocks are stored
        bits += sizeof(stored)*8;    // the size of each starting pos
        bits += stored*bits_per_pos; // starting pos of each mini block


        sb.miniblocks.resize(bits);

        // write
        uint64_t* data = (uint64_t*) sb.miniblocks.data();
        uint8_t offset = 0;
        size_type written = 0;

        // write header
        bit_magic::write_int_and_move(data,mask,offset,sizeof(mask)*8);
        written += sizeof(mask)*8;
        bit_magic::write_int_and_move(data,stored,offset,sizeof(stored)*8);
        written += sizeof(stored)*8;
        bit_magic::write_int_and_move(data,bits_per_pos,offset,sizeof(stored)*8);
        written += sizeof(stored)*8;

        size_type cur_offset = sizeof(mask)*8 + sizeof(stored)*8 + sizeof(stored)*8
                               + stored*bits_per_pos;

        //fprintf(stderr, "mask = %lu\n",mask);
        //fprintf(stderr, "offsetwidth = %lu\n",bits_per_pos);

        // write length
        for (size_type i = 0; i < sb.block_sizes.size(); ++i) {
            if (miniblocks[i].size() > 0) {
                bit_magic::write_int_and_move(data,cur_offset,offset,bits_per_pos);
                ////fprintf(stderr, "start_pos = %lu\n",cur_offset);
                cur_offset += miniblocks[i].size()*miniblocks[i].get_int_width();
                written += bits_per_pos;
            }
        }

        // write data
        for (size_type i = 0; i < sb.block_sizes.size(); ++i) {
            if (miniblocks[i].size() > 0) {
                for (size_type j = 0; j < miniblocks[i].size(); ++j) {
                    bit_magic::write_int_and_move(data,miniblocks[i][j],offset,
                                                  miniblocks[i].get_int_width());
                    written += miniblocks[i].get_int_width();
                }
            }
        }
        assert(written == bits);
    }
}


template<uint8_t b, uint8_t pattern_len>
void select_support_clarks<b,pattern_len>::init(const int_vector<1>* v)
{
    set_vector(v);
    initData();
    if (m_v==NULL)
        return;

    // Count the number of elements in the bit vector
    m_count = select_support_clarks_trait<b,pattern_len>::arg_cnt(*v);
    if (m_count==0) // if there are no elements in the vector we are done...
        return;

    m_numsuperblocks = (m_count+(m_elems_per_superblock-1))/m_elems_per_superblock;
    m_superblocks = int_vector<0>(m_numsuperblocks+1, 0, m_logn);
    m_block_pos = int_vector<0>(m_numsuperblocks,0,64); // we bitcompress later

    // resize those later
    m_longblocks = int_vector<0>(4096, 0, m_logn);
    m_blocks_and_mini.resize(4096*m_logn);

    size_type cnt = 0;
    size_type cur_block_cnt = 0;
    size_type cur_superblock = 0;

    //fprintf(stderr, "CREATING CLARK FOR:\n");
    //fprintf(stderr, "BV SIZE %zu\n",v->size());
    //fprintf(stderr, "num superblocks %zu\n",m_numsuperblocks);

    /* find the first elem for the first superblock */
    superblock_t sb;
    sb.id = 0;
    sb.elements.resize(m_elems_per_superblock);
    for (size_type i=0; i < v->size(); ++i) {
        if (select_support_clarks_trait<b,pattern_len>::found_arg(i, *v)) {
            m_superblocks[cur_superblock] = i;
            sb.elements[0] = i;
            cur_block_cnt = 1;
            cnt = 1;
            break;
        }
    }

    // after we found the first element of the next superblock we can process
    // the previous superblock.
    size_type start_long = 0;
    size_type start_block = 0;
    for (size_type i=m_superblocks[0]+1; i < v->size(); ++i) {
        if (select_support_clarks_trait<b,pattern_len>::found_arg(i, *v)) {

            // have we found the beginning of the next superblock
            if (cur_block_cnt == m_elems_per_superblock) {

                // calc r and log(r)
                sb.r = i - sb.elements[0];
                sb.logr = bit_magic::l1BP(sb.r)+1;
                process_superblock(sb,false,start_long,start_block);

                // start the next superblock
                cur_superblock++;
                m_superblocks[cur_superblock] = i;
                sb.id = cur_superblock;
                cur_block_cnt = 0;
            }

            // add new element to the superblock
            sb.elements[cur_block_cnt] = i;
            cnt++;
            cur_block_cnt++;
            if (cnt == m_count) break; // last element found
        }
    }
    // process the last superblock if there are elements left over
    if (cur_block_cnt!=0) {
        sb.elements.resize(cur_block_cnt);
        sb.r = v->size() - sb.elements[0];
        sb.logr = bit_magic::l1BP(sb.r)+1;
        process_superblock(sb,true,start_long,start_block);

        // we explicitly store this so we can determine the size of the last
        // superblock faster during select()
        m_superblocks[m_numsuperblocks] = v->size();
    }

    // bit compress the values if we can
    util::bit_compress(m_block_pos);
    util::bit_compress(m_superblocks);
    m_longblocks.resize(m_num_long_elems);
    util::bit_compress(m_longblocks);
    m_blocks_and_mini.resize(m_size_blocks+m_size_miniblocks);

    assert(m_numsuperblocks == cur_superblock);
}

template<uint8_t b, uint8_t pattern_len>
inline const typename select_support_clarks<b,pattern_len>::size_type select_support_clarks<b,pattern_len>::select(size_type i) const
{
    assert(i > 0);
    assert(i <= m_count);

    // numbering starts at 0 in the data structure
    i--;

    //fprintf(stderr, "=================== SELECT(%zu)\n",i);

    /* determine superblock */
    size_type superblock = i/m_elems_per_superblock;
    size_type sbv = m_superblocks[superblock];
    size_type sbn = m_superblocks[superblock+1];
    size_type block_offset = i%m_elems_per_superblock;
    /* calc superblock size */
    size_type r = sbn - sbv;
    //fprintf(stderr, "superblock = %zu\n",superblock);
    //fprintf(stderr, "sbv = %zu\n",sbv);
    //fprintf(stderr, "sbn = %zu\n",sbn);
    //fprintf(stderr, "block_offset = %zu\n",block_offset);
    //fprintf(stderr, "r = %zu\n",r);

    if (block_offset==0) { // we store the position explicitly as it is the
        // first element of a superblock
        return sbv;
    }
    if (superblock == m_numsuperblocks-1 || r >= m_longthreshold) {
        // we are inside a long block. get the answer directly from
        // the pre-stored positions in the long block
        return m_longblocks[ m_block_pos[superblock] + block_offset-1];
    } else {
        // we are in a subdivided superblock. find the correct block
        size_type logr = bit_magic::l1BP(r)+1;
        size_type elements_per_block = logr*m_loglogn;
        size_type block = block_offset/elements_per_block;

        //fprintf(stderr, "logr = %zu\n",logr);
        //fprintf(stderr, "elements_per_block = %zu\n",elements_per_block);
        //fprintf(stderr, "block = %zu\n",block);

        // get the block value
        size_type block_store_offset = m_block_pos[superblock];
        size_type block_value = 0;
        if (block) {
            block_store_offset += (block-1)*logr;
            const uint64_t* data = (const uint64_t*)(m_blocks_and_mini.data()+(block_store_offset>>6));
            uint8_t offset = block_store_offset%64;
            block_value = bit_magic::read_int(data,offset,logr);
        }

        // check if the answer is one of the ones stored in the blocks directly
        if (block_offset % (elements_per_block) == 0) {
            // if this is true, block always >= 1
            assert(block >= 1);
            size_type block_store_offset = m_block_pos[superblock];
            block_store_offset += (block-1)*logr;
            const uint64_t* data = (const uint64_t*)(m_blocks_and_mini.data()+(block_store_offset>>6));
            uint8_t offset = block_store_offset%64;
            size_type block_prev = bit_magic::read_int(data,offset,logr);
            return sbv + block_prev;
        }

        // determine the size of the block
        size_type rprime = block_value; // if we are in the first block
        //fprintf(stderr, "block_value = %zu\n",block_value);
        if (block!=0) {
            if (m_elems_per_superblock <= elements_per_block*(block+1)) {  // last block
                rprime = sbn-(sbv+block_value);
            } else {
                block_store_offset += logr;
                const uint64_t* data = (const uint64_t*)(m_blocks_and_mini.data()+(block_store_offset>>6));
                uint8_t offset = block_store_offset%64;
                size_type block_next = bit_magic::read_int(data,offset,logr);
                rprime = block_next - block_value;
            }
        }

        size_type logrprime = bit_magic::l1BP(rprime)+1;
        size_type threshold = logrprime*logr*m_loglogn*m_loglogn;
        //fprintf(stderr, "rprime = %zu\n",rprime);
        //fprintf(stderr, "logrprime = %zu\n",logrprime);
        if (rprime >= threshold) {
            // we have miniblocks.
            size_type num_blocks = m_elems_per_superblock/elements_per_block;
            if (m_elems_per_superblock%elements_per_block==0) num_blocks--;

            //fprintf(stderr, "num_blocks = %zu\n",num_blocks);


            size_type block_store_offset = m_block_pos[superblock];
            block_store_offset += num_blocks*logr; // skip overall block values
            const uint64_t* data = (const uint64_t*)(m_blocks_and_mini.data()+(block_store_offset>>6));
            uint8_t offset = block_store_offset%64;

            // determine which of the stored mini blocks we want to access
            uint64_t mask = bit_magic::read_int_and_move(data,offset,64);
            size_type id = bit_magic::b1Cnt(mask&bit_magic::Li1Mask[block+1]);

            //fprintf(stderr, "mask = %zu\n",mask);
            //fprintf(stderr, "masked = %zu\n",mask&bit_magic::Li1Mask[block+1]);

            // get the starting pos of that block
            uint64_t offsetwidth = bit_magic::read_int_and_move(data,offset,8);
            offsetwidth = bit_magic::read_int_and_move(data,offset,8);

            //fprintf(stderr, "offsetwidth = %zu\n",offsetwidth);

            size_t start_offset = offset+((id-1)*offsetwidth);
            data = (const uint64_t*)(data+(start_offset>>6));
            offset = start_offset&63;
            uint64_t start_pos = bit_magic::read_int(data,offset,offsetwidth);

            // skip to the correct element
            size_type subblock_offset = block_offset % elements_per_block;
            start_pos += ((subblock_offset-1)*logrprime);

            // now retrieve the number we are looking for
            block_store_offset += start_pos;
            data = (const uint64_t*)(m_blocks_and_mini.data()+(block_store_offset>>6));
            offset = block_store_offset&63;
            uint64_t elem_pos = bit_magic::read_int(data,offset,logrprime);

            // finally add up all the numbers
            if (block==0) return sbv + elem_pos;
            else return sbv + block_value + elem_pos;
        } else {
            // we scan
            size_type finalpos = sbv;
            // perform a scan through the bitvector

            // determine the starting position in the bitvector
            const uint64_t* data = m_v->data() + (sbv>>6);
            uint8_t offset = sbv&63;
            if (block != 0) {
                finalpos += block_value;
                data += (block_value>>6);
                offset += (block_value&63);
            }
            if (offset>=64) {
                data++;
                offset -= 64;
            }

            // substract the number of elements we skipped
            i -= superblock*m_elems_per_superblock;
            i -= block*elements_per_block;

            // first align to 64 bits so we can scan faster after
            uint64_t carry_unused;
            size_type cnt = select_support_clark_trait<b,pattern_len>::args_in_the_first_word(*data,offset+1,carry_unused);
            if (cnt >= i) {
                return finalpos +
                       select_support_clark_trait<b,pattern_len>::ith_arg_pos_in_the_first_word(*data,i,offset+1,carry_unused) - offset;
            }
            finalpos += 64-offset;
            // scan
            do {
                finalpos += 64;
                i -= cnt;
                data++;
                cnt = select_support_clark_trait<b,pattern_len>::args_in_the_word(*data,carry_unused);
            } while (cnt < i);
            finalpos += select_support_clark_trait<b,pattern_len>::ith_arg_pos_in_the_word(*data,i,carry_unused);
            finalpos -= 64;

            // we found the integer we are looking for...
            return finalpos;
        }
    }

    return i;
}

// alias for select(i)
template<uint8_t b, uint8_t pattern_len>
inline const typename select_support_clarks<b,pattern_len>::size_type select_support_clarks<b,pattern_len>::operator()(size_type i)const
{
    return select(i);
}

template<uint8_t b, uint8_t pattern_len>
void select_support_clarks<b,pattern_len>::initData()
{
    m_count = 0;
    m_size_miniblocks = 0;
    m_size_blocks = 0;
    m_num_long_elems = 0;
    if (m_v==NULL) {
        m_logn = 0;
        m_loglogn = 0;
        m_elems_per_superblock = 0;
        m_longthreshold = 0;
    } else {
        m_logn = bit_magic::l1BP(m_v->size())+1;
        m_loglogn = bit_magic::l1BP(bit_magic::l1BP(m_v->size()))+1;
        m_elems_per_superblock = m_logn*m_loglogn;
        m_longthreshold = (m_logn*m_logn) * (m_loglogn*m_loglogn);
    }
}

template<uint8_t b, uint8_t pattern_len>
void select_support_clarks<b,pattern_len>::set_vector(const int_vector<1>* v)
{
    m_v = v;
}


/*
 *  return the number of long blocks in the data structure
 */
template<uint8_t b, uint8_t pattern_len>
void select_support_clarks<b,pattern_len>::block_bytes(size_t& sb,size_t& lb,size_t& bb,size_t& mb) const
{
    sb = util::get_size_in_bytes(m_superblocks);
    lb = util::get_size_in_bytes(m_longblocks);
    bb = m_size_blocks/8;
    mb = m_size_miniblocks/8;
}


/*
 * store the support structure. use two small bitvectors marking how we stored
 * each superblock.
 */
template<uint8_t b, uint8_t pattern_len>
typename select_support_clarks<b,pattern_len>::size_type select_support_clarks<b,pattern_len>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;

    written_bytes += util::write_member(m_count,out,child,"count");
    written_bytes += util::write_member(m_elems_per_superblock,out,child,"elem_per_block");
    written_bytes += util::write_member(m_longthreshold,out,child,"longthreshold");
    written_bytes += util::write_member(m_numsuperblocks,out,child,"numsuperblocks");
    written_bytes += util::write_member(m_logn,out,child,"logn");
    written_bytes += util::write_member(m_loglogn,out,child,"loglogn");
    written_bytes += util::write_member(m_size_blocks,out,child,"size_blocks");
    written_bytes += util::write_member(m_size_miniblocks,out,child,"size_miniblocks");
    written_bytes += util::write_member(m_num_long_elems,out,child,"num_long_elems");

    written_bytes += m_superblocks.serialize(out, child, "superblock");
    written_bytes += m_longblocks.serialize(out, child, "longblocks");
    written_bytes += m_blocks_and_mini.serialize(out, child, "blocks_and_mini");
    written_bytes += m_block_pos.serialize(out, child, "block_pos");

    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<uint8_t b, uint8_t pattern_len>
void select_support_clarks<b,pattern_len>::load(std::istream& in, const int_vector<1>* v)
{
    set_vector(v);
    initData();

    util::read_member(m_count,in);
    util::read_member(m_elems_per_superblock,in);
    util::read_member(m_longthreshold,in);
    util::read_member(m_numsuperblocks,in);
    util::read_member(m_logn,in);
    util::read_member(m_loglogn,in);
    util::read_member(m_size_blocks,in);
    util::read_member(m_size_miniblocks,in);
    util::read_member(m_num_long_elems,in);

    m_superblocks.load(in);
    m_longblocks.load(in);
    m_blocks_and_mini.load(in);
    m_block_pos.load(in);
}

}

#endif
