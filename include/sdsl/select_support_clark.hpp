/*! \file select_support_clark.hpp
    \brief select_support_clark.hpp contains classes that support a sdsl::bit_vector with constant time select information.
	\author Matthias Petri
*/
#ifndef INCLUDED_SDSL_SELECT_SUPPORT_CLARK
#define INCLUDED_SDSL_SELECT_SUPPORT_CLARK

#include "int_vector.hpp"
#include "select_support.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{


template<uint8_t bit_pattern, uint8_t pattern_len>
struct select_support_clark_trait {
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
struct select_support_clark_trait<0,1> {
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
struct select_support_clark_trait<1,1> {
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
struct select_support_clark_trait<10,2> {
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
struct select_support_clark_trait<01,2> {
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

template<uint8_t b=1, uint8_t pattern_len=1>
class select_support_clark : public select_support
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
        int_vector<0>* m_longblocks;        // log(n) for each position
        int_vector<0>* m_blocks;            // log(r) for each relative position
        int_vector<0>* m_miniblocks;        // log(r') for each relative position
    private:
        void copy(const select_support_clark<b, pattern_len>& ss);
        void construct();
        void initData();
        inline size_t superblock_width(size_type superblock) const;
        inline size_t block_width(size_type superblock,size_type block) const;
        inline size_t elements_per_block(size_type superblock,size_type block) const; 

        void process_block(size_type superblock,size_type block,int_vector<64>& pos,int_vector<0>*& miniblocks);
        void store_longblock(int_vector<64>& pos,size_type superblock);
        void store_blocks(size_type superblock,int_vector<64>& pos,size_type width);
        void store_miniblock(size_type superblock,size_type block,int_vector<64>& pos,int_vector<0>*& miniblocks,
                            size_type skip,size_type num_elems);
        void process_superblock(int_vector<64>& pos,size_type superblock,bool last=false);
        void compress_mini_blocks(size_type superblock,int_vector<0>*& miniblocks);
    public:
        explicit select_support_clark(const int_vector<1>* v=NULL);
        select_support_clark(const select_support_clark<b,pattern_len>& ss);
        ~select_support_clark();
        void init(const int_vector<1>* v=NULL);

        inline const size_type select(size_type i) const;
        //! Alias for select(i).
        inline const size_type operator()(size_type i)const;
        
        size_type serialize(std::ostream& out, structure_tree_node* v=NULL, std::string name="")const;
        void load(std::istream& in, const int_vector<1>* v=NULL);
        void set_vector(const int_vector<1>* v=NULL);
        select_support_clark<b, pattern_len>& operator=(const select_support_clark& ss);

        //! Swap operator
        /*! This swap operator swaps two select_support_clarks in constant time.
         */
        void swap(select_support_clark<b, pattern_len>& ss);

};


template<uint8_t b, uint8_t pattern_len>
void select_support_clark<b,pattern_len>::construct()
{
    m_longblocks = NULL;
    m_blocks = NULL;
    m_miniblocks = NULL;
    m_count	= 0;
    m_numsuperblocks = 0;
}

template<uint8_t b, uint8_t pattern_len>
select_support_clark<b,pattern_len>::select_support_clark(const int_vector<1>* f_v):select_support(f_v)
{
    construct();
    init(f_v);
}

template<uint8_t b, uint8_t pattern_len>
select_support_clark<b,pattern_len>::select_support_clark(const select_support_clark& ss):select_support(ss.m_v)
{
    init();
    copy(ss);
}

template<uint8_t b, uint8_t pattern_len>
select_support_clark<b, pattern_len>& select_support_clark<b,pattern_len>::operator=(const select_support_clark& ss)
{
    if (this != &ss) {
        copy(ss);
    }
    return *this;
}

template<uint8_t b, uint8_t pattern_len>
void select_support_clark<b,pattern_len>::swap(select_support_clark& ss)
{
    std::swap(m_logn, ss.m_logn);
    std::swap(m_loglogn, ss.m_loglogn);
    std::swap(m_longthreshold, ss.m_longthreshold);
    std::swap(m_elems_per_superblock, ss.m_elems_per_superblock);
    std::swap(m_longblocks, ss.m_longblocks);
    std::swap(m_blocks, ss.m_blocks);
    std::swap(m_miniblocks, ss.m_miniblocks);
    std::swap(m_count, ss.m_count);
    std::swap(m_numsuperblocks, ss.m_numsuperblocks);
    m_superblocks.swap(ss.m_superblocks);
}

template<uint8_t b, uint8_t pattern_len>
void select_support_clark<b,pattern_len>::copy(const select_support_clark<b, pattern_len>& ss)
{
    m_logn			 = ss.m_logn;
    m_loglogn        = ss.m_loglogn;
    m_longthreshold  = ss.m_longthreshold;
    m_elems_per_superblock = ss.m_elems_per_superblock;
    m_superblocks    = ss.m_superblocks;
    m_count          = ss.m_count;
    m_numsuperblocks = ss.m_numsuperblocks;

    // copy the pointers
    if( m_longblocks ) delete [] m_longblocks;
    if( m_blocks ) delete [] m_blocks;
    if( m_miniblocks ) delete [] m_miniblocks;

    if( ss.m_longblocks != NULL ) {
        m_longblocks = new int_vector<0>[m_numsuperblocks];
        for (size_type i=0; i<m_numsuperblocks; ++i) {
            m_longblocks[i] = ss.m_longblocks[i];
        }
    }

    if( ss.m_blocks != NULL ) {
        m_blocks = new int_vector<0>[m_numsuperblocks];
        for (size_type i=0; i<m_numsuperblocks; ++i) {
            m_blocks[i] = ss.m_blocks[i];
        }
    }

    if( ss.m_miniblocks!=NULL ) {
        m_miniblocks = new int_vector<0>[m_numsuperblocks];
        for (size_type i=0; i<m_numsuperblocks; ++i) {
            m_miniblocks[i] = ss.m_miniblocks[i];
        }
    }
}

template<uint8_t b, uint8_t pattern_len>
select_support_clark<b,pattern_len>::~select_support_clark()
{
    if( m_longblocks != NULL ) delete [] m_longblocks;
    if( m_blocks != NULL ) delete [] m_blocks;
    if( m_miniblocks != NULL ) delete [] m_miniblocks;
}

/**
 *  returns the size of a specific superblock
 */
template<uint8_t b, uint8_t pattern_len>
size_t select_support_clark<b,pattern_len>::superblock_width(size_type superblock)  const
{
    if(superblock == m_numsuperblocks-1 ) 
        return v->size() - m_superblocks[superblock];
    else 
        return m_superblocks[superblock+1] - m_superblocks[superblock];
}

/**
 *  returns the size of a specific block in a superblock.
 */
template<uint8_t b, uint8_t pattern_len>
size_t select_support_clark<b,pattern_len>::block_width(size_type superblock,size_type block) const 
{
    if(block==0)
        return m_blocks[superblock][0];

    if(block==m_blocks[superblock].size()) {
        return m_superblocks[superblock+1] - 
               ( m_superblocks[superblock] + m_blocks[superblock][block-1]);
    }
    return m_blocks[superblock][block] - m_blocks[superblock][block-1];
}

/**
 *  returns the number of elements in a specific block.
 */
template<uint8_t b, uint8_t pattern_len>
size_t select_support_clark<b,pattern_len>::elements_per_block(size_type superblock,size_type block) const
{
    size_type r = m_superblocks[superblock+1] - m_superblocks[superblock];
    size_type logr = bit_magic::l1BP(r)+1;
    if(block==m_blocks[superblock].size()) {
        return m_elems_per_superblock - m_blocks[superblock].size()*(logr*m_loglogn);
    }
    return (logr*m_loglogn);
}

/**
 *  store all positions in a superblock explicitly using log(n) bits each.
 *  we create a "long" block.
 */
template<uint8_t b, uint8_t pattern_len>
void select_support_clark<b,pattern_len>::store_longblock(int_vector<64>& pos,
                                            size_type superblock) 
{
    // we store all positions explicitly using log(n) bits
    if (m_longblocks == NULL) 
        m_longblocks = new int_vector<0>[m_numsuperblocks];
    m_longblocks[superblock] = int_vector<0>( pos.size()-1 , 0 , m_logn );

    // copy the positions into the long block
    for (size_type j=0; j < pos.size()-1; j++) {
        m_longblocks[superblock][j] = pos[j+1];
    }
}

/**
 *   we subdivide a superblock and store position samples using log(r) bits each.
 */
template<uint8_t b, uint8_t pattern_len>
void select_support_clark<b,pattern_len>::store_blocks(size_type superblock,
    int_vector<64>& pos,size_type width) 
{
    if (m_blocks == NULL) 
        m_blocks = new int_vector<0>[m_numsuperblocks];

    // we don't count the first position as we store it explicitly
    size_type subblock_size = width * m_loglogn;
    size_type num_blocks = pos.size()/subblock_size;
    if(pos.size() % subblock_size == 0) num_blocks--;

    if(num_blocks) {
        m_blocks[superblock] = int_vector<0>( num_blocks , 0 , width );

        // we sample every blocksize-th position and store the relative
        // position in logr bits
        size_type start = pos[0];
        for (size_type i = 0; i < num_blocks; ++i) {
            m_blocks[superblock][i] = pos[(i+1)*subblock_size]-start;
        }
    }

}

/**
 *   decide if we have to store miniblocks for a given block or not
 */
template<uint8_t b, uint8_t pattern_len>
void 
select_support_clark<b,pattern_len>::process_block(size_type superblock,
    size_type block,int_vector<64>& pos,int_vector<0>*& miniblocks) 
{
    // check how big the block is to decide if we store something
    size_type sb_width = superblock_width(superblock); // r
    size_type blocksize = block_width(superblock,block); // r'
    size_type bits_block = bit_magic::l1BP(blocksize)+1; // log(r')
    size_type bits_super = bit_magic::l1BP(sb_width)+1;  // log(r)
    size_type threshold = bits_block*bits_super*m_loglogn*m_loglogn;

    // we only store mini blocks if we are larger than the threshold
    // of log(r')*log(r)*loglog(n)
    if( blocksize >= threshold ) {
        size_type num_elems = bits_super*m_loglogn;
        if(block == m_blocks[superblock].size()) {
            // last block might have less elements
            num_elems = pos.size()%num_elems;
        }
        size_type skip_positions = block*bits_super*m_loglogn;
        store_miniblock(superblock,block,pos,miniblocks,skip_positions,num_elems);
    }
}

/**
 *   store all positions in a given block explicitly relative to the subblock
 *   using log(r') bits each.
 */
template<uint8_t b, uint8_t pattern_len>
void 
select_support_clark<b,pattern_len>::store_miniblock(size_type superblock,
    size_type block,int_vector<64>& pos,int_vector<0>*& miniblocks,
    size_type skip,size_type num_elems) 
{
    size_type blocksize = block_width(superblock,block); // r'
    size_type bits_block = bit_magic::l1BP(blocksize)+1; // log(r')

    // store each relative position in log(r') bits
    miniblocks[block] = int_vector<0>( num_elems , 0 , bits_block );

    size_type start = m_superblocks[superblock];
    if(block!=0)
        start += m_blocks[superblock][block-1];

    for (size_type i = 1; i < num_elems; ++i) {
        miniblocks[block][i-1] = pos[skip+i] - start;
    }
}

/**
 *   process a superblock. first decide if we have a long block or if we
 *   divide into subblocks. 
 */
template<uint8_t b, uint8_t pattern_len>
void select_support_clark<b,pattern_len>::process_superblock(int_vector<64>&
    pos,size_type superblock,bool last) 
{
    // the "width" of the block
    size_type r = superblock_width(superblock);

    if(  last || r >= m_longthreshold  ) {
        store_longblock(pos,superblock);
    } else {
        // subdivide the super block
        size_type logr = bit_magic::l1BP(r)+1;
        store_blocks(superblock,pos,logr);

        // for each block we might have to subdivide again
        int_vector<0>* mini_blocks = new int_vector<0>[m_blocks[superblock].size()+1];
        for (size_type i = 0; i <= m_blocks[superblock].size(); ++i) {
            process_block(superblock,i,pos,mini_blocks);
        }

        // if we created miniblocks compress them to save space
        compress_mini_blocks(superblock,mini_blocks);
        delete [] mini_blocks;
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
select_support_clark<b,pattern_len>::compress_mini_blocks(size_type superblock,
    int_vector<0>*& miniblocks)
{

    uint64_t mask = 0;
    uint8_t stored = 0;
    // check if we have to do something first and count the size of the
    // resulting vector
    size_type bits = 0;
    for (size_type i = 0; i < m_blocks[superblock].size()+1; ++i) {
        if( miniblocks[i].size() > 0 ) {
            bits += miniblocks[i].size()*miniblocks[i].get_int_width();
            mask |= (1ULL << i);
            stored++;
        }
    }

    if(stored) {
        size_type header_size = sizeof(mask) + sizeof(stored) + stored*32;
        size_type bits_per_pos = bit_magic::l1BP(header_size+bits)+1;

        // calc the total bitvector size which includes the header
        bits += sizeof(mask)*8;      // mask
        bits += sizeof(stored)*8;    // how many blocks are stored
        bits += sizeof(stored)*8;    // the size of each starting pos
        bits += stored*bits_per_pos; // starting pos of each mini block


        if(m_miniblocks== NULL)
            m_miniblocks = new int_vector<0>[m_numsuperblocks];
        m_miniblocks[superblock] = int_vector<0>(bits,0,1);

        // write 
        uint64_t* data = (uint64_t*) m_miniblocks[superblock].data();
        uint8_t offset = 0;
        size_type written = 0;

        // write header 
        bit_magic::write_int_and_move(data,mask,offset,sizeof(mask)*8);
        written += sizeof(mask)*8;
        bit_magic::write_int_and_move(data,stored,offset,sizeof(stored)*8 );
        written += sizeof(stored)*8;
        bit_magic::write_int_and_move(data,bits_per_pos,offset,sizeof(stored)*8);
        written += sizeof(stored)*8;

        size_type cur_offset = sizeof(mask)*8 + sizeof(stored)*8 + sizeof(stored)*8
                    + stored*bits_per_pos;


        // write length
        for (size_type i = 0; i < m_blocks[superblock].size()+1; ++i) {
            if( miniblocks[i].size() > 0 ) {
                bit_magic::write_int_and_move(data,cur_offset,offset,bits_per_pos);
                cur_offset += miniblocks[i].size()*miniblocks[i].get_int_width();
                written += bits_per_pos;
            }
        }

        // write data
        for (size_type i = 0; i < m_blocks[superblock].size()+1; ++i) {
            if( miniblocks[i].size() > 0 ) {
                for (size_type j = 0; j < miniblocks[i].size(); ++j) {
                    bit_magic::write_int_and_move(data,miniblocks[i][j],offset,
                                        miniblocks[i].get_int_width());
                    written += miniblocks[i].get_int_width();
                }
            }
        }
        assert( written == bits );
    }
}


template<uint8_t b, uint8_t pattern_len>
void select_support_clark<b,pattern_len>::init(const int_vector<1>* v)
{
    set_vector(v);
    initData();
    if (m_v==NULL)
        return;
    
    // Count the number of elements in the bit vector
    m_count = select_support_clark_trait<b,pattern_len>::arg_cnt(*v);
    if (m_count==0) // if there are no elements in the vector we are done...
        return;

    m_numsuperblocks = (m_count+(m_elems_per_superblock-1))/m_elems_per_superblock;
    m_superblocks = int_vector<0>(m_numsuperblocks, 0, m_logn);

    size_type cnt = 0;
    size_type cur_block_cnt = 0;
    size_type cur_superblock = 0;

    /* find the first elem for the first superblock */
    int_vector<64> element_positions(m_elems_per_superblock);
    for (size_type i=0; i < v->size(); ++i) {
        if (select_support_clark_trait<b,pattern_len>::found_arg(i, *v)) {
            m_superblocks[cur_superblock] = i;
            element_positions[0] = i;
            cur_block_cnt = 1;
            cnt = 1;
            break;
        }
    }
    cur_superblock++;

    // after we found the first element of the next superblock we can process
    // the previous superblock.
    for (size_type i=m_superblocks[0]+1; i < v->size(); ++i) {
        if (select_support_clark_trait<b,pattern_len>::found_arg(i, *v)) {
            if( cur_block_cnt == m_elems_per_superblock ) {
                m_superblocks[cur_superblock] = i;
                process_superblock(element_positions,cur_superblock-1);
                cur_superblock++;
                cur_block_cnt = 0;
            }
            element_positions[cur_block_cnt] = i;
            cnt++;
            cur_block_cnt++;
            if(cnt == m_count) break;
        }
    }
    // process the last superblock if there are elements left over
    if(cur_block_cnt!=0) {
        element_positions.resize(cur_block_cnt);
        process_superblock(element_positions,cur_superblock-1,true);
    }

    assert( m_numsuperblocks == cur_superblock );
}

template<uint8_t b, uint8_t pattern_len>
inline const typename select_support_clark<b,pattern_len>::size_type select_support_clark<b,pattern_len>::select(size_type i) const
{
    assert(i > 0);
    assert(i <= m_count);

    // numbering starts at 0 in the data structure
    i--;

    /* determine superblock */
    size_type superblock = i/m_elems_per_superblock;
    size_type block_offset = i%m_elems_per_superblock;
    size_type r = superblock_width(superblock);

    if(block_offset==0) { // we store the position explicitly as it is the
                          // first element of a superblock
        return m_superblocks[superblock];
    }

    if( superblock == m_numsuperblocks-1 || r >= m_longthreshold ) {
        // get the answer directly from the pre-stored positions in the long block
        return m_longblocks[superblock][block_offset-1];
    } else {
        // we are in a subdivided superblock. find the correct block
        size_type logr = bit_magic::l1BP(r)+1;
        size_type elements_per_block = logr*m_loglogn;
        size_type block = block_offset/elements_per_block;
        size_type block_size = block_width(superblock,block);
        size_type bits_block = bit_magic::l1BP(block_size)+1; // log(r')
        size_type threshold = bits_block*logr*m_loglogn*m_loglogn;

        // check if the answer is one of the ones stored in the blocks
        // anyway
        if( block_offset % (elements_per_block) == 0) {
            // if this is true, block always >= 1
            assert(block >= 1);
            return m_superblocks[superblock] + m_blocks[superblock][block-1];
        }

        if(block_size >= threshold) {   
            // process miniblock to get the answer
            const uint64_t* data = m_miniblocks[superblock].data();
            // determine which of the stored mini blocks we want to access
            size_type id = bit_magic::b1Cnt(data[0]&bit_magic::Li1Mask[block+1]);

            // get the starting pos of that block
            data++;
            uint64_t offsetwidth = bit_magic::read_int(data,8,8);
            uint64_t start_pos = bit_magic::read_int(data,16+((id-1)*offsetwidth)
                                                     ,offsetwidth);

            // skip to the correct element
            size_type subblock_offset = block_offset % elements_per_block;
            start_pos += ((subblock_offset-1)*bits_block);

            // now retrieve the number we are looking for
            const uint64_t* mb_data = (m_miniblocks[superblock].data() +
                                        (start_pos >> 6));
            uint8_t offset = start_pos&63;
            uint64_t elem_pos = bit_magic::read_int(mb_data,offset,bits_block);

            // finally add up all the numbers
            if(block==0) return m_superblocks[superblock] + elem_pos;
            else return m_superblocks[superblock] + m_blocks[superblock][block-1] + elem_pos;
        } else {
            size_type finalpos = m_superblocks[superblock];
            size_type word = (m_superblocks[superblock]>>6);
            // perform a scan through the bitvector
            

            // determine the starting position in the bitvector
            const uint64_t* data = m_v->data() + (m_superblocks[superblock]>>6);
            uint8_t offset = m_superblocks[superblock]&63;
            if(block != 0) { 
                finalpos += m_blocks[superblock][block-1];
                data += (m_blocks[superblock][block-1]>>6);
                offset += (m_blocks[superblock][block-1]&63);
                word += (m_blocks[superblock][block-1]>>6);
            }
            if(offset>=64) {
                data++;
                offset -= 64;
                word++;
            } 

            // substract the number of elements we skipped
            i -= superblock*m_elems_per_superblock;
            i -= block*elements_per_block;

            // first align to 64 bits so we can scan faster after
            uint64_t carry_unused;
            size_type cnt = select_support_clark_trait<b,pattern_len>::args_in_the_first_word(*data,offset+1,carry_unused);
            if(cnt >= i) {
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
            } while( cnt < i ); 
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
inline const typename select_support_clark<b,pattern_len>::size_type select_support_clark<b,pattern_len>::operator()(size_type i)const
{
    return select(i);
}

template<uint8_t b, uint8_t pattern_len>
void select_support_clark<b,pattern_len>::initData()
{
    m_count = 0;
    if (m_v==NULL) {
        m_logn = 0;
        m_loglogn = 0;
        m_elems_per_superblock = 0;
        m_longthreshold = 0;
    } else {
        m_logn = bit_magic::l1BP(m_v->size())+1; 
        m_loglogn = bit_magic::l1BP( bit_magic::l1BP(m_v->size()) )+1;
        m_elems_per_superblock = m_logn*m_loglogn;
        m_longthreshold = (m_logn*m_logn) * (m_loglogn*m_loglogn);
    }

    if(m_longblocks) delete [] m_longblocks;
    if(m_blocks) delete [] m_blocks;
    if(m_miniblocks) delete [] m_miniblocks;

    m_longblocks = m_blocks = m_miniblocks = NULL;
}

template<uint8_t b, uint8_t pattern_len>
void select_support_clark<b,pattern_len>::set_vector(const int_vector<1>* v)
{
    m_v = v;
}

/*
 * store the support structure. use two small bitvectors marking how we stored
 * each superblock.
 */
template<uint8_t b, uint8_t pattern_len>
typename select_support_clark<b,pattern_len>::size_type select_support_clark<b,pattern_len>::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;

    written_bytes += util::write_member(m_count,out,child,"count");
    written_bytes += util::write_member(m_elems_per_superblock,out,child,"elem_per_block");
    written_bytes += util::write_member(m_longthreshold,out,child,"longthreshold");
    written_bytes += util::write_member(m_numsuperblocks,out,child,"numsuperblocks");
    written_bytes += util::write_member(m_logn,out,child,"logn");
    written_bytes += util::write_member(m_loglogn,out,child,"loglogn");

    written_bytes += m_superblocks.serialize(out, child, "superblock");

    int_vector<1> block_type_long(m_numsuperblocks);
    int_vector<1> block_type_mini(m_numsuperblocks);

    if (m_longblocks!=NULL) {
        for (size_type i=0; i< m_numsuperblocks; ++i)
            block_type_long[i] = !m_longblocks[i].empty();
    }

    if (m_miniblocks!=NULL) {
        for (size_type i=0; i< m_numsuperblocks; ++i)
            block_type_mini[i] = !m_miniblocks[i].empty();
    }
    written_bytes += block_type_long.serialize(out, child, "block_type_long");
    written_bytes += block_type_mini.serialize(out, child, "block_type_mini");

    size_type written_long = 0;
    size_type written_block = 0;
    size_type written_mini = 0;
    for (size_type i=0; i< m_numsuperblocks; ++i) {
        if(block_type_long[i])
            written_long+= m_longblocks[i].serialize(out);
        else {
            written_block+= m_blocks[i].serialize(out);
            if(block_type_mini[i]) {
                written_mini+= m_miniblocks[i].serialize(out);
            }
        }
    }
    written_bytes += written_long;
    written_bytes += written_block;
    written_bytes += written_mini;

    structure_tree_node* child_long = structure_tree::add_child(child, "longblock", util::class_name(m_longblocks));
    structure_tree::add_size(child_long, written_long);
    structure_tree_node* child_block = structure_tree::add_child(child, "block", util::class_name(m_blocks));
    structure_tree::add_size(child_block, written_block);
    structure_tree_node* child_mini = structure_tree::add_child(child, "miniblocks", util::class_name(m_miniblocks));
    structure_tree::add_size(child_mini, written_mini);

    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<uint8_t b, uint8_t pattern_len>
void select_support_clark<b,pattern_len>::load(std::istream& in, const int_vector<1>* v)
{
    set_vector(v);
    initData();

    util::read_member(m_count,in);
    util::read_member(m_elems_per_superblock,in);
    util::read_member(m_longthreshold,in);
    util::read_member(m_numsuperblocks,in);
    util::read_member(m_logn,in);
    util::read_member(m_loglogn,in);

    m_superblocks.load(in); 

    int_vector<1> block_type_long;
    int_vector<1> block_type_mini;

    block_type_long.load(in);
    block_type_mini.load(in);

    m_longblocks = new int_vector<0>[m_numsuperblocks];
    m_blocks = new int_vector<0>[m_numsuperblocks];
    m_miniblocks = new int_vector<0>[m_numsuperblocks];

    for (size_type i=0; i< m_numsuperblocks; ++i) {
        if(block_type_long[i])
            m_longblocks[i].load(in);
        else {
            m_blocks[i].load(in);
            if(block_type_mini[i]) m_miniblocks[i].load(in);
        }
    }
}

}

#endif
