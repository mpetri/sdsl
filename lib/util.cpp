/* sdsl - succinct data structures library
    Copyright (C) 2009 Simon Gog

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

#include "sdsl/util.hpp"
#include "sdsl/structure_tree.hpp"
#include "cxxabi.h"
#include <vector>

namespace sdsl
{

namespace util
{

off_t get_file_size(const char* file_name)
{
    struct stat filestatus;
    stat(file_name, &filestatus);
    return filestatus.st_size;
}


void file::write_text(const char* file_name, const char* c, uint64_t len)
{
    std::ofstream out(file_name);
    if (out) {
        out.write(c, len);
        out.close();
    }
}

uint64_t file::read_text(const char* file_name, char*& c, bool trunc, uint64_t lim)
{
    if (c != NULL) {
        delete [] c;
        c = NULL;
    }
    uint64_t n = get_file_size(file_name) + 1; // add one for the 0 byte
    if (trunc and lim+1 < n) {
        n = lim+1;
    }
    std::ifstream in;
    in.open(file_name);
    if (in) {
        c = new char[n];
        c[n-1] = '\0';
        char* cp = c;
        in.read(cp, n-1);
        if (n > 1 and c[n-2] == 0)
            return n-1; // last byte was already a null byte
        else
            return n; // added 0 byte
    }
    return 0;
}

std::string encode_base64(const std::string& txt)
{
    return encode_base64(txt.c_str(),txt.length());
}


std::string decode_base64(const std::string& txt)
{
    return decode_base64(txt.c_str(),txt.length());
}

std::string encode_base64(const char* txt,size_t n)
{
    // calc result size by adding padding
    size_t res_size = ((n/3)*4)+4;
    std::string base64_str;
    base64_str.resize(res_size);
    assert(res_size % 4 = 0);

    // process input. 3 ascii symbols become 4 base64 symbols
    size_t i,j;
    uint32_t first,second,third,octets;
    for (i=0,j=0; i<(n-2); i+=3) {
        first = txt[i]; second = txt[i+1]; third = txt[i+2];
        octets = (first << 16) + (second << 8) + third;
        base64_str[j++] = encode_table_base64[(octets >> 18) & 0x3F ];
        base64_str[j++] = encode_table_base64[(octets >> 12) & 0x3F ];
        base64_str[j++] = encode_table_base64[(octets >> 6) & 0x3F ];
        base64_str[j++] = encode_table_base64[ octets & 0x3F ];
    }

    // padding?
    first = second = third = 0;
    switch (n%3) {
        case 1:
            first = txt[n-1];
            octets = (first << 16) + (second << 8) + third;
            base64_str[j++] = encode_table_base64[(octets >> 18) & 0x3F ];
            base64_str[j++] = encode_table_base64[(octets >> 12) & 0x3F ];
            base64_str[j++] = base64_padding_sym;
            base64_str[j] = base64_padding_sym;
            break;
        case 2:
            first = txt[n-2];
            second = txt[n-1];
            octets = (first << 16) + (second << 8) + third;
            base64_str[j++] = encode_table_base64[(octets >> 18) & 0x3F ];
            base64_str[j++] = encode_table_base64[(octets >> 12) & 0x3F ];
            base64_str[j++] = encode_table_base64[(octets >> 6) & 0x3F ];
            base64_str[j] = base64_padding_sym;
            break;
        case 0:
            // no padding needed
            j--; // we added one too many symbols before
            break;
    }
    base64_str.resize(j+1);
    return base64_str;
}


std::string decode_base64(const char* txt,size_t n)
{
    assert(n%4 == 0);
    size_t res_size = n/4 * 3;
    std::string ascii_str;
    ascii_str.resize(res_size);

    uint32_t first,second,third,forth;
    size_t i,j;
    // decode 4 base64 symbols into 3 ascii symbols
    for (i=0,j=0; i<n; i+=4) {
        first = decode_table_base64[(uint8_t)txt[i]];
        second = decode_table_base64[(uint8_t)txt[i+1]];
        third = decode_table_base64[(uint8_t)txt[i+2]];
        forth = decode_table_base64[(uint8_t)txt[i+3]];
        uint32_t recover = (first << 18) + (second << 12) + (third << 6) + forth;

        ascii_str[j++] = (recover >> 16) & 0xFF;
        if (txt[i+2] != '=') ascii_str[j++] = (recover >> 8) & 0xFF;
        if (txt[i+3] != '=') ascii_str[j++] =  recover & 0xFF;
    }
    ascii_str.resize(j); // fix the size
    return ascii_str;
}


uint64_t _id_helper::id = 0;

std::string basename(const std::string& file_name)
{
    char* c = strdup((const char*)file_name.c_str());
    std::string res = std::string(::basename(c));
    free(c);
    return res;
}

std::string dirname(const std::string& file_name)
{
    char* c = strdup((const char*)file_name.c_str());
    std::string res = std::string(::dirname(c));
    free(c);
    return res;
}

uint64_t get_pid()
{
    return getpid();
}

std::string demangle(const char* name)
{
#ifndef HAVE_CXA_DEMANGLE
    char buf[4096];
    size_t size = 4096;
    int status = 0;
    abi::__cxa_demangle(name, buf, &size, &status);
    if (status==0)
        return std::string(buf);
    return std::string(name);
#else
    return std::string(name);
#endif
}

std::string demangle2(const char* name)
{
    std::string result = demangle(name);
    std::vector<std::string> words_to_delete;
    words_to_delete.push_back("sdsl::");
    words_to_delete.push_back("(unsigned char)");
    words_to_delete.push_back(", unsigned long");

    for (size_t k=0; k<words_to_delete.size(); ++k) {
        std::string w = words_to_delete[k];
        for (size_t i = result.find(w); i != std::string::npos; i = result.find(w, i)) {
            result.erase(i, w.length());
            ++i;
        }
    }
    size_t index = 0;
    std::string to_replace = "int_vector<1>";
    while ((index = result.find(to_replace, index)) != string::npos) {
        result.replace(index, to_replace.size(), "bit_vector");
    }
    return result;
}

void delete_all_files(tMSS& file_map)
{
    for (tMSS::iterator file_it=file_map.begin(); file_it!=file_map.end(); ++file_it) {
        std::remove(file_it->second.c_str());
    }
    file_map.clear();
}

std::string to_latex_string(unsigned char c)
{
    if (c == '_')
        return "\\_";
    else if (c == '\0')
        return "\\$";
    else
        return to_string(c);
}

template<>
size_t write_member<std::string>(const std::string& t, std::ostream& out, structure_tree_node* v, std::string name)
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(t));
    size_t written_bytes = 0;
    written_bytes += write_member(t.size(), out, child, "length");
    out.write(t.c_str(), t.size());
    written_bytes += t.size();
    structure_tree::add_size(v, written_bytes);
    return written_bytes;
}

template<>
void read_member<std::string>(std::string& t, std::istream& in)
{
    std::string::size_type size;
    read_member(size, in);
    char* buf = new char[size];
    in.read(buf, size);
    std::string temp(buf, size);
    delete [] buf;
    t.swap(temp);
}

template<>
bool load_from_file(void*& v, const char* file_name)
{
    return true;
}

bool load_from_file(char*& v, const char* file_name)
{
    if (v != NULL) {
        delete [] v;
        v = NULL;
    }
    std::ifstream in;
    in.open(file_name, std::ios::binary | std::ios::in);
    if (in) {
        const uint64_t SDSL_BLOCK_SIZE = (1<<20);
        uint64_t n=0, read = 0;
        char buf[SDSL_BLOCK_SIZE], *cp;
        do {
            in.read(buf, SDSL_BLOCK_SIZE);
            read = in.gcount();
            n+=read;
        } while (SDSL_BLOCK_SIZE == read);
        if (n==0)
            return false;
        v = new char[n+1];
        in.close();
        in.open(file_name);
        if (!in) {
            delete [] v;
            v = NULL;
            return false;
        }
        cp=v;
        do {
            in.read(cp, SDSL_BLOCK_SIZE);
            read = in.gcount();
            cp+= read;
        } while (SDSL_BLOCK_SIZE == read);
        *(v+n) = '\0';
        return true;
    } else
        return false;
}

void set_verbose()
{
    verbose = true;
}


}// end namespace util
}// end namespace sdsl

