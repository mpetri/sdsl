#include "sdsl/int_vector.hpp"
#include "sdsl/bitmagic.hpp"
#include "sdsl/util.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <cstdlib> // for rand()
#include <string>

using namespace sdsl;

TEST(UtilTest, base64_encode)
{
    ASSERT_STREQ("YW55IGNhcm5hbCBwbGVhc3VyZS4=" , util::encode_base64("any carnal pleasure.").c_str());
    ASSERT_STREQ("YW55IGNhcm5hbCBwbGVhc3VyZQ==" , util::encode_base64("any carnal pleasure").c_str());
    ASSERT_STREQ("YW55IGNhcm5hbCBwbGVhc3Vy" , util::encode_base64("any carnal pleasur").c_str());
    ASSERT_STREQ("YW55IGNhcm5hbCBwbGVhc3U=" , util::encode_base64("any carnal pleasu").c_str());
    ASSERT_STREQ("YW55IGNhcm5hbCBwbGVhcw==" , util::encode_base64("any carnal pleas").c_str());


}

TEST(UtilTest, base64_decode)
{
    ASSERT_STREQ(util::decode_base64("YW55IGNhcm5hbCBwbGVhc3VyZS4=").c_str() , "any carnal pleasure.");
    ASSERT_STREQ(util::decode_base64("YW55IGNhcm5hbCBwbGVhc3VyZQ==").c_str() , "any carnal pleasure");
    ASSERT_STREQ(util::decode_base64("YW55IGNhcm5hbCBwbGVhc3Vy").c_str() , "any carnal pleasur");
    ASSERT_STREQ(util::decode_base64("YW55IGNhcm5hbCBwbGVhc3U=").c_str() , "any carnal pleasu");
    ASSERT_STREQ(util::decode_base64("YW55IGNhcm5hbCBwbGVhcw==").c_str() , "any carnal pleas");
}

TEST(UtilTest, base64_encodedecode)
{
    for (size_t i=0; i<2000; i++) {
        std::string str;
        str.resize(10+ (rand()%256));
        for (size_t j=0; j<str.size(); j++) {
            int sym = rand() % 256;
            while (isprint(sym) == 0) sym = rand() % 256;
            str[j] = (unsigned char) sym;
        }
        ASSERT_STREQ(str.c_str() ,  util::decode_base64(util::encode_base64(str)).c_str());
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
