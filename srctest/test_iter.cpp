#include "catch.hpp"
#include <iter.hpp>

TEST_CASE("iterator over unmappable regions") {
  //  1  2  3  4  5
  //  3  4  5  6  8
  
  //  1  0  0  1  0
  //  0  0  1  1  1
  
  double data[10] = { 1, 3, 2, 4, 3, 5, 4, 6, 5, 8 };
  int missing[10] = { 1, 0, 0, 0, 0, 1, 1, 1, 0, 1 };
  int slot_dim[2] = { 1, 1 };
  int length = 5;
  int n_slots = 2;
  
  Iter * iter = new Iter(length, n_slots, slot_dim, data, 0, NULL, NULL, missing);

  SECTION("check first slot subiterators") {
    std::vector<Iter> * sub_iters = iter->sub_iterators(0);
    
    // appropriate length
    REQUIRE( sub_iters->size() == 2 );
    
    // values
    Iter iter_slot0_0 = sub_iters->at(0);
    
    REQUIRE( iter_slot0_0.emission(0) == 2 );
    CHECK( iter_slot0_0.next() );
    REQUIRE( iter_slot0_0.emission(0) == 3 );
    CHECK_FALSE( iter_slot0_0.next() );

    Iter iter_slot0_1 = sub_iters->at(1);
    
    REQUIRE( iter_slot0_1.emission(0) == 5 );
    CHECK_FALSE( iter_slot0_1.next() );
    
    delete sub_iters;
  }
  
  SECTION("check second slot subiterators") {
    std::vector<Iter> * sub_iters = iter->sub_iterators(1);
    
    // appropriate length
    REQUIRE( sub_iters->size() == 1 );

    // values
    Iter iter_slot1_0 = sub_iters->at(0);
    
    REQUIRE( iter_slot1_0.emission(1) == 3);
    CHECK( iter_slot1_0.next());
    REQUIRE( iter_slot1_0.emission(1) == 4);
    CHECK_FALSE( iter_slot1_0.next());
  }
  
  delete iter;
}
