/**
 * @file ht_pair_ZZ_long.hpp
 * @author Michael Jacobson
 * @version $Header$
 */

#ifndef ANTL_HT_PAIR_ZZ_LONG_H
#define ANTL_HT_PAIR_ZZ_LONG_H

#include <ANTL/HashTable/HashTable.hpp>
#include <ANTL/PairZZLong.hpp>

namespace ANTL {

// NTL_hentry_decl (pair_ZZ_long,hentry_pair_ZZ_long)
typedef HEntry<PairZZLong> hentry_pair_ZZ_long;
// NTL_hash_table_decl
// (pair_ZZ_long,hentry_pair_ZZ_long,hash_table_pair_ZZ_long)
typedef HashTable<PairZZLong> hash_table_pair_ZZ_long;

} // namespace ANTL

#endif // guard
