/*
 * variant_file.h
 *
 *  Created on: Dec 12, 2012
 *      Author: amarcketta
 */

#ifndef VARIANT_FILE_H_
#define VARIANT_FILE_H_

#include <algorithm>
#include <bitset>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <functional>
#include <limits>
#include <set>
#include <sstream>
#include <map>
#include <numeric>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <utility>
#include <vector>
#include <zlib.h>

#include "parameters.h"
#include "entry.h"
#include "vcf_entry.h"
#include "bcf_entry.h"
#include "header.h"

extern output_log LOG;

using namespace std;

class variant_file
{
public:
	string filename;
	bool compressed;
	istream *file_in;
	ifstream file_tmp;
	unsigned int gzMAX_LINE_LEN;
	gzFile gzfile_in;

	header meta_data;
	vector<bool> include_indv;
	unsigned int N_entries;
	unsigned int N_kept_entries;

	int N_kept_individuals() const;
	int N_kept_sites() const;
	int N_total_sites() const;

	virtual void open() = 0;
	virtual void open_gz() = 0;
	virtual void close() = 0;
	virtual bool eof() = 0;

	virtual void get_entry(vector<char> &out) = 0;
	virtual entry* get_entry_object() = 0;
	void ByteSwap(unsigned char *b, int n) const;
	static inline bool is_big_endian() { long one= 1; return !(*((char *)(&one))); };

	void apply_filters(const parameters &params);
	void filter_individuals(const set<string> &indv_to_keep, const set<string> &indv_to_exclude, const string &indv_to_keep_filename, const string &indv_to_exclude_filename, bool keep_then_exclude=true);
	void filter_individuals_by_keep_list(const set<string> &indv_to_keep, const string &indv_to_keep_filename);
	void filter_individuals_by_exclude_list(const set<string> &indv_to_exclude, const string &indv_to_exclude_filename);
	void filter_individuals_randomly(int max_N_indv);

	void read_file(const parameters &params, vector<vector<bool> > &haplotypes, vector<unsigned int> &positions, vector<string> &haplotype_names, vector<string> &alleles);
	int add_file(const parameters &params, vector<vector<bool> > &haplotypes, vector<unsigned int> &positions, vector<string> &haplotype_names, vector<string> &alleles);

	void read_temp_site(ifstream &tmp_file, string &CHROM, int &POS, vector< pair<int,int> > &GTs);
	void read_big_temp_site(ifstream &tmp_file, string &CHROM, int &POS, int &alleles, vector< pair<int,int> > &GTs);
	void return_indv_union(variant_file &file2, map<string, pair< int, int> > &combined_individuals, const string &indv_ID_map_file="");

	void get_contigs(const std::string &contigs_file, vector<string> &contig_vector);
	virtual ~variant_file() = 0;
};

#endif /* VARIANT_FILE_H_ */
