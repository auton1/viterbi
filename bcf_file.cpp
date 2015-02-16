/*
 * bcf_file.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: amarcketta
 */

#include "bcf_file.h"

char bcf_21[5] = {'B','C','F',(int8_t)2,(int8_t)1};
char bcf_22[5] = {'B','C','F',(int8_t)2,(int8_t)2};

bcf_file::bcf_file(const parameters &p, string &vcf_filename, int idx)
{
	filename = vcf_filename;

	big_endian = is_big_endian();
	is_BGZF = false; stream = p.stream_in;
	N_entries = 0; N_kept_entries = 0;
	meta_data = header();

	if (stream)
	{
		is_BGZF = true;
		open_gz();
	}
	else
		open();

	check_bcf();
	read_header();

	include_indv = vector<bool>(meta_data.N_indv,true);
}

bcf_file::~bcf_file()
{
	close();
}

void bcf_file::open()
{
	int ret;

	if (filename.substr(filename.size()-4) == ".vcf")
			LOG.error("Filename ends in '.vcf'. Shouldn't you be using --vcf?\n");

	if (filename.substr(filename.size()-7) == ".vcf.gz")
			LOG.error("Filename ends in '.vcf.gz'. Shouldn't you be using --gzvcf?\n");

	ret = bgzf_is_bgzf(filename.c_str());

	if (ret == 1)
		is_BGZF = true;
	else
		is_BGZF = false;

	if (is_BGZF)
		open_gz();
	else
	{
		file_tmp.open(filename.c_str(), ios::in);
		if (!file_tmp.is_open())
			LOG.error("Could not open VCF file: " + filename, 0);
		file_in = &file_tmp;
	}
}

void bcf_file::open_gz()
{
	int ret;
	gzMAX_LINE_LEN = 1024*1024;

	if (stream)
		gzfile_in = gzdopen(fileno(stdin), "r");
	else
		gzfile_in = gzopen(filename.c_str(), "rb");

	if (gzfile_in == NULL)
		LOG.error("Could not open BGZF BCF file: " + filename, 0);
	#ifdef ZLIB_VERNUM
		string tmp(ZLIB_VERSION);
		LOG.printLOG("Using zlib version: " + tmp + "\n");
		#if (ZLIB_VERNUM >= 0x1240)
			ret = gzbuffer(gzfile_in, gzMAX_LINE_LEN); // Included in zlib v1.2.4 and makes things MUCH faster
			if (ret != 0)
				LOG.warning("Unable to change zlib buffer size.");
		#else
			LOG.printLOG("Versions of zlib >= 1.2.4 will be *much* faster when reading compressed BCF files.\n");
		#endif
	#endif
}

void bcf_file::check_bcf()
{
	char magic[5];

	read(magic, 5, 1);
	if ((strcmp(magic, bcf_21) == 0) && (strcmp(magic, bcf_22) == 0))
		LOG.error("Does not appear to be a BCF file\n");
}

void bcf_file::close()
{
	if (!stream && is_BGZF)
		gzclose(gzfile_in);
}

void bcf_file::get_entry(vector<char> &out)
{
	uint32_t size_int[2];
	int ret, read_size = 0;

	ret = read(&size_int[0], 2, sizeof(uint32_t) );
	read_size = size_int[0] + size_int[1];

	if (ret)
	{
		out.resize(read_size+2*sizeof(uint32_t));
		memcpy(&out[0], size_int, 2*sizeof(uint32_t));
		read(&out[2*sizeof(uint32_t)], 1, read_size);
	}
	else
		out.resize(0);
}

entry* bcf_file::get_entry_object()
{
	return new bcf_entry(meta_data, include_indv);
}

int bcf_file::read(void *buffer, unsigned int len, size_t size)
{
	int ret;

	if (is_BGZF)
	{
		ret = gzread(gzfile_in, buffer, size*len);
		ret = (ret == (int)(len*size) );
	}
	else
	{
		file_in->read((char*)buffer, size*len);
		ret = !file_in->eof();
	}

	if ((big_endian) && (size > 1)) // Note: don't both swapping character arrays - BCF is defined as little endian.
	{
		unsigned int ui;
		for (ui=0; ui<len; ui++)
			ByteSwap((unsigned char *)buffer+(size*ui), size);
	}
	return ret;
}

void bcf_file::read_header()
{
	uint32_t len_text;

	read(&len_text, 1, sizeof(uint32_t));
	char *header_array = new char[(uint32_t)len_text];
	read(header_array, len_text, 1);

	string header(header_array);
	delete [] header_array;
	string line;
	unsigned int line_index = 0;
	header.erase( header.find_last_not_of(" \f\n\r\t\v\0" ) + 1 );

	istringstream iss(header);
	line_index += meta_data.add_FILTER_descriptor("ID=PASS,Description=PASS", line_index);

	while (getline(iss, line))
	{
		if (line[0] == '#')
		{
			if (line[1] == '#')
				meta_data.parse_meta(line, line_index);
			else
				meta_data.parse_header(line);
		}
	}
}

bool bcf_file::eof()
{
	if (is_BGZF)
		return gzeof(gzfile_in);
	else
		return(file_in->eof());
}
