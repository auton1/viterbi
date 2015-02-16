/*
 * vcf_file.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: amarcketta
 */

#include "vcf_file.h"

vcf_file::vcf_file(const parameters &p, string &vcf_filename, int idx)
{
	filename = vcf_filename;
	compressed = p.vcf_compressed[idx];
	stream = p.stream_in;

	gzMAX_LINE_LEN = 0;
	N_entries = 0; N_kept_entries = 0;
	meta_data = header();

	if (stream && compressed)
		open_gz();
	else if (stream)
	{
		char first = cin.peek();
		if (first == 0x1f)
			LOG.error("File starts with gzip magic string. Shouldn't you be using --gzvcf?\n");

		file_in = &std::cin;
	}
	else
		open();

	read_header();
	include_indv = vector<bool>(meta_data.N_indv,true);
}

vcf_file::~vcf_file()
{
	close();
}

void vcf_file::read_header()
{
	string line;
	unsigned int line_index = 0;
	line_index += meta_data.add_FILTER_descriptor("ID=PASS,Description=PASS", line_index);

	while (!eof())
	{
		read_line(line);
		if (line[0] == '#')
			if (line[1] == '#')
				meta_data.parse_meta(line, line_index);
			else
			{
				meta_data.parse_header(line);
				return;
			}
		else
			return;
	}
}

void vcf_file::open()
{
	struct stat buf;

	int i = stat(filename.c_str(), &buf);
	if (i != 0)
	{
		perror("stat error");
		LOG.error("Can't determine file type of " + filename, 0);
	}
	if (!S_ISREG(buf.st_mode))
		LOG.error("Does not appear to be a regular file: " + filename, 0);

	if (filename.substr(filename.size()-4) == ".bcf")
		LOG.error("Filename ends in '.bcf'. Shouldn't you be using --bcf?\n");

	if (!compressed)
	{
		if (filename.substr(filename.size()-3) == ".gz")
			LOG.error("Filename ends in '.gz'. Shouldn't you be using --gzvcf or --gzdiff?\n");
		file_tmp.open(filename.c_str(), ios::in);
		if (!file_tmp.is_open())
			LOG.error("Could not open VCF file: " + filename, 0);

		file_in = &file_tmp;
	}
	else
		open_gz();
}

void vcf_file::open_gz()
{
	gzMAX_LINE_LEN = 1024*1024;
	gz_readbuffer = new char[gzMAX_LINE_LEN];

	if (stream)
		gzfile_in = gzdopen(fileno(stdin), "r");
	else
		gzfile_in = gzopen(filename.c_str(), "rb");

	if (gzfile_in == NULL)
		LOG.error("Could not open GZVCF file: " + filename, 0);
	#ifdef ZLIB_VERNUM
		string tmp(ZLIB_VERSION);
		LOG.printLOG("Using zlib version: " + tmp + "\n");
		#if (ZLIB_VERNUM >= 0x1240)
			gzbuffer(gzfile_in, gzMAX_LINE_LEN); // Included in zlib v1.2.4 and makes things MUCH faster
		#else
			LOG.printLOG("Versions of zlib >= 1.2.4 will be *much* faster when reading zipped VCF files.\n");
		#endif
	#endif
}

void vcf_file::close()
{
	if (compressed)
	{
		gzclose(gzfile_in);
		delete [] gz_readbuffer;
	}
}

bool vcf_file::eof()
{
	bool out;
	if (!compressed)
		out = file_in->eof();
	else
		out = gzeof(gzfile_in);	// Returns 1 when EOF has previously been detected reading the given input stream, otherwise zero.
	return out;
}

void vcf_file::get_entry(vector<char> &out)
{
	out.resize(0);
	read_line(out);
}

entry* vcf_file::get_entry_object()
{
	return new vcf_entry(meta_data, include_indv);
}

void vcf_file::read_line(string &out)
{
	char * tmp;
	out = "";
	if (!compressed)
	{
		getline(*file_in, out);
		out.erase( out.find_last_not_of(" \t\n\r") + 1);	// Trim whitespace at end of line
	}
	else
	{
		bool again = true;
		while (again == true)
		{
			tmp = gzgets(gzfile_in, gz_readbuffer, gzMAX_LINE_LEN);
			if (tmp == NULL)
				return;

			out.append(gz_readbuffer);
			if (strlen(gz_readbuffer) != gzMAX_LINE_LEN-1)
				again = false;
		}
		out.erase( out.find_last_not_of(" \t\n\r") + 1);	// Trim whitespace at end of line (required in gzipped case!)
	}
}

void vcf_file::read_line(vector<char> &out)
{
	static string tmp;
	tmp="";
	out.resize(0);
	read_line(tmp);
	vector<char> tmp_char(tmp.begin(),tmp.end());
	out = tmp_char;
}
