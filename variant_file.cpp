/*
 * variant_file.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: amarcketta
 */

#include "variant_file.h"

variant_file::~variant_file() {}

// Return the number of individuals that have not been filtered out
int variant_file::N_kept_individuals() const
{
	int N_kept = 0;
	for (unsigned int ui=0; ui<include_indv.size(); ui++)
		if (include_indv[ui] == true)
			N_kept++;
	return N_kept;
}

// Return the number of sites that have not been filtered out
int variant_file::N_kept_sites() const
{
	return N_kept_entries;
}

// Return the total number of sites in the file
int variant_file::N_total_sites() const
{
	return N_entries;
}

void variant_file::ByteSwap(unsigned char *b, int n) const
{
   register int i = 0;
   register int j = n-1;
   while (i<j)
   {
      std::swap(b[i], b[j]);
      i++, j--;
   }
}

void variant_file::get_contigs(const string &contigs_file, vector<string> &contig_vector)
{
	if (contigs_file == "")
		LOG.error("Contig declarations in header are necessary for BCF conversion. Use --contigs <filename> to add contigs to the header.");

	ifstream contigs(contigs_file.c_str());
	if (!contigs.is_open())
		LOG.error("Could not open contigs file: " + contigs_file);

	string line;
	int contig_lines = 0;
	contig_vector.resize(0);

	while (getline(contigs, line))
	{
		if (line.find("##contig=")==string::npos)
			LOG.error("Contigs file must contain only contig header lines.");

		contig_vector.push_back(line);
		contig_lines++;
	}

	contigs.close();
	LOG.printLOG("Including "+header::int2str(contig_lines)+" header lines from the contig file.\n");
}

void variant_file::read_temp_site(ifstream &tmp_file, string &CHROM, int &POS, vector< pair<int,int> > &GTs)
{
	stringstream chr;
	char tmp_char;
	while(true)
	{
		tmp_file.read(&tmp_char,sizeof(char));
		if (tmp_char == '\n')
			break;
		chr << tmp_char;
	}
	CHROM = chr.str();

	tmp_file.read((char*)&POS,sizeof(POS));

	char in_byte, tmp_gt;
	for(unsigned int ui=0; ui<GTs.size(); ui++)
	{
		tmp_file.read(&in_byte,sizeof(in_byte));
		tmp_gt = in_byte & 0x03;
		if (tmp_gt == 0x02)
			GTs[ui].second = -1;
		else
			GTs[ui].second = (int)tmp_gt;
		in_byte = in_byte >> 4;
		tmp_gt = in_byte & 0x03;
		if (tmp_gt == 0x02)
			GTs[ui].first = -1;
		else
			GTs[ui].first = (int)tmp_gt;
	}
}

void variant_file::read_big_temp_site(ifstream &tmp_file, string &CHROM, int &POS, int &alleles, vector< pair<int,int> > &GTs)
{
	stringstream chr;
	char tmp_char;
	while(true)
	{
		tmp_file.read(&tmp_char,sizeof(char));
		if (tmp_char == '\n')
			break;
		chr << tmp_char;
	}
	CHROM = chr.str();

	tmp_file.read((char*)&POS,sizeof(POS));

	int8_t tmp_alleles;
	tmp_file.read((char*)&tmp_alleles,sizeof(tmp_alleles));
	alleles = (int)tmp_alleles;

	char in_byte = 0xFF;
	for(unsigned int ui=0; ui<GTs.size(); ui++)
	{
		tmp_file.read(&in_byte,sizeof(in_byte));
		if (in_byte == (char)0xFF)
			GTs[ui].first = -1;
		else
			GTs[ui].first = (int)in_byte;

		tmp_file.read(&in_byte,sizeof(in_byte));
		if (in_byte == (char)0xFF)
			GTs[ui].second = -1;
		else
			GTs[ui].second = (int)in_byte;
	}
}

void variant_file::read_file(const parameters &params, vector<vector<bool> > &haplotypes, vector<unsigned int> &positions, vector<string> &haplotype_names, vector<string> &alleles)
{
	unsigned int N_chr = 2*N_kept_individuals();
	haplotypes.resize(N_chr);

	haplotype_names.resize(0);
	for (unsigned int ui=0; ui<N_kept_individuals(); ui++)
	{
		if (include_indv[ui] == false)
			continue;
		haplotype_names.push_back(meta_data.indv[ui] + "_1");
		haplotype_names.push_back(meta_data.indv[ui] + "_2");
	}

	vector<char> variant_line;
	entry *e = get_entry_object();
	unsigned int N_alleles;

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;
		e->parse_basic_entry(true);

		e->parse_genotype_entries(true);
		N_alleles = e->get_N_alleles();

		if (N_alleles != 2)
		{
			LOG.one_off_warning("\tWarning: Only using biallelic sites.");
			continue;
		}
		if (e->is_diploid() == false)
		{
			LOG.one_off_warning("\tWarning: Only using fully diploid sites.");
			continue;
		}

		vector<bool> tmp;
		tmp.reserve(meta_data.N_indv*2);
		int allele_count = 0;

		pair<int, int> geno;
		bool skip = false;
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if (e->include_indv[ui] == false)
				continue;

			if (e->include_genotype[ui] == true)
			{
				e->get_indv_GENOTYPE_ids(ui, geno);
				if ((geno.first > 1) || (geno.second > 1))
					LOG.error("Only supporting biallelics... something bad happened", 1);
				if ((geno.first != -1) && (geno.second != -1))
				{
					tmp.push_back((bool)geno.first);
					tmp.push_back((bool)geno.second);
					allele_count += geno.first;
					allele_count += geno.second;
				}
				else
				{
					skip = true;
					LOG.one_off_warning("\tWarning: Only using sites without missing data.");
					break;
				}
			}
			else
			{
				skip = true;
				LOG.one_off_warning("\tWarning: Only using sites without missing data.");
				break;
			}
		}

		if (skip == true)
			continue;

		if (N_chr != tmp.size())
			LOG.error("Something unexpected happened.... oh dear.");

		string chr = e->get_CHROM();

		alleles.push_back(e->get_REF());
		alleles.push_back(e->get_ALT());
		positions.push_back(e->get_POS());
		
		for (unsigned int ui=0; ui<tmp.size(); ui++)
			haplotypes[ui].push_back(tmp[ui]);
	}
	delete e;
}

int variant_file::add_file(const parameters &params, vector<vector<bool> > &haplotypes, vector<unsigned int> &positions, vector<string> &haplotype_names, vector<string> &alleles)
{
	int passed_sites = 0;
	unsigned int N_chr = 2*N_kept_individuals();
	vector<vector<bool> > new_haplotypes;
	vector<unsigned int> new_positions;
	vector<string> new_alleles;
	
	new_haplotypes.resize(N_chr);
	
	for (unsigned int ui=0; ui<N_kept_individuals(); ui++)
	{
		if (include_indv[ui] == false)
			continue;
		haplotype_names.push_back(meta_data.indv[ui] + "_1");
		haplotype_names.push_back(meta_data.indv[ui] + "_2");
	}

	vector<char> variant_line;
	entry *e = get_entry_object();
	unsigned int N_alleles;

	while(!eof())
	{
		get_entry(variant_line);
		e->reset(variant_line);
		N_entries += e->apply_filters(params);

		if(!e->passed_filters)
			continue;
		N_kept_entries++;
		e->parse_basic_entry(true);

		e->parse_genotype_entries(true);
		N_alleles = e->get_N_alleles();

		if (N_alleles != 2)
		{
			LOG.one_off_warning("\tWarning: Only using biallelic sites.");
			continue;
		}
		if (e->is_diploid() == false)
		{
			LOG.one_off_warning("\tWarning: Only using fully diploid sites.");
			continue;
		}

		vector<bool> tmp;
		tmp.reserve(meta_data.N_indv*2);
		int allele_count = 0;

		pair<int, int> geno;
		bool skip = false;
		for (unsigned int ui=0; ui<meta_data.N_indv; ui++)
		{
			if (e->include_indv[ui] == false)
				continue;

			if (e->include_genotype[ui] == true)
			{
				e->get_indv_GENOTYPE_ids(ui, geno);
				if ((geno.first > 1) || (geno.second > 1))
					LOG.error("Only supporting biallelics... something bad happened", 1);
				if ((geno.first != -1) && (geno.second != -1))
				{
					tmp.push_back((bool)geno.first);
					tmp.push_back((bool)geno.second);
					allele_count += geno.first;
					allele_count += geno.second;
				}
				else
				{
					skip = true;
					LOG.one_off_warning("\tWarning: Only using sites without missing data.");
					break;
				}
			}
			else
			{
				skip = true;
				LOG.one_off_warning("\tWarning: Only using sites without missing data.");
				break;
			}
		}

		if (skip == true)
			continue;

		if (N_chr != tmp.size())
			LOG.error("Something unexpected happened.... oh dear.");

		string chr = e->get_CHROM();
		passed_sites++;

		new_alleles.push_back(e->get_REF());
		new_alleles.push_back(e->get_ALT());
		new_positions.push_back(e->get_POS());

		for (unsigned int ui=0; ui<tmp.size(); ui++)
			new_haplotypes[ui].push_back(tmp[ui]);
		
	}
	vector<vector<bool> > merged_haplotypes;
	vector<unsigned int> merged_positions;
	vector<string> merged_alleles;
	
	merged_haplotypes.resize(haplotypes.size()+new_haplotypes.size());
	for (unsigned int ui=0; ui<positions.size(); ui++)
	{
		int idx = find(new_positions.begin(), new_positions.end(), positions[ui]) - new_positions.begin();

		if (idx >= 0 && idx < new_positions.size())
		{
			if (alleles[2*ui] == new_alleles[2*idx] && alleles[2*ui+1] == new_alleles[2*idx+1])
			{
				merged_positions.push_back(positions[ui]);
				for (unsigned int uj=0; uj<haplotypes.size(); uj++)
					merged_haplotypes[uj].push_back(haplotypes[uj][ui]);

				for (unsigned int uj=0; uj<new_haplotypes.size(); uj++)
					merged_haplotypes[haplotypes.size()+uj].push_back(new_haplotypes[uj][idx]);
			}
		}
	}
	haplotypes = merged_haplotypes;
	positions = merged_positions;
	
	delete e;
	return passed_sites;
}
