/*
 * parameters.cpp
 *
 *  Created on: Nov 11, 2009
 *      Author: Adam Auton
 *      ($Revision: 249 $)
 */

// Class for reading in, checking and storing user parameters
#include "parameters.h"

parameters::parameters(int argc, char *argv[])
{
	if (isatty(STDERR_FILENO))
		stream_err = false;
	else
		stream_err = true;

	string tmp;
	for (int i=0; i<argc; i++)
	{
		tmp = argv[i];
		this->argv.push_back(tmp);
	}

	BED_exclude = false;
	BED_file = "";
	contigs_file = "";
	end_pos = numeric_limits<int>::max();
	exclude_positions_file = "";
	exclude_positions_overlap_file = "";
	indv_exclude_file = "";
	indv_keep_file = "";
	invert_mask = false;
	keep_only_indels = false;
	min_mac = -1;
	min_maf = -1.0;
	mask_file = "";
	max_alleles = numeric_limits<int>::max();
	max_genotype_depth = numeric_limits<int>::max();
	max_mac = numeric_limits<int>::max();
	max_maf = numeric_limits<double>::max();
	max_mean_depth = numeric_limits<double>::max();
	max_missing_call_count = numeric_limits<int>::max();
	max_non_ref_ac = numeric_limits<int>::max();
	max_non_ref_af = numeric_limits<double>::max();
	max_N_indv = -1;
	min_alleles = -1;
	min_genotype_depth = -1;
	min_genotype_quality = -1.0;
	min_HWE_pvalue = -1.0;
	min_interSNP_distance = -1;
	min_kept_mask_value = 0;
	min_mean_depth = -1.0;
	min_quality = -1.0;
	min_r2 = -1.0;
	min_site_call_rate = 0;
	min_non_ref_ac = -1;
	min_non_ref_af = -1.0;
	num_outputs = 0;
	output_prefix = "out";
	phased_only = false;
	positions_file = "";
	positions_overlap_file = "";
	remove_all_filtered_genotypes = false;
	remove_all_filtered_sites = false;
	remove_indels = false;
	snps_to_exclude_file = "";
	snps_to_keep_file = "";
	start_pos = -1;
	stream_in = false;
	stream_out = false;
	temp_dir = "/tmp/";

	Ne = 10000;
	recomb_rate = 1.0;	// cM/Mb
	map_filename = "";
	p_error = 0.0;
	max_threads = 64;

	time_t start;
	time(&start);
	seed = start;
	generator.seed(seed);

	run_viterbi = true;
}

void parameters::read_parameters()
{
	unsigned int i=1;
	string in_str;
	while (i<argv.size())
	{
		in_str = argv[i];
		if (in_str == "--vcf") // VCF file to process
		{
			vcf_format.push_back( true );
			vcf_compressed.push_back( false );
			if (!stream_in)
			{
				vcf_filenames.push_back(get_arg(i+1));
				i++;
			}
		}
		else if (in_str == "--bcf") // BCF file to process
		{
			vcf_format.push_back( false );
			vcf_compressed.push_back( false );
			if (!stream_in)
			{
				vcf_filenames.push_back(get_arg(i+1));
				i++;
			}
		}
		else if (in_str == "--bed") {
			if (BED_file == "")
			{
				BED_file = get_arg(i+1); i++; BED_exclude=false;
			}
			else
				LOG.error(" Multiple --bed/--exclude-bed options can not be used together.");
		}
		else if (in_str == "-c") {stream_out = true;}						// Write output to stream
		else if (in_str == "--chr") { chrs_to_keep.insert(get_arg(i+1)); i++; }					// Chromosome to process
		else if (in_str == "--contigs") { contigs_file = get_arg(i+1); i++;}	// Contigs file for header
		else if (in_str == "--exclude-bed") {
			if (BED_file == "")
			{
				BED_file = get_arg(i+1); i++; BED_exclude=true;
			}
			else
				LOG.error(" Multiple --bed/--exclude-bed options can not be used together.");
		}
		else if (in_str == "--exclude") { snps_to_exclude_file = get_arg(i+1); i++; }				// List of SNPs to exclude
		else if (in_str == "--exclude-positions") { exclude_positions_file = get_arg(i+1); i++; }
		else if (in_str == "--exclude-positions-overlap") { exclude_positions_overlap_file = get_arg(i+1); i++; }
		else if (in_str == "--from-bp") { start_pos = atoi(get_arg(i+1).c_str()); i++; }					// Start position
		else if (in_str == "--get-INFO") {
			if (INFO_to_extract.empty())
				num_outputs++;
			INFO_to_extract.push_back(get_arg(i+1)); i++;}	// Add to list of INFO fields to extract
		else if (in_str == "--gzvcf") // Compressed VCF file to process
		{
			vcf_format.push_back( true );
			vcf_compressed.push_back( true );
			if (!stream_in)
			{
				vcf_filenames.push_back(get_arg(i+1));
				i++;
			}
		}
		else if (in_str == "--hwe") { max_alleles = 2; min_HWE_pvalue = atof(get_arg(i+1).c_str()); i++; }					// Minimum per-site HWE p-value
		else if (in_str == "--indv") { indv_to_keep.insert(get_arg(i+1)); i++; }						// List of individuals to keep
		else if (in_str == "--invert-mask") { mask_file = get_arg(i+1); i++; invert_mask = true; }
		else if (in_str == "--keep-filtered") { site_filter_flags_to_keep.insert(get_arg(i+1)); i++; }	// Remove a specific filter flag
		else if (in_str == "--keep") { indv_keep_file = get_arg(i+1); i++; }						// List of individuals to keep
		else if (in_str == "--keep-only-indels") { keep_only_indels = true; }
		else if (in_str == "--keep-INFO") { site_INFO_flags_to_keep.insert(get_arg(i+1)); i++; }	// Filter sites by INFO flags
		else if (in_str == "--mac") { min_mac = atoi(get_arg(i+1).c_str()); i++; }								// Minimum Site MAC
		else if (in_str == "--maf") { min_maf = atof(get_arg(i+1).c_str()); i++; }								// Minimum Site MAF
		else if (in_str == "--mask-min") { min_kept_mask_value = atoi(get_arg(i+1).c_str()); i++; }
		else if (in_str == "--mask") { mask_file = get_arg(i+1); i++; invert_mask = false; }
		else if (in_str == "--max-alleles") { max_alleles = atoi(get_arg(i+1).c_str()); i++; }				// Maximum number of alleles per-site
		else if (in_str == "--max-mac") { max_mac = atoi(get_arg(i+1).c_str()); i++; }						// Maximum site MAC
		else if (in_str == "--max-maf") { max_maf = atof(get_arg(i+1).c_str()); i++; }						// Maximum Site MAF
		else if (in_str == "--max-meanDP") { max_mean_depth = atof(get_arg(i+1).c_str()); i++; }			// Site Maximum mean depth across individuals
		else if (in_str == "--max-missing") { min_site_call_rate = atof(get_arg(i+1).c_str()); i++; }
		else if (in_str == "--max-missing-count") { max_missing_call_count = atoi(get_arg(i+1).c_str()); i++; } // Site maximum missing genotypes
		else if (in_str == "--max-non-ref-ac") { max_non_ref_ac = atoi(get_arg(i+1).c_str()); i++; }		// Minimum Site non-ref AC
		else if (in_str == "--max-non-ref-af") { max_non_ref_af = atof(get_arg(i+1).c_str()); i++; }		// Minimum Site non-ref AF
		else if (in_str == "--maxDP") { max_genotype_depth = atoi(get_arg(i+1).c_str()); i++; }				// Maximum genotype depth
		else if (in_str == "--max-indv") {max_N_indv = atoi(get_arg(i+1).c_str()); i++; }
		else if (in_str == "--min-alleles") { min_alleles = atoi(get_arg(i+1).c_str()); i++; }				// Minimum number of alleles per-site
		else if (in_str == "--min-meanDP") { min_mean_depth = atof(get_arg(i+1).c_str()); i++; }			// Site Minimum mean depth
		else if (in_str == "--min-r2") { min_r2 = atof(get_arg(i+1).c_str()); i++; }					// Min r^2 for LD output
		else if (in_str == "--minDP") { min_genotype_depth = atoi(get_arg(i+1).c_str()); i++; }				// Minimum genotype depth
		else if (in_str == "--minGQ") { min_genotype_quality = atof(get_arg(i+1).c_str()); i++; }			// Minimum genotype quality
		else if (in_str == "--minQ") { min_quality = atof(get_arg(i+1).c_str()); i++; }						// Minimum per-site quality
		else if (in_str == "--non-ref-ac") { min_non_ref_ac = atoi(get_arg(i+1).c_str()); i++; }				// Minimum Site non-ref AC
		else if (in_str == "--non-ref-af") { min_non_ref_af = atof(get_arg(i+1).c_str()); i++; }				// Minimum Site non-ref AF
		else if (in_str == "--not-chr") { chrs_to_exclude.insert(get_arg(i+1)); i++; }					// Chromosome to process
		else if (in_str == "--out") { output_prefix = get_arg(i+1); i++; }							// Output file prefix
		else if (in_str == "--phased") phased_only = true;								// Keep only phased individuals / sites
		else if (in_str == "--positions") { positions_file = get_arg(i+1); i++; }
		else if (in_str == "--positions-overlap") { positions_overlap_file = get_arg(i+1); i++; }
		else if (in_str == "--remove-filtered-all") remove_all_filtered_sites = true;							// Remove sites flagged as filtered
		else if (in_str == "--remove-filtered-geno-all") remove_all_filtered_genotypes = true;			// Remove genotypes flagged as filtered
		else if (in_str == "--remove-filtered-geno") { geno_filter_flags_to_exclude.insert(get_arg(i+1)); i++; }		// Remove genotypes flagged as filtered
		else if (in_str == "--remove-filtered") { site_filter_flags_to_exclude.insert(get_arg(i+1)); i++; }	// Remove a specific filter flag
		else if (in_str == "--remove-indels") { remove_indels = true; }
		else if (in_str == "--remove-indv") { indv_to_exclude.insert(get_arg(i+1)); i++; }						// List of individuals to keep
		else if (in_str == "--remove-INFO") { site_INFO_flags_to_remove.insert(get_arg(i+1)); i++; }	// Filter sites by INFO flags
		else if (in_str == "--remove") { indv_exclude_file = get_arg(i+1); i++; }					// List of individuals to exclude
		else if (in_str == "--seed") { seed = atoi(get_arg(i+1).c_str()); generator.seed(seed); i++; }
		else if (in_str == "--snp") { snps_to_keep.insert(get_arg(i+1)); i++; }						// SNP to keep
		else if (in_str == "--snps") { snps_to_keep_file = get_arg(i+1); i++; }						// List of SNPs to keep
		else if (in_str == "--stdout") {stream_out = true; }						// Write output to stream
		else if (in_str == "--temp") { temp_dir = get_arg(i+1); i++;}	// Directory for vcftools temporary files
		else if (in_str == "--to-bp") { end_pos = atoi(get_arg(i+1).c_str()); i++; }						// End position
		else if (in_str == "--thin") { min_interSNP_distance = atoi(get_arg(i+1).c_str()); i++; }
		else if (in_str == "--Ne") { Ne = atof(get_arg(i+1).c_str()); i++; }
		else if (in_str == "--map") { map_filename = get_arg(i+1); i++; }
		else if (in_str == "--recomb") { recomb_rate = atof(get_arg(i+1).c_str()); i++; }
		else if (in_str == "--test-indv") { viterbi_indv.insert(get_arg(i+1)); i++; }
		else if (in_str == "--error") { p_error = atof(get_arg(i+1).c_str()); i++; }
		else if (in_str == "--maxthreads") { max_threads = atoi(get_arg(i+1).c_str()); i++; }
		else if (in_str == "--fwdbck") { run_viterbi = false; }
		else
			error("Unknown option: " + string(in_str), 0);
		i++;
	}
	check_parameters();
}

string parameters::get_arg(unsigned int i)
{
	if (i>=argv.size())
		error("Requested Missing Argument",76);
	return argv[i];
}

void parameters::print_params()
{
	parameters defaults(0, 0);

	LOG.printLOG("Parameters as interpreted:\n");

	for (unsigned int ui=0; ui<vcf_filenames.size(); ui++)
	{
		string tmp_name = vcf_filenames[ui];
		if (vcf_format[ui] == false)
			LOG.printLOG("\t--bcf " + tmp_name + "\n");
		else if (vcf_format[ui] == true && vcf_compressed[ui] == false)
			LOG.printLOG("\t--vcf " + tmp_name + "\n");
		else if (vcf_format[ui] == true && vcf_compressed[ui] == true)
			LOG.printLOG("\t--gzvcf " + tmp_name + "\n");
	}

	if (chrs_to_keep.size() > 0)
	{
		for (set<string>::iterator it=chrs_to_keep.begin(); it != chrs_to_keep.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--chr " + tmp + "\n");
		}
	}

	if (chrs_to_exclude.size() > 0)
	{
		for (set<string>::iterator it=chrs_to_exclude.begin(); it != chrs_to_exclude.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--not-chr " + tmp + "\n");
		}
	}
	if (contigs_file != defaults.contigs_file) LOG.printLOG("\t--contigs " + contigs_file + "\n");
	if (end_pos != defaults.end_pos) LOG.printLOG("\t--to-bp " + output_log::int2str(end_pos) + "\n");
	if (exclude_positions_file != defaults.exclude_positions_file) LOG.printLOG("\t--exclude-positions " + exclude_positions_file + "\n");
	if (exclude_positions_overlap_file != defaults.exclude_positions_overlap_file) LOG.printLOG("\t--exclude-positions-overlap " + exclude_positions_overlap_file + "\n");
	if (indv_exclude_file != defaults.indv_exclude_file) LOG.printLOG("\t--exclude " + indv_exclude_file + "\n");
	if (indv_keep_file != defaults.indv_keep_file) LOG.printLOG("\t--keep " + indv_keep_file + "\n");
	if (keep_only_indels != defaults.keep_only_indels) LOG.printLOG("\t--keep-only-indels\n");
	if (min_mac != defaults.min_mac) LOG.printLOG("\t--mac " + output_log::dbl2str(min_mac, 3) + "\n");
	if (min_maf != defaults.min_maf) LOG.printLOG("\t--maf " + output_log::dbl2str(min_maf, 3) + "\n");
	if (max_alleles != defaults.max_alleles) LOG.printLOG("\t--max-alleles " + output_log::int2str(max_alleles) + "\n");
	if (max_genotype_depth != defaults.max_genotype_depth) LOG.printLOG("\t--maxDP " + output_log::dbl2str(max_genotype_depth, 3) + "\n");
	if (max_mac != defaults.max_mac) LOG.printLOG("\t--max-mac " + output_log::dbl2str(max_mac, 3) + "\n");
	if (max_maf != defaults.max_maf) LOG.printLOG("\t--max-maf " + output_log::dbl2str(max_maf, 3) + "\n");
	if (max_missing_call_count != defaults.max_missing_call_count) LOG.printLOG("\t--max-missing-count " + output_log::dbl2str(max_missing_call_count, 3) + "\n");
	if (max_mean_depth != defaults.max_mean_depth) LOG.printLOG("\t--max-meanDP " + output_log::dbl2str(max_mean_depth, 3) + "\n");
	if (max_non_ref_ac != defaults.max_non_ref_ac) LOG.printLOG("\t--max-non-ref-ac " + output_log::dbl2str(max_non_ref_ac, 3) + "\n");
	if (max_non_ref_af != defaults.max_non_ref_af) LOG.printLOG("\t--max-non-ref-af " + output_log::dbl2str(max_non_ref_af, 3) + "\n");
	if (max_N_indv != defaults.max_N_indv) LOG.printLOG("\t--max-indv " + output_log::int2str(max_N_indv) + "\n");
	if (min_alleles != defaults.min_alleles) LOG.printLOG("\t--min-alleles " + output_log::int2str(min_alleles) + "\n");
	if (min_genotype_depth != defaults.min_genotype_depth) LOG.printLOG("\t--minDP " + output_log::dbl2str(min_genotype_depth, 3) + "\n");
	if (min_genotype_quality != defaults.min_genotype_quality) LOG.printLOG("\t--minGQ " + output_log::dbl2str(min_genotype_quality, 3) + "\n");
	if (min_HWE_pvalue != defaults.min_HWE_pvalue) LOG.printLOG("\t--hwe " + output_log::dbl2str(min_HWE_pvalue, 3) + "\n");
	if (min_interSNP_distance != defaults.min_interSNP_distance) LOG.printLOG("\t--thin " + output_log::int2str(min_interSNP_distance) + "\n");
	if (min_kept_mask_value != defaults.min_kept_mask_value) LOG.printLOG("\t--mask-min " + output_log::int2str(min_kept_mask_value) + "\n");
	if (min_mean_depth != defaults.min_mean_depth) LOG.printLOG("\t--min-meanDP " + output_log::dbl2str(min_mean_depth, 3) + "\n");
	if (min_quality != defaults.min_quality) LOG.printLOG("\t--minQ " + output_log::dbl2str(min_quality, 3) + "\n");
	if (min_r2 != defaults.min_r2) LOG.printLOG("\t--min-r2 " + output_log::dbl2str(min_r2, 3) + "\n");
	if (min_site_call_rate != defaults.min_site_call_rate) LOG.printLOG("\t--max-missing " + output_log::dbl2str(min_site_call_rate, 3) + "\n");
	if (min_non_ref_ac != defaults.min_non_ref_ac) LOG.printLOG("\t--non-ref-ac " + output_log::dbl2str(min_non_ref_ac, 3) + "\n");
	if (min_non_ref_af != defaults.min_non_ref_af) LOG.printLOG("\t--non-ref-af " + output_log::dbl2str(min_non_ref_af, 3) + "\n");
	if (output_prefix != defaults.output_prefix) LOG.printLOG("\t--out " + output_prefix + "\n");
	if (phased_only) LOG.printLOG("\t--phased\n");
	if (positions_file != defaults.positions_file) LOG.printLOG("\t--positions " + positions_file + "\n");
	if (positions_overlap_file != defaults.positions_overlap_file) LOG.printLOG("\t--positions-overlap " + positions_overlap_file + "\n");
	if (remove_all_filtered_genotypes) LOG.printLOG("\t--remove-filtered-geno-all\n");
	if (remove_all_filtered_sites) LOG.printLOG("\t--remove-filtered-all\n");
	if (remove_indels != defaults.remove_indels) LOG.printLOG("\t--remove-indels\n");
	if (snps_to_exclude_file != defaults.snps_to_exclude_file) LOG.printLOG("\t--exclude " + snps_to_exclude_file + "\n");
	if (snps_to_keep_file != defaults.snps_to_keep_file) LOG.printLOG("\t--snps " + snps_to_keep_file + "\n");
	if (start_pos != defaults.start_pos) LOG.printLOG("\t--from-bp " + output_log::int2str(start_pos) + "\n");
	if (stream_out != defaults.stream_out) LOG.printLOG("\t--stdout\n");
	if (temp_dir != defaults.temp_dir) LOG.printLOG("\t--temp " + temp_dir + "\n");

	if (site_filter_flags_to_exclude.size() > 0)
		for (set<string>::iterator it=site_filter_flags_to_exclude.begin(); it != site_filter_flags_to_exclude.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--remove-filtered " + tmp + "\n");
		}

	if (site_filter_flags_to_keep.size() > 0)
		for (set<string>::iterator it=site_filter_flags_to_keep.begin(); it != site_filter_flags_to_keep.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--keep-filtered " + tmp + "\n");
		}

	if (geno_filter_flags_to_exclude.size() > 0)
		for (set<string>::iterator it=geno_filter_flags_to_exclude.begin(); it != geno_filter_flags_to_exclude.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--remove-filtered-geno " + tmp + "\n");
		}

	if (site_INFO_flags_to_remove.size() > 0)
		for (set<string>::iterator it=site_INFO_flags_to_remove.begin(); it != site_INFO_flags_to_remove.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--remove-INFO " + tmp + "\n");
		}

	if (site_INFO_flags_to_keep.size() > 0)
		for (set<string>::iterator it=site_INFO_flags_to_keep.begin(); it != site_INFO_flags_to_keep.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--keep-INFO " + tmp + "\n");
		}

	if (BED_file != defaults.BED_file)
	{
		if (BED_exclude == false)
			LOG.printLOG("\t--bed " + BED_file + "\n");
		else
			LOG.printLOG("\t--exclude-bed " + BED_file + "\n");
	}

	if (mask_file != defaults.mask_file)
	{
		if (invert_mask == false)
			LOG.printLOG("\t--mask " + mask_file + "\n");
		else
			LOG.printLOG("\t--invert-mask " + mask_file + "\n");
	}

	if (snps_to_keep.size() > 0)
		for (set<string>::iterator it=snps_to_keep.begin(); it != snps_to_keep.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--snp " + tmp + "\n");
		}

	if (indv_to_keep.size() > 0)
		for (set<string>::iterator it=indv_to_keep.begin(); it != indv_to_keep.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--indv " + tmp + "\n");
		}

	if (indv_to_exclude.size() > 0)
		for (set<string>::iterator it=indv_to_exclude.begin(); it != indv_to_exclude.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--remove-indv " + tmp + "\n");
		}

	LOG.printLOG("\t--seed " + output_log::int2str(seed) + "\n");
	LOG.printLOG("\t--maxthreads " + output_log::int2str(max_threads) + "\n");

	if (run_viterbi != defaults.run_viterbi)
		LOG.printLOG("\t--fwdbck\n");

	if (map_filename != defaults.map_filename) LOG.printLOG("\t--map " + map_filename + "\n");
	if (recomb_rate != defaults.recomb_rate) LOG.printLOG("\t--recomb " + output_log::dbl2str(recomb_rate, 3) + "\n");

	if (viterbi_indv != defaults.viterbi_indv)
	{
		for (auto it = viterbi_indv.begin(); it != viterbi_indv.end(); ++it)
		{
			LOG.printLOG("\t--test-indv " + *it + "\n");
		}
	}

	if (p_error != defaults.p_error) LOG.printLOG("\t--error " + output_log::dbl2str(p_error, 3) + "\n");

	LOG.printLOG("\n");
}

void parameters::print_help()
{
	unsigned int i;
	string in_str;

	if (argv.size() <= 1)
	{	// If there are no user parameters, display help.
		argv.push_back("--?");
		print_help();
	}

	for(i = 0; i < argv.size(); i++)
	{
		in_str = argv[i];
		if ((in_str == "-h") || (in_str == "-?") || (in_str == "-help") || (in_str == "--?") || (in_str == "--help") || (in_str == "--h"))
		{
			cout << endl << "Viterbi (" << VCFTOOLS_VERSION << ")" << endl;
			cout << "\u00A9 Adam Auton 2014" << endl << endl;
			cout << "Identify shared haplotypes using Viterbi algorithm." << endl;
			cout << endl;
			cout << endl;

			exit(0);
		}
	}
}

void parameters::check_parameters()
{
	parameters defaults(0, 0);

	if (isatty(STDIN_FILENO) && stream_in)
		LOG.error("No input detected via stream.");

	if (chrs_to_keep.size() != 1)
		LOG.error("Require a chromosome.");

	if ((p_error < 0.0) || (p_error > 1.0)) LOG.error("Error probability must be between 0 and 1.", 1);

	if (num_outputs > 1) error("Only one output function may be called.",0);
	if (chrs_to_keep.size() > 0 && chrs_to_exclude.size() > 0) error("Cannot specify chromosomes to keep and to exclude", 1);
	if (end_pos < start_pos) error("End position must be greater than Start position.", 1);
	if (((end_pos != numeric_limits<int>::max()) || (start_pos != -1)) && (chrs_to_keep.size() != 1)) error("Require a single chromosome when specifying a range.", 2);
	if (max_maf < min_maf) error("Maximum MAF must be not be less than Minimum MAF.", 4);
	if (max_mac < min_mac) error("Maximum MAC must be not be less than Minimum MAC.", 4);
	if (min_maf != defaults.min_maf)
	{
		if ((min_maf < 0.0) || (min_maf > 1.0)) error("MAF must be between 0 and 1.", 4);
	}
	if (max_maf != defaults.max_maf)
	{
		if ((max_maf < 0.0) || (max_maf > 1.0)) error("Maximum MAF must be between 0 and 1.", 4);
	}
	if (min_non_ref_af != defaults.min_non_ref_af)
	{
		if ((min_non_ref_af < 0.0) || (min_non_ref_af > 1.0)) error("Non-Ref Allele Frequency must be between 0 and 1.", 4);
	}
	if (max_non_ref_af < min_non_ref_af) error("Maximum Non-Ref Allele Frequency must not be less that Minimum Non-Ref AF.", 4);
	if (max_non_ref_ac < min_non_ref_ac) error("Maximum Non-Ref Allele Count must not be less that Minimum Non-Ref AC.", 4);
	if (min_site_call_rate > 1) error("Minimum Call rate cannot be greater than 1.", 5);
	if (max_alleles < min_alleles) error("Max Number of Alleles must be greater than Min Number of Alleles.", 6);
	if (max_mean_depth < min_mean_depth) error("Max Mean Depth must be greater the Min Mean Depth.", 7);
	if (max_genotype_depth < min_genotype_depth) error("Max Genotype Depth must be greater than Min Genotype Depth.", 9);
	if (min_kept_mask_value > 9) error("Min Mask value must be between 0 and 9.", 14);
}

void parameters::error(string err_msg, int code)
{
	LOG.printLOG("\n\nError: " + err_msg + "\n\n");
	exit(code);
}
