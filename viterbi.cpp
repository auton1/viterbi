/*
 * epihapgen.cpp
 */

#include "viterbi.h"

output_log LOG;

void set_genetic_map_at_positions(vector<unsigned int> &positions, double rate, vector<double> &map_out)
{
	map_out.resize(positions.size());
	map_out[0] = 0;
	for (unsigned int ui=1; ui<map_out.size(); ui++)
	{
		map_out[ui] = map_out[ui-1] + rate * (double(positions[ui] - positions[ui-1]) / 1000000);
	}
}

void read_genetic_map_at_positions(string &filename, vector<unsigned int> &positions, vector<double> &map_out)
{
	// FILE FORMAT IS CHROM POS(BP) RATE(cM/Mb) MAP(cM)
	ifstream in(filename.c_str());
	if (!in.good())
		LOG.error("Could not open map file: " + filename);
	string line;
	stringstream ss;
	getline(in, line);	// Header
	string CHROM;
	int POS;
	double rate, map_pos;
	vector<int> pos;
	vector<double> map;

	while (!in.eof())
	{
		getline(in, line);
		if (line == "")
			continue;
		ss.clear(); ss.str("");
		ss.str(line);

		ss >> CHROM >> POS >> rate >> map_pos;

		pos.push_back(POS);
		map.push_back(map_pos);
	}
	// Interpolate map at positions
	map_out.resize(positions.size());

	for (unsigned int ui=0; ui<positions.size(); ui++)
	{
		vector<int>::iterator it = upper_bound(pos.begin(), pos.end(), positions[ui]);
		if (it == pos.end())
		{
			map_out[ui] = map[map.size()-1];
			continue;
		}
		else if (it != pos.begin())
			it--;

		int idx = distance(pos.begin(), it);
		int idx2 = idx+1;

		if (pos[idx] == positions[ui])
			map_out[ui] = map[idx];
		else
			map_out[ui] = map[idx] + ((map[idx2] - map[idx]) * (positions[ui] - pos[idx]) / (pos[idx2] - pos[idx]));
	}
}

void inner_loop(const vector<double> &v, double p_rec, double p_no_rec,
		const vector<vector<bool> > &haps, double em_match, double em_mismatch,
		int i, int start, int end, const vector<int> &v_idx_to_hap_idx,
		bool allele,
		vector<vector<unsigned short> > &ptr, vector<double> &vb)
{
	vector<double> va(v.size());
	vector<double>::iterator it;

	for (int j=start; j < end; j++)
	{
		transform(v.begin(), v.end(), va.begin(), bind2nd(std::plus<double>(), p_rec));
		va[j] = v[j] + p_no_rec;
		it = max_element(va.begin(), va.end());
		ptr[i-1][j] = distance(va.begin(), it);

		if (haps[ v_idx_to_hap_idx[j] ][i] == allele)
		{	// Match
			vb[j] = em_match + *it;
		}
		else
		{	// Mismatch
			vb[j] = em_mismatch + *it;
		}
	}
}

double diff(timespec start, timespec end)
{
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp.tv_sec + (temp.tv_nsec / 1000000000.0);
}

void viterbi(const vector<vector<bool> > &haps, const vector<double> &map, int indv_idx, const parameters &params, vector<int> &out_path)
{
	int N_hap = haps.size()-1;	// The minus 1 accounts for the current haplotype.
	int N_pos = haps[0].size();
	vector<double> v(N_hap, 0.0);
	//vector<double> va(N_hap, -numeric_limits<double>::max());
	vector<double> vb(N_hap, -numeric_limits<double>::max());
	vector<double>::iterator it;
	vector<vector<unsigned short> > ptr(N_pos-1, vector<unsigned short>(N_hap, 0));	// Lots of memory used here

	double em_match = log(1.0 - params.p_error);
	double em_mismatch = max(log(params.p_error), -9999.9);

	vector<int> v_idx_to_hap_idx(N_hap);
	int idx_i=0;
	for (int i=0; i<haps.size(); i++)
	{
		if (i == indv_idx)
			continue;
		v_idx_to_hap_idx[idx_i] = i;

		if (haps[indv_idx][0] == haps[i][0])
			v[idx_i] = em_match;
		else
			v[idx_i] = em_mismatch;

		idx_i++;
	}

	int numthreads = min(min((int)thread::hardware_concurrency(), params.max_threads), 8);
	LOG.printLOG("Using " + output_log::int2str(numthreads) + " threads\n");
	int rows = N_hap / numthreads;
	int extra = N_hap % numthreads;
	int output_count = 0;

	struct timespec start_time;
	struct timespec end_time;
	double prev_elapsed = numeric_limits<double>::max();
	bool increase = true;

	#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
		clock_serv_t cclock;
		mach_timespec_t mts;
		host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
		clock_get_time(cclock, &mts);
		mach_port_deallocate(mach_task_self(), cclock);
		start_time.tv_sec = mts.tv_sec;
		start_time.tv_nsec = mts.tv_nsec;
	#else
		clock_gettime(CLOCK_REALTIME, &start_time);
	#endif

	for (int i=1; i < N_pos; i++)
	{
		bool allele = haps[indv_idx][i];
		double rho = 4.0*params.Ne*((map[i] - map[i-1])/100.0);

		double tmp = exp(-rho / N_hap);
		double p_rec = (1.0-tmp) / N_hap;
		double p_no_rec = tmp + p_rec;

		p_rec = log(p_rec);
		p_no_rec = log(p_no_rec);

		/*
		for (int j=0; j < N_hap; j++)	// Could parallelize this loop?
		{	// Fast Version
			transform(v.begin(), v.end(), va.begin(), bind2nd(std::plus<double>(), p_rec));
			va[j] = v[j] + p_no_rec;
			it = max_element(va.begin(), va.end());
			ptr[i-1][j] = distance(va.begin(), it);

			if (haps[ v_idx_to_hap_idx[j] ][i] == allele)
			{	// Match
				vb[j] = em_match + *it;
			}
			else
			{	// Mismatch
				vb[j] = em_mismatch + *it;
			}
		}
		*/

		int start = 0; int end = rows;
		vector<std::thread> workers;

		for (int t=0; t<numthreads; t++)
		{
			if (t == (numthreads-1))
				end += extra;	// Last thread does extra work

			workers.push_back(std::thread(inner_loop, ref(v), p_rec, p_no_rec, ref(haps),
					em_match, em_mismatch, i, start, end,
					ref(v_idx_to_hap_idx), allele, ref(ptr), ref(vb)));
			start = end;
			end = end + rows;
		}
		for (thread& t : workers)
			t.join();

		v = vb;	// Replace old state with new state.

		if ((i % 100) == 0)
		{
			// Adaptive threading
			#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
				host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
				clock_get_time(cclock, &mts);
				mach_port_deallocate(mach_task_self(), cclock);
				end_time.tv_sec = mts.tv_sec;
				end_time.tv_nsec = mts.tv_nsec;
			#else
				clock_gettime(CLOCK_REALTIME, &end_time);
			#endif

			double elapsed = diff(start_time, end_time);

			LOG.printLOG(output_log::dbl2str(double(i)/N_pos*100, 3) + "% " + output_log::int2str(numthreads) + " threads " + output_log::dbl2str(elapsed, 3) + " seconds\n");

			if (((elapsed > prev_elapsed) && (increase == false)) || ((elapsed < prev_elapsed) && (increase == true)))
			{	// Try increasing the number of threads
				numthreads++;
				numthreads = min(numthreads, min((int)thread::hardware_concurrency(), params.max_threads));
				rows = N_hap / numthreads;
				extra = N_hap % numthreads;
				increase = true;
			}
			else if (((elapsed > prev_elapsed) && (increase == true)) || ((elapsed < prev_elapsed) && (increase == false)))
			{	// Try decreasing the number of threads
				numthreads -= 2;
				numthreads = max(numthreads, 1);
				rows = N_hap / numthreads;
				extra = N_hap % numthreads;
				increase = false;
			}
			prev_elapsed = elapsed;

			#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
				host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
				clock_get_time(cclock, &mts);
				mach_port_deallocate(mach_task_self(), cclock);
				start_time.tv_sec = mts.tv_sec;
				start_time.tv_nsec = mts.tv_nsec;
			#else
				clock_gettime(CLOCK_REALTIME, &start_time);
			#endif
		}
	}

	double max_v = -numeric_limits<double>::max();
	int max_v_idx = -1;
	for (int j=0; j<N_hap; j++)
	{
		if (v[j] > max_v)
		{
			max_v = v[j];
			max_v_idx = j;
		}
	}

	LOG.printLOG("Best Path Prob: " + output_log::dbl2str(max_v,6) + "\n");

	out_path.resize(N_pos);
	out_path[N_pos-1] = v_idx_to_hap_idx[max_v_idx];
	for (int i=N_pos-2; i >= 0; i--)
	{
		max_v_idx = ptr[i][max_v_idx];
		out_path[i] = v_idx_to_hap_idx[max_v_idx];
	}
}

void get_haplotype_count(const vector<vector<bool> > &haps, int indv_idx, int start_idx, int end_idx, int &N_exact_match, int &N_one_prct_mismatch)
{
	N_exact_match = 1;	// Obviously the haplotype itself should count towards the total!
	N_one_prct_mismatch = 1;

	int N_hap = haps.size();
	for (int i=0; i < N_hap; i++)
	{
		if (i == indv_idx)
			continue;

		int N_mismatch = 0;
		double frac_mismatch = 0;
		for (int j=start_idx; j<=end_idx; j++)
		{
			if (haps[indv_idx][j] != haps[i][j])
			{
				N_mismatch++;
				frac_mismatch = ((double)N_mismatch) / (end_idx - start_idx + 1.0);
				if (frac_mismatch > 0.01)
					break;
			}
		}

		if (N_mismatch == 0)
			N_exact_match++;
		if (frac_mismatch <= 0.01)
			N_one_prct_mismatch++;
	}
}

int main(int argc, char *argv[])
{
	time_t start,end;
	time(&start);

	// The following turns off sync between C and C++ streams.
	// Apparently it's faster to turn sync off, and as I don't use C streams, it's okay to turn off.
	ios_base::sync_with_stdio(false);

	parameters params(argc, argv);
	params.print_help();
	params.read_parameters();
	params.max_alleles = 2;	// Only use biallelic sites

	LOG.open(params.stream_out, params.stream_err, params.output_prefix);
	LOG.printLOG("\nViterbi - " + VCFTOOLS_VERSION + "\n");
	LOG.printLOG("(C) Adam Auton 2014\n\n");

	params.print_params();

	vector<vector<bool> > input_haplotypes;
	vector<unsigned int> positions;
	vector<string> alleles;
	vector<string> haplotype_names;
	
	// Read input haplotypes from VCF
	for (unsigned int ui=0; ui<params.vcf_filenames.size(); ui++)
	{
		LOG.printLOG("Reading " + params.vcf_filenames[ui] + "\n");
		variant_file *vf;
		if (params.vcf_format[ui])
			vf = new vcf_file(params, params.vcf_filenames[ui], ui);
		else
			vf = new bcf_file(params, params.vcf_filenames[ui], ui);

		vf->apply_filters(params);
		LOG.printLOG("After filtering, kept " + output_log::int2str(vf->N_kept_individuals()) + " out of " + output_log::int2str(vf->meta_data.N_indv) + " Individuals\n");

		if (ui==0)
		{
			vf->read_file(params, input_haplotypes, positions, haplotype_names, alleles);
			LOG.printLOG("After filtering, kept " + header::int2str(positions.size()) + " out of a possible " + header::int2str(vf->N_total_sites()) + " Sites\n");
		}
		else
		{
			int passed_sites = vf->add_file(params, input_haplotypes, positions, haplotype_names, alleles);
			LOG.printLOG("After filtering, kept " + header::int2str(passed_sites) + " out of a possible " + header::int2str(vf->N_total_sites()) + " Sites\n");
		}

		if (vf->N_total_sites() <= 0)
			LOG.warning("File does not contain any sites");
		else if (vf->N_kept_sites() <= 0)
			LOG.warning("No data left for analysis!");
		delete vf;
	}
	
	vector<vector<bool> > new_haplotypes;
	vector<unsigned int> new_positions;
	new_haplotypes.resize(haplotype_names.size());
	
	for (unsigned int ui=0; ui<positions.size(); ui++)
	{
		int allele_count = 0;
		for (unsigned int uj=0; uj<haplotype_names.size(); uj++)
			if (input_haplotypes[uj][ui])
				allele_count++;

		if (allele_count > 1 && allele_count < (haplotype_names.size()-1))
		{
			new_positions.push_back(positions[ui]);

			for (unsigned int uj=0; uj<haplotype_names.size(); uj++)
				new_haplotypes[uj].push_back(input_haplotypes[uj][ui]);
		}
		else
			LOG.one_off_warning("\tWarning: Skipping singletons.");
	}
	
	positions = new_positions;
	input_haplotypes = new_haplotypes;

	if (params.vcf_filenames.size() > 1)
		LOG.printLOG("Kept "+ header::int2str(positions.size()) + " sites that are contained in all files\n");
	if (positions.size() == 0)
		LOG.error("No data left for analysis!");
	
	string CHROM = *(params.chrs_to_keep.begin());

	vector<int> alt_allele_counts;
	alt_allele_counts.resize(input_haplotypes[0].size(), 0);
	for (unsigned int uj=0; uj<input_haplotypes[0].size(); uj++)
	{
		for (unsigned int uk=0; uk<input_haplotypes.size(); uk++)
			alt_allele_counts[uj] += input_haplotypes[uk][uj];
//		cout << positions[uj] << " " << alt_allele_counts[uj] << endl;
	}

	vector<double> map;
	if (params.map_filename == "")
	{
		set_genetic_map_at_positions(positions, params.recomb_rate, map);
	}
	else
	{
		read_genetic_map_at_positions(params.map_filename, positions, map);
	}

	cout << "#HAP\tCOPY\tCHR\tSTART\tEND\tN_SNPS\tMIN_ALLELE_COUNT\tMIN_ALLELE_POS\tEXACT_HAP_COUNT\tHAP_COUNT_1_PRCT_MISMATCH" << endl;
	for (int i=0; i<input_haplotypes.size(); i++)
	{
		string indv = haplotype_names[i].substr(0, haplotype_names[i].size()-2);
		if ((params.viterbi_indv.size() > 0) && (params.viterbi_indv.find(indv) == params.viterbi_indv.end()))
		{
			continue;
		}

		LOG.printLOG("Running Viterbi on " + haplotype_names[i] + "\n");
		vector<int> viterbi_path;
		viterbi(input_haplotypes, map, i, params, viterbi_path);

		// Process the viterbi path
		int N_exact_match;
		int N_one_prct_mismatch;
		int min_pos = positions[0];
		int max_pos = positions[0];
		int start_idx = 0;
		int current_haplotype = viterbi_path[0];
		int N_SNPs = 1;
		int min_freq = numeric_limits<int>::max();
		int min_freq_idx = 0;
		if (input_haplotypes[i][0] == 1)
			min_freq = alt_allele_counts[0];
		else
			min_freq = (int)input_haplotypes.size() - alt_allele_counts[0];

		for (int j=1; j<viterbi_path.size(); j++)
		{
			if (viterbi_path[j] != current_haplotype)
			{	// Output haplotype

				// Get haplotype frequency
				get_haplotype_count(input_haplotypes, i, start_idx, j-1, N_exact_match, N_one_prct_mismatch);

				cout << haplotype_names[i] << "\t" << haplotype_names[current_haplotype];
				cout << "\t" << CHROM;
				cout << "\t" << min_pos << "\t" << max_pos << "\t" << N_SNPs;
				cout << "\t" << min_freq << "\t" << positions[min_freq_idx];
				cout << "\t" << N_exact_match << "\t" << N_one_prct_mismatch << endl;
				// Reset
				min_pos = positions[j];	// Haplotype is left shifted
				max_pos = positions[j];
				start_idx = j;
				current_haplotype = viterbi_path[j];
				N_SNPs = 1;
				min_freq = numeric_limits<int>::max();
				min_freq_idx = j;
				if (input_haplotypes[i][j] == 1)
					min_freq = alt_allele_counts[j];
				else
					min_freq = (int)input_haplotypes.size() - alt_allele_counts[j];
			}
			else
			{
				max_pos = positions[j];
				N_SNPs++;
				if (input_haplotypes[i][j] == 1)
				{
					if (alt_allele_counts[j] < min_freq)
					{
						min_freq = alt_allele_counts[j];
						min_freq_idx = j;
					}
				}
				else
				{
					if (((int)input_haplotypes.size() - alt_allele_counts[j]) < min_freq)
					{
						min_freq = (int)input_haplotypes.size() - alt_allele_counts[j];
						min_freq_idx = j;
					}
				}

			}
		}
		// Output last haplotype
		get_haplotype_count(input_haplotypes, i, start_idx, viterbi_path.size()-1, N_exact_match, N_one_prct_mismatch);
		cout << haplotype_names[i] << "\t" << haplotype_names[viterbi_path[viterbi_path.size()-1]];
		cout << "\t" << CHROM;
		cout << "\t" << min_pos << "\t" << max_pos << "\t" << N_SNPs;
		cout << "\t" << min_freq << "\t" << positions[min_freq_idx];
		cout << "\t" << N_exact_match << "\t" << N_one_prct_mismatch << endl;
	}

	time(&end);
	double running_time = difftime(end,start);
	LOG.printLOG("Run Time = " + output_log::dbl2str_fixed(running_time, 2) + " seconds\n");
	LOG.close();
	return 0;
}



