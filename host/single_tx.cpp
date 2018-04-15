#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/usrp/multi_usrp.hpp>

#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/thread.hpp>
#include <iostream>
#include <complex>
#include <cmath>

#include <csignal>
#include <stdio.h>
#include <stdlib.h>
#include "single_tx.h"

namespace po = boost::program_options;

//System parameters 
double freq, gain, thres;
double inter;
double rate;

// USRP
uhd::usrp::multi_usrp::sptr usrp;
string usrp_ip;

// Tx/Rx metadata
uhd::rx_metadata_t rx_md;
uhd::tx_metadata_t tx_md;

// Buffer
gr_complex pkt[MAX_PKT_LEN];
gr_complex zeros[SYM_LEN];

// File
FILE* in_file;
FILE* out_file;
string in_name, out_name; 

// Evaluation
size_t r_cnt;
static bool stop_signal = false;

void init_usrp() {
	usrp = uhd::usrp::multi_usrp::make(usrp_ip);
	usrp->set_rx_rate(rate);
	usrp->set_tx_rate(rate);

	usrp->set_rx_freq(freq);
	usrp->set_tx_freq(freq);

	usrp->set_rx_gain(gain);
}

void sync_clock() {
	cout << "SYNC Clock" << endl;
	usrp->set_clock_config(uhd::clock_config_t::external());
	usrp->set_time_next_pps(uhd::time_spec_t(0.0));
}


void init_sys() {
	// Buffer initialize	
	memset(pkt, 0, sizeof(pkt));
	
	// TODO: put generated signal in pkt
	while (sample_cnt < MAX_PKT_LEN) {
		pkt[sample_cnt++] = (gr_complex)(rand() & 255);
	}

	// Write into BIN file
	in_file = fopen(in_name.c_str(), "rb");
	size_t in_pkt = fread(pkt, sizeof(gr_complex), MAX_PKT_LEN, in_file);
	fprintf(stderr, "%zu bytes\n", in_pkt);
	fclose(in_file);

	/*
	for( int i = 0; i < sample_cnt; i++ )
		printf("%3d: %.4lf %.4lf\n", i, pkt[i].real(), pkt[i].imag());
	*/

	// USRP setting
	rate = 1e8 / inter;
	freq = 1e9 * freq;

	init_usrp();
}

void sig_int_handler(int) {
	stop_signal = true;
}

int UHD_SAFE_MAIN(int argc, char *argv[]) {
	uhd::set_thread_priority_safe();
	uhd::time_spec_t refer;

	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "help message")
		("r0", po::value<string>(&usrp_ip)->default_value("addr=192.168.10.2"), "usrp's IP")
		("in", po::value<string>(&in_name)->default_value("wcs_trace/tx_data.bin"), "binary samples file")
		("out", po::value<string>(&out_name)->default_value("rx_log.dat"), "signal file")
		("i", po::value< double >(&inter)->default_value(128), "interval of two sampling")
		("f", po::value< double >(&freq)->default_value(2.49), "RF center frequency in Hz")
		("g", po::value< double >(&gain)->default_value(30.0), "gain for the RF chain")
		("c", po::value< size_t >(&r_cnt)->default_value(90), "round count");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
		cout << boost::format("UHD TX samples from file %s") % desc << endl;
		return ~0;
	}

	// Init
	init_sys();

	// Setup time
	boost::this_thread::sleep(boost::posix_time::milliseconds(WARM_UP_TIME));

	std::signal(SIGINT, &sig_int_handler);
	std::cout << "Press Ctrl + C to stop streaming..." << std::endl;

	//tx_md.time_spec = usrp->get_time_now() + uhd::time_spec_t(0, SYM_CNT*start, 1e8/inter);
	
	tx_md.start_of_burst = true;
	tx_md.end_of_burst = false;
	tx_md.has_time_spec = false;

	usrp->get_device()->send(zeros, SYM_LEN, tx_md, C_FLOAT32, S_ONE_PKT);

	tx_md.start_of_burst = false;
	tx_md.end_of_burst = false; 

	// Send Signals until press ^C
	// HINT: You have to send signals here
	// How many symbols you have to send? Ans: sym_cnt
	// pkt: records the samples which we want to send
	
	// remove content of within while loop

	while (!stop_signal) {
		size_t sym_cnt = sample_cnt/SYM_LEN;
		size_t offset  = 0;
		for (size_t s = 0; s < sym_cnt; ++s) {
			// printf("%3d: %.4lf %.4lf\n", offset, pkt[offset].real(), pkt[offset].imag());
			// printf("%3d: %.4lf %.4lf\n", offset+1, pkt[offset+1].real(), pkt[offset+1].imag());
			offset += usrp->get_device()->send(pkt+offset, SYM_LEN, tx_md, C_FLOAT32, S_ONE_PKT);
		}
		// Clean the buffer of USRP
		for (size_t j = 0; j < 20; ++j)
			usrp->get_device()->send(zeros, SYM_LEN, tx_md, C_FLOAT32, S_ONE_PKT);
	}

	tx_md.start_of_burst    = false;
	tx_md.end_of_burst		= true; 
	
	usrp->get_device()->send(zeros, SYM_LEN, tx_md, C_FLOAT32, S_ONE_PKT);

    boost::this_thread::sleep(boost::posix_time::seconds(1));
	cout << "Terminate systems ... " << endl;
	return 0;
}
