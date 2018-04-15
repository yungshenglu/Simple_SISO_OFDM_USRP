#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/usrp/multi_usrp.hpp>

#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/thread.hpp>
#include <iostream>
#include <complex>
#include <stdio.h>
#include <cmath>

#include "single_rx.h"

namespace po = boost::program_options;

//System parameters 
double freq, gain, thres;
double inter;
double rate;

// USRP
uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);

uhd::usrp::multi_usrp::sptr usrp;
string usrp_ip;

// Tx/Rx metadata
uhd::rx_metadata_t rx_md;
// uhd::tx_metadata_t tx_md;

// Buffer
gr_complex *pkt;

// File
string out_name; 

// Evaluation
size_t r_cnt;
double r_sec;
size_t s_cnt;

void init_usrp() {
	cout << "Initial USRP" << endl;
	
	usrp = uhd::usrp::multi_usrp::make(usrp_ip);
	usrp->set_rx_rate(rate);
	usrp->set_tx_rate(rate);

	usrp->set_rx_freq(freq);
	usrp->set_tx_freq(freq);

	usrp->set_rx_gain(gain);
	uhd::meta_range_t rx_range = usrp->get_rx_gain_range();
}

void sync_clock() {
	cout << "SYNC Clock" << endl;
	usrp->set_clock_config(uhd::clock_config_t::external());
	usrp->set_time_next_pps(uhd::time_spec_t(0.0));
}

void init_stream() {
    boost::this_thread::sleep(boost::posix_time::milliseconds(WARM_UP_TIME));
    stream_cmd.time_spec = time_start_recv  = uhd::time_spec_t(2.0) + usrp->get_time_now();
	cout << "Time to start receiving: " << time_start_recv.get_real_secs() << endl;
    stream_cmd.stream_now = false;
    usrp->issue_stream_cmd(stream_cmd);
}


void init_sys() {
	//Buffer initialize	
	printf("New receiving signal buffer\n");
	s_cnt = (size_t)(1e8/SYM_LEN * r_sec);
	pkt = new gr_complex[s_cnt];
	if (pkt!= NULL) {
		memset(pkt, 0, sizeof(gr_complex) * s_cnt);
	}

	rate = 1e8 / inter;
	freq = 1e9 * freq;

	init_usrp();
	//sync_clock();
	init_stream();
}

void dump_signals() {
	cout << endl << "Dump signals" << endl;
	FILE* out_file;
	char tmp_n[1000];
	sprintf(tmp_n, "%s", out_name.c_str());
	out_file = fopen(tmp_n, "wb");
	fwrite(pkt, sizeof(gr_complex), s_cnt, out_file);
	fclose(out_file);
}

void end_sys() {
	printf("Delete receiving signal buffer\n");
	if (pkt != NULL) {
		delete [] pkt;
	}
}

int UHD_SAFE_MAIN(int argc, char *argv[]) {
	size_t rx_cnt;
	rx_cnt = 0;

	uhd::set_thread_priority_safe();
	uhd::time_spec_t refer;

	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "help message")
		("r0", po::value<string>(&usrp_ip)->default_value("addr=192.168.20.2" ), "usrp's IP")
		("out", po::value<string>(&out_name)->default_value("wcs_trace/rx_signals.bin"), "signal file")
		("i", po::value< double >(&inter)->default_value(128), "interval of two sampling")
		("f", po::value< double >(&freq)->default_value(2.49), "RF center frequency in Hz")
		("g", po::value< double >(&gain)->default_value(30.0), "gain for the RF chain")
		("s", po::value< double >(&r_sec)->default_value(RECV_SEC), "recording seconds")
		("c", po::value< size_t >(&r_cnt)->default_value(90), "round count");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);


	if (vm.count("help")) {
		cout << boost::format("UHD TX samples from file %s") % desc << endl;
		return ~0;
	}
	// Initial systems
	init_sys();

	size_t cleaning, done_cleaning;
	done_cleaning = 0;
	cleaning = 0;
	while (cleaning < ANT_CNT) {
		if (!done_cleaning) {
			usrp->get_device()->recv(pkt, SYM_LEN, rx_md, C_FLOAT32, R_ONE_PKT);
			if(rx_md.time_spec.get_real_secs() >= time_start_recv.get_real_secs()) {
				done_cleaning = 1;
				cleaning++;
			}
		}
		// cout << cleaning << "-" << done_cleaning << " Clean ant" << i << " buff:" << rx_md.time_spec.get_real_secs() << endl;
	}
	
	// remove content of within while loop
	cout << endl << "# of recv samples: " << s_cnt << endl;
	while (rx_cnt < s_cnt) {
		size_t  read_cnt = 80;
		//At last recv(), modify read_cnt to receive the remaining samples
		if (s_cnt - rx_cnt < read_cnt)
			read_cnt = s_cnt - rx_cnt;
		rx_cnt += usrp->get_device()->recv(pkt+rx_cnt, read_cnt, rx_md, C_FLOAT32, R_ONE_PKT);

		if (rx_cnt < 100)
			cout << "Ant" << " recving at " << rx_md.time_spec.get_real_secs() << endl;
	}
	

	stream_cmd.stream_mode = uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS;
	usrp->issue_stream_cmd(stream_cmd);
	
	boost::this_thread::sleep(boost::posix_time::seconds(1)); //allow for some setup time
	
	
	// End systems
	dump_signals();
	end_sys();	
	
	cout << "Terminate systems ... " << endl << endl;
	return 0;
}

