#include <stdio.h>
#include <complex>
#include <cmath>
#include <fstream>

#include <stdint.h>
#include <gnuradio/gr_complex.h>

using namespace std;

#ifndef FLEXNEMO
#define FLEXNEMO

// USRP
#define C_FLOAT32 	uhd::io_type_t::COMPLEX_FLOAT32
#define R_ONE_PKT	uhd::device::RECV_MODE_ONE_PACKET
#define S_ONE_PKT	uhd::device::SEND_MODE_ONE_PACKET
#define STOP 		uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS
#define START 		uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS

// System setting 
const size_t WARM_UP_TIME = 1000;	//ms
const size_t NOTIFY_CNT = 5;	    //sample
const size_t WAIT_CNT = 10;

// USRP
const size_t ANT_CNT = 1;

size_t sample_cnt = 0;
#endif

#ifndef FFT
#define FFT
const double PI = 4 * atan(1);

const size_t SC_LEN = 64;
const size_t CP_LEN	= 16;
const size_t SYM_LEN = SC_LEN + CP_LEN; 
const size_t NUM_SYMS = 50;
const size_t STS_LEN = 480;
const size_t LTS_LEN = 160;
const size_t MAX_PKT_LEN = STS_LEN + LTS_LEN + LTS_LEN + NUM_SYMS * SYM_LEN;
const size_t ELEM_SIZE = sizeof(complex< float >);
const size_t SYM_SIZE = ELEM_SIZE * SYM_LEN;
#endif
