#ifndef MERSENNETWISTER
#define MERSENNETWISTER
#include <ctime>
//Note: Not standard.
#include <stdint.h>

unsigned long initseed() {
	unsigned long c = time(0);
	bool d = c % 2;
	c = c >> 2;
	c = c << 1;
	c += d;
	return c << 16 + c;
}

class MersenneTwister {
	unsigned long vector[624];
	unsigned short state;
public:
	MersenneTwister() {
		state = 0xFFFF;
		vector[0] = initseed();
		for (unsigned short a = 1; a<624; a++)
			vector[a] = (0x10DCD * vector[a - 1]) & 0xFFFFFFFF;
	}
	MersenneTwister(unsigned long seed) {
		state = 0xFFFF;
		vector[0] = seed;
		for (unsigned short a = 1; a<624; a++)
			vector[a] = (0x10DCD * vector[a - 1]) & 0xFFFFFFFF;
	}
	void gen() {
		unsigned long y = 0;
		for (unsigned short a = 0; a <= 622; a++) {
			y = this->vector[a] & 0x7FFFFFFF;
			y += unsigned long(double(this->vector[a + 1] & 0xFFFFFFFF) / 0xFFFFFFFF);
			if (!(y % 2))
				this->vector[a] = this->vector[(a + 367) % 624] ^ (y >> 1);
			else
				this->vector[a] = this->vector[(a + 397) % 624] ^ (y >> 1) ^ 0x9908B0DF;
		}
		y = this->vector[623] & 0x7FFFFFFF;
		y += (unsigned long)(double(this->vector[0] & 0xFFFFFFFF) / 0xFFFFFFFF);
		if (!(y % 2))
			this->vector[623] = this->vector[396] ^ (y >> 1);
		else
			this->vector[623] = this->vector[396] ^ (y >> 1) ^ 0x9908B0DF;
	}
	//0<=uint32_t<256^4-1
	unsigned long rnd() {
		unsigned long y;
		if (this->state == 0xFFFF) {
			this->state = 1;
			this->gen();
			y = this->vector[0];
		}
		else {
			y = this->vector[this->state];
			this->state++;
			y = y ^ (y >> 11);
			y = y ^ (y << 7) & 0x9D2C5680;
			y = y ^ (y << 15) & 0xEFC60000;
			y = y ^ (y >> 18);
			if (this->state == 623) {
				this->state = 0;
				this->gen();
			}
		}
		return y - 1;
	}
	//0<=double<1
	double rnd2() {
		unsigned long y = this->rnd();
		return y == 0xFFFFFFFF ? double(y - 1) / 0xFFFFFFFF : double(y) / 0xFFFFFFFF;
	}
};
#endif