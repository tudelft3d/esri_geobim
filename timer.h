#ifndef TIMER_H
#define TIMER_H

#include <map>
#include <string>
#include <chrono>

class timer {
	typedef std::map<std::string, double> timings_map;
	timings_map timings_;

public:
	// @todo something like a scoped_stopwatch using RAII
	class stopwatch {
		typedef std::chrono::time_point<std::chrono::steady_clock> timepoint;
		std::string key;
		timepoint start_, stop_;
		timings_map& timings;

	public:
		stopwatch(timings_map& timings) : timings(timings) {}
		stopwatch& start(const std::string& s) {
			key = s;
			start_ = std::chrono::steady_clock::now();
			return *this;
		}
		void stop() {
			stop_ = std::chrono::steady_clock::now();
			timings[key] += std::chrono::duration<double>(stop_ - start_).count();
		}
	};

	static timer& i() {
		static timer t;
		return t;
	}

	static stopwatch measure(const std::string& s) {
		return stopwatch(i().timings_).start(s);
	}

	static void print(std::ostream& s) {
		size_t max_key_length = 0;
		for (auto& p : i().timings_) {
			if (p.first.size() > max_key_length) {
				max_key_length = p.first.size();
			}
		}
		for (auto& p : i().timings_) {
			s << p.first << std::string(max_key_length - p.first.size(), ' ') << ":" << p.second << std::endl;
		}
	}
};

#endif