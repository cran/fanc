/*  profiling.h
 *  - 2013/03/01
 *  - Haruhisa Nagata 
 */

#ifndef PROFILING_H

/* ------------------------------------------------ *
 *  Settings                                        *
 * ------------------------------------------------ */

// Enable time measurement: comment out this to be disable
#define PROFILING


/* ------------------------------------------------ *
 *  Code                                            *
 * ------------------------------------------------ */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
//#include <time.h>
#include <sys/time.h>
#include <R.h>

struct profiling {
	const char *name;
	int64_t time;
	uint64_t start;
} *g_prof = NULL;
int g_prof_n = 0;

uint64_t get_nanotime() {
//	struct timespec ts;
//	return clock_gettime(CLOCK_MONOTONIC, &ts) == 0 ?
//		(uint64_t)ts.tv_sec * 1000000000 + ts.tv_nsec : (uint64_t)-1;
	struct timeval t;
	return gettimeofday(&t, NULL) == 0 ?
		(uint64_t)t.tv_sec * 1000000000 + t.tv_usec * 1000 : (uint64_t)-1;
}

void prof_init(int n) {
#ifdef PROFILING
	if (g_prof == NULL) {
		g_prof = (struct profiling *)calloc(n, sizeof(struct profiling));
		if (g_prof == NULL) return;
	}
	g_prof_n = n;
	for (int i = 0; i < n; i++) {
		g_prof[i].name = NULL;
	}
#endif
}

void prof_start(int index, const char *name) {
#ifdef PROFILING
	if (g_prof == NULL || index >= g_prof_n) return;
	g_prof[index].name = name;
	g_prof[index].start = get_nanotime();
#endif
}

void prof_starti(int index) {
#ifdef PROFILING
	prof_start(index, "$INDEX");
#endif
}

void prof_stop(int index) {
#ifdef PROFILING
	if (g_prof == NULL || index >= g_prof_n) return;
	uint64_t t = get_nanotime();
	g_prof[index].time += t - g_prof[index].start;
	g_prof[index].start = t;
#endif
}

void prof_print() {
#ifdef PROFILING
	if (g_prof == NULL) return;
	for (int i = 0; i < g_prof_n; i++) {
		if (g_prof[i].name == NULL) continue;
		double t = (double)(g_prof[i].time / 1000000) / 1000;
		if (strcmp(g_prof[i].name, "$INDEX") == 0) {
			Rprintf("%8.3fsec: %d\n", t, i);
		} else {
			Rprintf("%8.3fsec: %s\n", t, g_prof[i].name);
		}
	}
	free(g_prof);
	g_prof = NULL;
	g_prof_n = 0;
#endif
}

#endif
