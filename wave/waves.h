#ifndef WAVES_H
#define WAVES_H

void time_step(int n, int b,
		const double* u0,
		const double* u1,
		double* u2,
		double c, double dx, double dt);

#endif /* WAVES_H */
