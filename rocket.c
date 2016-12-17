#include <math.h>
#include <stdio.h>

/* SIMULATION PARAMETERS */
const double t0 = 0;		/* Initial time */
const double t1 = 8.849;	/* Time when thrust ceases */
const double tf = 150;		/* Time when simulation ends */
const double dt = 1e-2;		/* Timestep */

/* LAUNCH PARAMETERS */
const double h0 = 0;		/* Initial altitude */

/* LOCAL PHYSICS PARAMETERS */
const double g = -9.81;		/* Local gravity */
const double p0 = 101325;	/* Pressure at sea level */
const double k = -1.3024E-4;	/* Exponential constant factor for pressure calculation */
const double T0 = 300;		/* Temperature at sea level */

/* ROCKET PARAMETERS */
const double m0 = 2530;		/* Initial mass */
const double mf = 2530 - 1050;	/* Final mass */
const double d = 0.175;		/* Coefficient of drag */
const double A = 3.14159 / 4;	/* Cross-sectional area */
double dm_dt = -1050 / t1;      /* Mass change over time (when rockets firing) */
double dT_dI = dm_dt * g;	/* Thrust vs specific impulse */
const double Iv = 165;		/* Specific impulse in a vacuum */
const double Ia = 140;		/* Specific impulse at sea level */
double dI_dp = (Ia - Iv) / p0;  /* Specific impulse versus pressure */

int main(void) {
  double t;
  double h, P, Isp, T, Sf, m, a, dh, v = 0;
  double Fd;
  h = h0;

  int pct, old_pct = -1;
  printf("t (s),h (m),P (Pa),Isp (sec),Σf (N),m (kg),a (m/s^2),Δh (m),Vf (m/s)\n");
  fprintf(stderr, "Simulating...\n");
  for(t = t0; t < tf; t += dt) {
    P = p0 * exp(k * h);
    Isp = dI_dp * P + Iv;
    if(t < t1) {
      T = Isp * dT_dI;
      m = dm_dt * t + m0;
    } else {
      T = 0;
    }
    Fd = 0.5 * (P / (287.053 * T0)) * v * v * d * A;
    Sf = T + m * g - ((v > 0) ? Fd : -Fd);
    a = Sf / m;
    v += a * dt;
    dh = a * dt * dt / 2 + dt * v;
    printf("%.4lf,%.1lf,%.0lf,%.2lf,%.0lf,%.0lf,%.4lf,%.2lf,%g\n",
	   t, h, P, Isp, Sf, m, a, dh, v);
    h += dh;
    pct = (int) t * 100 / tf;
    if(pct > old_pct) {
      int i;
      for(i = 0; i < pct; i++) putc('#', stderr);
      fprintf(stderr, " %d%%\n\e[A", pct);
      old_pct = pct;
    }
    if(h < -10) break;
  }
  fprintf(stderr, "\n");
}
