/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright (C) 2026 Syota Sasaki
 *
 * Time-symmetric Brill wave initial data from a conformal factor table.
 * See the README for the table layout and the metric.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <hdf5.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define FORMAT_NAME "brill-psi-table"
#define FORMAT_VERSION 1
#define MAX_POINTS 12

/* The table, read once per process at PARAMCHECK. */
typedef struct {
  int loaded;
  int n_s, n_theta;
  double h_s, h_theta;
  double scale;           /* L in r = L tan(pi s / 2) */
  double parity_s_upper;  /* parity of u about s = 1 */
  double amplitude, sigma;
  int exponent_factor;
  double *u;              /* n_s * n_theta, C order (i_s, i_theta) */
} psi_table;

static psi_table table = {0};

/* ------------------------------------------------------------------ */
/* HDF5 helpers                                                        */
/* ------------------------------------------------------------------ */

static void check_h5(herr_t status, const char *what) {
  if (status < 0) {
    CCTK_VERROR("HDF5 error while reading %s from the conformal factor table", what);
  }
}

static double read_double_attr(hid_t file, const char *name) {
  double value;
  hid_t attr = H5Aopen(file, name, H5P_DEFAULT);
  if (attr < 0) {
    CCTK_VERROR("table attribute \"%s\" is missing", name);
  }
  check_h5(H5Aread(attr, H5T_NATIVE_DOUBLE, &value), name);
  H5Aclose(attr);
  return value;
}

static int read_int_attr(hid_t file, const char *name) {
  int value;
  hid_t attr = H5Aopen(file, name, H5P_DEFAULT);
  if (attr < 0) {
    CCTK_VERROR("table attribute \"%s\" is missing", name);
  }
  check_h5(H5Aread(attr, H5T_NATIVE_INT, &value), name);
  H5Aclose(attr);
  return value;
}

/* Fixed-length string attribute; the writer stores ASCII bytes. */
static void read_string_attr(hid_t file, const char *name, char *buffer, size_t size) {
  hid_t attr = H5Aopen(file, name, H5P_DEFAULT);
  if (attr < 0) {
    CCTK_VERROR("table attribute \"%s\" is missing", name);
  }
  hid_t type = H5Aget_type(attr);
  if (H5Tget_class(type) != H5T_STRING || H5Tis_variable_str(type)) {
    CCTK_VERROR("table attribute \"%s\" is not a fixed-length string", name);
  }
  size_t length = H5Tget_size(type);
  if (length >= size) {
    CCTK_VERROR("table attribute \"%s\" is too long", name);
  }
  /* One byte more than the stored length: the stored string is padded,
     not terminated, and HDF5 writes the terminating NUL into the last
     byte of the memory type. */
  hid_t memtype = H5Tcopy(H5T_C_S1);
  H5Tset_size(memtype, length + 1);
  memset(buffer, 0, size);
  check_h5(H5Aread(attr, memtype, buffer), name);
  H5Tclose(memtype);
  H5Tclose(type);
  H5Aclose(attr);
}

/* ------------------------------------------------------------------ */
/* Loading and validation                                              */
/* ------------------------------------------------------------------ */

static int matches(double expected, double actual, double tolerance) {
  return fabs(expected - actual) <= tolerance * fmax(1.0, fabs(expected));
}

void BrillPsiTable_Load(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  char text[256];

  if (table.loaded) {
    return;
  }
  if (strlen(table_file) == 0) {
    CCTK_PARAMWARN("BrillPsiTable::table_file is empty");
    return;
  }

  hid_t file = H5Fopen(table_file, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) {
    CCTK_VPARAMWARN("cannot open table \"%s\"", table_file);
    return;
  }

  read_string_attr(file, "format", text, sizeof text);
  if (strcmp(text, FORMAT_NAME) != 0) {
    CCTK_VERROR("\"%s\" is not a %s file (format \"%s\")", table_file, FORMAT_NAME, text);
  }
  if (read_int_attr(file, "format_version") != FORMAT_VERSION) {
    CCTK_VERROR("\"%s\": unsupported format_version", table_file);
  }
  read_string_attr(file, "family", text, sizeof text);
  if (strcmp(text, "holz") != 0) {
    CCTK_VERROR("\"%s\": unsupported seed family \"%s\"", table_file, text);
  }

  table.amplitude = read_double_attr(file, "amplitude");
  table.sigma = read_double_attr(file, "sigma");
  table.exponent_factor = read_int_attr(file, "exponent_factor");
  table.scale = read_double_attr(file, "scale_L");
  table.n_s = read_int_attr(file, "n_s");
  table.n_theta = read_int_attr(file, "n_theta");
  const double adm_mass = read_double_attr(file, "adm_mass_flux");

  if (!matches(amplitude, table.amplitude, match_tolerance)) {
    CCTK_VPARAMWARN("amplitude %.17g does not match the table (%.17g)", (double)amplitude,
                    table.amplitude);
  }
  if (!matches(sigma, table.sigma, match_tolerance)) {
    CCTK_VPARAMWARN("sigma %.17g does not match the table (%.17g)", (double)sigma, table.sigma);
  }
  if (exponent_factor != table.exponent_factor) {
    CCTK_VPARAMWARN("exponent_factor %d does not match the table (%d)", (int)exponent_factor,
                    table.exponent_factor);
  }
  if (read_int_attr(file, "parity_s_lower") != 1 || read_int_attr(file, "parity_theta_lower") != 1 ||
      read_int_attr(file, "parity_theta_upper") != 1) {
    CCTK_VERROR("\"%s\": unexpected table parities", table_file);
  }
  table.parity_s_upper = read_int_attr(file, "parity_s_upper");
  if (table.n_s < interpolation_points || table.n_theta < interpolation_points) {
    CCTK_VERROR("\"%s\": table (%d x %d) smaller than the interpolation stencil (%d)", table_file,
                table.n_s, table.n_theta, (int)interpolation_points);
  }

  hid_t dataset = H5Dopen2(file, "u", H5P_DEFAULT);
  if (dataset < 0) {
    CCTK_VERROR("\"%s\": dataset \"u\" is missing", table_file);
  }
  hid_t space = H5Dget_space(dataset);
  hsize_t dims[2];
  if (H5Sget_simple_extent_ndims(space) != 2) {
    CCTK_VERROR("\"%s\": dataset \"u\" is not two-dimensional", table_file);
  }
  H5Sget_simple_extent_dims(space, dims, NULL);
  if ((int)dims[0] != table.n_s || (int)dims[1] != table.n_theta) {
    CCTK_VERROR("\"%s\": dataset shape does not match n_s, n_theta", table_file);
  }
  table.u = malloc(sizeof(double) * table.n_s * table.n_theta);
  if (table.u == NULL) {
    CCTK_ERROR("out of memory for the conformal factor table");
  }
  check_h5(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, table.u), "u");
  H5Sclose(space);
  H5Dclose(dataset);
  H5Fclose(file);

  table.h_s = 1.0 / table.n_s;
  table.h_theta = 0.5 * M_PI / table.n_theta;
  table.loaded = 1;
  CCTK_VINFO("loaded \"%s\": A = %.17g, sigma = %g, k = %d, %d x %d cells, L = %g, "
             "M_ADM = %.12f", table_file, table.amplitude, table.sigma, table.exponent_factor,
             table.n_s, table.n_theta, table.scale, adm_mass);
}

/* ------------------------------------------------------------------ */
/* Interpolation                                                       */
/* ------------------------------------------------------------------ */

/* Lagrange weights at x for the integer nodes first, ..., first + n - 1. */
static void lagrange_weights(double x, int first, int n, double *w) {
  for (int k = 0; k < n; ++k) {
    double weight = 1.0;
    for (int m = 0; m < n; ++m) {
      if (m != k) {
        weight *= (x - (first + m)) / (double)(k - m);
      }
    }
    w[k] = weight;
  }
}

/* Cell-centred mirror: even at the lower face, `upper` at the upper face. */
static int mirror(int index, int n, double upper, double *parity) {
  *parity = 1.0;
  if (index < 0) {
    return -index - 1;
  }
  if (index >= n) {
    *parity = upper;
    return 2 * n - 1 - index;
  }
  return index;
}

/* u = psi - 1 at radius r and polar angle theta in [0, pi/2]. */
static double table_u(double r, double theta, int n_points) {
  const double s = 2.0 / M_PI * atan(r / table.scale);
  const double xs = s / table.h_s - 0.5;
  const double xt = theta / table.h_theta - 0.5;
  const int half = n_points / 2;
  const int fs = (int)floor(xs) - half + 1;
  const int ft = (int)floor(xt) - half + 1;
  double ws[MAX_POINTS], wt[MAX_POINTS];
  lagrange_weights(xs, fs, n_points, ws);
  lagrange_weights(xt, ft, n_points, wt);

  double value = 0.0;
  for (int a = 0; a < n_points; ++a) {
    double ps;
    const int i = mirror(fs + a, table.n_s, table.parity_s_upper, &ps);
    double row = 0.0;
    for (int b = 0; b < n_points; ++b) {
      double pt;
      const int j = mirror(ft + b, table.n_theta, 1.0, &pt);
      row += wt[b] * pt * table.u[i * table.n_theta + j];
    }
    value += ws[a] * ps * row;
  }
  return value;
}

/* ------------------------------------------------------------------ */
/* Initial data                                                        */
/* ------------------------------------------------------------------ */

void BrillPsiTable_InitialData(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (!table.loaded) {
    CCTK_ERROR("conformal factor table not loaded");
  }
  const double k = table.exponent_factor;
  const double inv_sigma2 = 1.0 / (table.sigma * table.sigma);
  const int n_points = interpolation_points;
  const int np = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

#pragma omp parallel for
  for (int idx = 0; idx < np; ++idx) {
    const double xx = x[idx], yy = y[idx], zz = z[idx];
    const double rho2 = xx * xx + yy * yy;
    const double r2 = rho2 + zz * zz;
    const double rr = sqrt(r2);
    const double theta = atan2(sqrt(rho2), fabs(zz));

    const double p = 1.0 + table_u(rr, theta, n_points);
    const double p4 = p * p * p * p;

    /* k q / rho^2 and k q, exact on the axis. */
    const double kq_over_rho2 = k * table.amplitude * inv_sigma2 * exp(-r2 * inv_sigma2);
    const double kq = kq_over_rho2 * rho2;
    const double f = kq == 0.0 ? kq_over_rho2 : kq_over_rho2 * expm1(kq) / kq;

    gxx[idx] = p4 * (1.0 + f * xx * xx);
    gxy[idx] = p4 * f * xx * yy;
    gxz[idx] = 0.0;
    gyy[idx] = p4 * (1.0 + f * yy * yy);
    gyz[idx] = 0.0;
    gzz[idx] = p4 * exp(kq);

    kxx[idx] = 0.0;
    kxy[idx] = 0.0;
    kxz[idx] = 0.0;
    kyy[idx] = 0.0;
    kyz[idx] = 0.0;
    kzz[idx] = 0.0;

    psi[idx] = p;
  }
}

void BrillPsiTable_Free(CCTK_ARGUMENTS) {
  free(table.u);
  table.u = NULL;
  table.loaded = 0;
}
