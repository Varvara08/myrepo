/* 
WARNING: This file was autogenerated by

   ../../../src/core/gen_featureconfig.py on Sun Jul 24 17:21:23 2016

   Do not modify it or your changes will be overwritten!
   Modify features.def instead.
*/
#ifndef _FEATURECONFIG_HPP
#define _FEATURECONFIG_HPP

/* Handle implications */
// COLLISION_DETECTION implies GHOST_FLAG
#if defined(COLLISION_DETECTION) && !defined(GHOST_FLAG)
#define GHOST_FLAG
#endif

// COLLISION_DETECTION implies GHOSTS_HAVE_BONDS
#if defined(COLLISION_DETECTION) && !defined(GHOSTS_HAVE_BONDS)
#define GHOSTS_HAVE_BONDS
#endif

// ROTATIONAL_INERTIA implies ROTATION
#if defined(ROTATIONAL_INERTIA) && !defined(ROTATION)
#define ROTATION
#endif

// ROTATION_PER_PARTICLE implies ROTATION
#if defined(ROTATION_PER_PARTICLE) && !defined(ROTATION)
#define ROTATION
#endif

// MOL_CUT implies VIRTUAL_SITES_COM
#if defined(MOL_CUT) && !defined(VIRTUAL_SITES_COM)
#define VIRTUAL_SITES_COM
#endif

// INTER_RF implies ELECTROSTATICS
#if defined(INTER_RF) && !defined(ELECTROSTATICS)
#define ELECTROSTATICS
#endif

// MMM1D_GPU implies PARTIAL_PERIODIC
#if defined(MMM1D_GPU) && !defined(PARTIAL_PERIODIC)
#define PARTIAL_PERIODIC
#endif

// MMM1D_GPU implies ELECTROSTATICS
#if defined(MMM1D_GPU) && !defined(ELECTROSTATICS)
#define ELECTROSTATICS
#endif

// VIRTUAL_SITES_RELATIVE implies ROTATION
#if defined(VIRTUAL_SITES_RELATIVE) && !defined(ROTATION)
#define ROTATION
#endif

// THERMOSTAT_IGNORE_NON_VIRTUAL implies VIRTUAL_SITES_THERMOSTAT
#if defined(THERMOSTAT_IGNORE_NON_VIRTUAL) && !defined(VIRTUAL_SITES_THERMOSTAT)
#define VIRTUAL_SITES_THERMOSTAT
#endif

// TRANS_DPD implies DPD
#if defined(TRANS_DPD) && !defined(DPD)
#define DPD
#endif

// DPD_MASS implies MASS
#if defined(DPD_MASS) && !defined(MASS)
#define MASS
#endif

// DPD_MASS implies DPD
#if defined(DPD_MASS) && !defined(DPD)
#define DPD
#endif

// TUNABLE_SLIP implies DPD
#if defined(TUNABLE_SLIP) && !defined(DPD)
#define DPD
#endif

// LB_BOUNDARIES implies LB
#if defined(LB_BOUNDARIES) && !defined(LB)
#define LB
#endif

// LB_BOUNDARIES implies CONSTRAINTS
#if defined(LB_BOUNDARIES) && !defined(CONSTRAINTS)
#define CONSTRAINTS
#endif

// LB_BOUNDARIES_GPU implies LB_GPU
#if defined(LB_BOUNDARIES_GPU) && !defined(LB_GPU)
#define LB_GPU
#endif

// LB_BOUNDARIES_GPU implies CONSTRAINTS
#if defined(LB_BOUNDARIES_GPU) && !defined(CONSTRAINTS)
#define CONSTRAINTS
#endif

// LB_ELECTROHYDRODYNAMICS implies LB
#if defined(LB_ELECTROHYDRODYNAMICS) && !defined(LB)
#define LB
#endif

// ELECTROKINETICS implies LB_GPU
#if defined(ELECTROKINETICS) && !defined(LB_GPU)
#define LB_GPU
#endif

// ELECTROKINETICS implies EXTERNAL_FORCES
#if defined(ELECTROKINETICS) && !defined(EXTERNAL_FORCES)
#define EXTERNAL_FORCES
#endif

// ELECTROKINETICS implies ELECTROSTATICS
#if defined(ELECTROKINETICS) && !defined(ELECTROSTATICS)
#define ELECTROSTATICS
#endif

// EK_BOUNDARIES implies ELECTROKINETICS
#if defined(EK_BOUNDARIES) && !defined(ELECTROKINETICS)
#define ELECTROKINETICS
#endif

// EK_BOUNDARIES implies LB_GPU
#if defined(EK_BOUNDARIES) && !defined(LB_GPU)
#define LB_GPU
#endif

// EK_BOUNDARIES implies LB_BOUNDARIES_GPU
#if defined(EK_BOUNDARIES) && !defined(LB_BOUNDARIES_GPU)
#define LB_BOUNDARIES_GPU
#endif

// EK_BOUNDARIES implies CONSTRAINTS
#if defined(EK_BOUNDARIES) && !defined(CONSTRAINTS)
#define CONSTRAINTS
#endif

// EK_BOUNDARIES implies EXTERNAL_FORCES
#if defined(EK_BOUNDARIES) && !defined(EXTERNAL_FORCES)
#define EXTERNAL_FORCES
#endif

// EK_BOUNDARIES implies ELECTROSTATICS
#if defined(EK_BOUNDARIES) && !defined(ELECTROSTATICS)
#define ELECTROSTATICS
#endif

// EK_REACTION implies ELECTROKINETICS
#if defined(EK_REACTION) && !defined(ELECTROKINETICS)
#define ELECTROKINETICS
#endif

// EK_REACTION implies LB_GPU
#if defined(EK_REACTION) && !defined(LB_GPU)
#define LB_GPU
#endif

// EK_REACTION implies EXTERNAL_FORCES
#if defined(EK_REACTION) && !defined(EXTERNAL_FORCES)
#define EXTERNAL_FORCES
#endif

// EK_REACTION implies ELECTROSTATICS
#if defined(EK_REACTION) && !defined(ELECTROSTATICS)
#define ELECTROSTATICS
#endif

// SHANCHEN implies LB_GPU
#if defined(SHANCHEN) && !defined(LB_GPU)
#define LB_GPU
#endif

// LENNARD_JONES_GENERIC implies LENNARD_JONES
#if defined(LENNARD_JONES_GENERIC) && !defined(LENNARD_JONES)
#define LENNARD_JONES
#endif

// GAY_BERNE implies ROTATION
#if defined(GAY_BERNE) && !defined(ROTATION)
#define ROTATION
#endif

// BOND_ANGLEDIST_HARMONIC implies BOND_ANGLEDIST
#if defined(BOND_ANGLEDIST_HARMONIC) && !defined(BOND_ANGLEDIST)
#define BOND_ANGLEDIST
#endif

// BOND_ANGLEDIST_HARMONIC implies CONSTRAINTS
#if defined(BOND_ANGLEDIST_HARMONIC) && !defined(CONSTRAINTS)
#define CONSTRAINTS
#endif

// BOND_ENDANGLEDIST_HARMONIC implies BOND_ENDANGLEDIST
#if defined(BOND_ENDANGLEDIST_HARMONIC) && !defined(BOND_ENDANGLEDIST)
#define BOND_ENDANGLEDIST
#endif

// BOND_ENDANGLEDIST_HARMONIC implies CONSTRAINTS
#if defined(BOND_ENDANGLEDIST_HARMONIC) && !defined(CONSTRAINTS)
#define CONSTRAINTS
#endif
/* Warn when derived switches are specified manually */
// P3M equals ELECTROSTATICS and FFTW
#ifdef P3M
#warning P3M is a derived switch and should not be set manually!
#elif defined(ELECTROSTATICS)  &&  defined(FFTW)
#define P3M
#endif

// DP3M equals DIPOLES and FFTW
#ifdef DP3M
#warning DP3M is a derived switch and should not be set manually!
#elif defined(DIPOLES)  &&  defined(FFTW)
#define DP3M
#endif

// VIRTUAL_SITES equals VIRTUAL_SITES_COM or VIRTUAL_SITES_RELATIVE
#ifdef VIRTUAL_SITES
#warning VIRTUAL_SITES is a derived switch and should not be set manually!
#elif defined(VIRTUAL_SITES_COM)  ||  defined(VIRTUAL_SITES_RELATIVE)
#define VIRTUAL_SITES
#endif

// DPD_MASS equals DPD_MASS_RED or DPD_MASS_LIN
#ifdef DPD_MASS
#warning DPD_MASS is a derived switch and should not be set manually!
#elif defined(DPD_MASS_RED)  ||  defined(DPD_MASS_LIN)
#define DPD_MASS
#endif

// LATTICE equals LB or LB_GPU
#ifdef LATTICE
#warning LATTICE is a derived switch and should not be set manually!
#elif defined(LB)  ||  defined(LB_GPU)
#define LATTICE
#endif

// USE_TEMPORARY equals LB or LB_GPU
#ifdef USE_TEMPORARY
#warning USE_TEMPORARY is a derived switch and should not be set manually!
#elif defined(LB)  ||  defined(LB_GPU)
#define USE_TEMPORARY
#endif

// BOND_ANGLE_OLD equals BOND_ANGLE_HARMONIC or BOND_ANGLE_COSINE or BOND_ANGLE_COSSQUARE
#ifdef BOND_ANGLE_OLD
#warning BOND_ANGLE_OLD is a derived switch and should not be set manually!
#elif defined(BOND_ANGLE_HARMONIC)  ||  defined(BOND_ANGLE_COSINE)  ||  defined(BOND_ANGLE_COSSQUARE)
#define BOND_ANGLE_OLD
#endif

extern const char* FEATURES[];
extern const int NUM_FEATURES;

#endif /* of _FEATURECONFIG_HPP */