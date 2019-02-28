// π itself
pub const pi : f64 = 3.1415926535897932384626433832795029;
// Boltzmann constant in Eh/K
pub const kB : f64 = 3.166808578545117e-06;
// speed of light c in vacuum in a.u.
pub const lightspeed : f64 = 137.0359990740;
// convert bohr (a.u.) to Ångström and back
pub const autoaa : f64 = 0.52917726;
pub const aatoau : f64 = 1.0/autoaa;
// convert Hartree to eV and back
pub const autoev : f64 = 27.21138505;
pub const evtoau : f64 = 1.0/autoev;
// convert Hartree to kcal/mol and back
pub const autokcal : f64 = 627.50947428;
pub const kcaltoau : f64 = 1.0/autokcal;
// convert Hartree to kJ/mol and back
pub const autokj : f64 = 2625.49964038;
pub const kjtoau : f64 = 1.0/autokj;
// convert Hartree to reciproce centimeters/wavenumbers and back
pub const autorcm : f64 = 219474.63067;
pub const autowav : f64 = autorcm;
pub const rcmtoau : f64 = 1.0/autorcm;
pub const wavtoau : f64 = 1.0/autowav;
// convert Hartree to nanometers and back
pub const autonm : f64 = 45.56335266;
pub const nmtoau : f64 = 1.0/autonm;
// masses
// amu -> kg :: conversion from atomic mass units to kg
// me  -> kg :: electron mass (a.u.) in kg
// amu -> au :: conversion from a.u. to amu
pub const amutokg : f64 = 1.660539040e-27;
pub const kgtoamu : f64 = 1.0/amutokg;
pub const metokg  : f64 = 9.10938356e-31;
pub const kgtome  : f64 = 1.0/metokg;
pub const amutoau : f64 = amutokg*kgtome;
pub const autoamu : f64 = kgtoamu*metokg;
// femtosectons to atomic time units
pub const fstoau : f64 = 41.3413733365614;
pub const autofs : f64 = 1.0/fstoau;
// Coulomb to atomic charge units (electrons)
pub const autoc : f64 = 1.6021766208e-19;
pub const ctoau : f64 = 1.0/autoc;
// Debye to atomic units
pub const autod : f64 = autoc * lightspeed * autoaa * autoaa * fstoau * 1.0e+16;
pub const dtoau : f64 = 1.0/autod;
