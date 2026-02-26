#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Adopted from PostCactus (author: Wolfgang Kastaun)
Standalone unit converter for SI / CGS / Cactus geometric units (CU).

- Includes:
  * Units class and derived units
  * Physical constants (SI)
  * Geometric-unit constructors
  * Predefined unit systems: SI, CGS, CACTUS (CU)
  * PostCactus-style conversion factor report writer
  * CLI for converting individual values

Notes
-----
- "CU" (Cactus Units) are geometric units with mass unit = 1 solar mass.
- For magnetic field conversion:
    * CGS magnetic field is interpreted as Gauss
    * SI magnetic field is interpreted as Tesla
  The conversion constant is retained from your original script for compatibility.

Examples
--------
1) Write the conversion report:
   python unit_converter.py --write-report

2) Convert 10 km to CU:
   python unit_converter.py convert --value 10 --quantity length_km --from-unit si --to-unit cu

3) Convert density from CU to CGS:
   python unit_converter.py convert --value 1.0 --quantity density --from-unit cu --to-unit cgs

4) Convert magnetic field from CU to Gauss (CGS):
   python unit_converter.py convert --value 1.0 --quantity magnetic_field --from-unit cu --to-unit cgs
"""

import math
import argparse
from typing import Dict


# =============================================================================
# Core unit system
# =============================================================================

class Units:
    """
    Class representing unit conversion.

    A unit system is specified by base units (length, time, mass).
    Derived units are computed automatically.

    The division operator returns a Units object that expresses "self" in terms
    of "base" (e.g., CU/CGS gives factors to convert CU -> CGS).
    """

    def __init__(self, ulength, utime, umass):
        self.length   = float(ulength)
        self.time     = float(utime)
        self.mass     = float(umass)

        self.freq     = 1.0 / self.time
        self.velocity = self.length / self.time
        self.accel    = self.velocity / self.time

        self.force    = self.accel * self.mass
        self.area     = self.length ** 2
        self.volume   = self.length ** 3
        self.density  = self.mass / self.volume
        self.pressure = self.force / self.area
        self.power    = self.force * self.velocity
        self.energy   = self.force * self.length
        self.edens    = self.energy / self.volume
        self.angmom   = self.energy * self.time
        self.minertia = self.mass * self.area

    def __truediv__(self, base):
        """Express this unit system in terms of another unit system."""
        return Units(
            self.length / base.length,
            self.time / base.time,
            self.mass / base.mass
        )

    # For legacy compatibility
    __div__ = __truediv__


# =============================================================================
# SI constants (same values as your original script for compatibility)
# =============================================================================

C_SI          = 299792458.0            # Vacuum speed of light [m/s]
G_SI          = 6.6743e-11             # Gravitational constant [m^3 kg^-1 s^-2]
M_SOL_SI      = 1.98892e30             # Solar mass [kg]
EV_SI         = 1.602176565e-19        # Electron volt [J]
MEV_SI        = 1e6 * EV_SI
UAMU_SI       = 931.494061 * MEV_SI / (C_SI ** 2)  # Unified atomic mass unit [kg]
KB_SI         = 1.3806488e-23          # Boltzmann constant [J/K]
M_ELECTRON_SI = 9.10938291e-31         # Electron mass [kg]
M_PROTON_SI   = 1.672621777e-27        # Proton mass [kg]
M_NEUTRON_SI  = 1.674927351e-27        # Neutron mass [kg]
LIGHTYEAR_SI  = 9460730472580800.0     # Lightyear [m]
PARSEC_SI     = 30.856776e15           # Parsec [m]
N_AVOGADRO    = 6.02214129e23          # [1/mol]
H_SI          = 6.62606957e-34         # Planck constant [J s]
HBAR_SI       = 1.054571726e-34        # Reduced Planck constant [J s]


# =============================================================================
# Geometric-unit constructors
# =============================================================================

def geom_ulength(ulength_si: float) -> Units:
    """
    Create a geometric unit system (in SI) where the geometric length unit
    is ulength_si [m].
    """
    return Units(
        ulength_si,
        ulength_si / C_SI,
        ulength_si * (C_SI ** 2) / G_SI
    )


def geom_udensity(udensity_si: float) -> Units:
    """
    Create a geometric unit system (in SI) where the density unit is given
    by udensity_si [kg/m^3].
    """
    return geom_ulength(C_SI / math.sqrt(G_SI * udensity_si))


def geom_umass(umass_si: float) -> Units:
    """
    Create a geometric unit system (in SI) where the mass unit is given
    by umass_si [kg].
    """
    return geom_ulength(umass_si * G_SI / (C_SI ** 2))


# =============================================================================
# Predefined unit systems
# =============================================================================

SI_UNITS      = Units(1.0, 1.0, 1.0)        # SI base units: m, s, kg
CGS_UNITS     = Units(1e-2, 1.0, 1e-3)      # CGS base units: cm, s, g
CACTUS_UNITS  = geom_umass(M_SOL_SI)        # Geometric units, 1 mass unit = 1 solar mass


# Aliases (kept for convenience / familiarity)
SI  = SI_UNITS
CGS = CGS_UNITS
CU  = CACTUS_UNITS


# =============================================================================
# Common abbreviations / handy factors (formerly unit_abbrev.py)
# =============================================================================

CU_CGS = CU / CGS
CGS_CU = CGS / CU

KM_CU      = 1e3 / CU.length
PARSEC_CU  = PARSEC_SI / CU.length
MS_CU      = 1e-3 / CU.time
HZ_CU      = 1.0 / CU.freq
KHZ_CU     = 1e3 / CU.freq

# Retained exactly from your original code for compatibility
# This is "1 Gauss in CU", so "1 CU in Gauss" = 1 / GAUSS_CU
GAUSS_CU = 1.0 / (8.351966458327302e+15 * 1e4)

TEX_KM        = r'\mathrm{km}'
TEX_MS        = r'\mathrm{ms}'
TEX_GAUSS     = r'\mathrm{G}'
TEX_MEV       = r'\mathrm{MeV}'
TEX_DENS_CGS  = r'\mathrm{g}/\mathrm{cm}^{3}'


# =============================================================================
# Conversion helpers
# =============================================================================

SYSTEMS: Dict[str, Units] = {
    "si": SI,
    "cgs": CGS,
    "cu": CU,
}

# Supported derived quantities directly from Units attributes
UNITS_QUANTITIES = {
    "length": "length",
    "time": "time",
    "mass": "mass",
    "freq": "freq",
    "frequency": "freq",
    "velocity": "velocity",
    "speed": "velocity",
    "accel": "accel",
    "acceleration": "accel",
    "force": "force",
    "area": "area",
    "volume": "volume",
    "density": "density",
    "pressure": "pressure",
    "power": "power",
    "energy": "energy",
    "edens": "edens",
    "energy_density": "edens",
    "angmom": "angmom",
    "angular_momentum": "angmom",
    "minertia": "minertia",
    "moment_of_inertia": "minertia",
}

# Convenience physical-unit aliases
# These convert first to SI base form, then into target system.
SPECIAL_QUANTITIES = {
    # length
    "length_m":   ("length", 1.0, "si"),
    "length_cm":  ("length", 1e-2, "si"),
    "length_km":  ("length", 1e3, "si"),
    "parsec":     ("length", PARSEC_SI, "si"),
    "lightyear":  ("length", LIGHTYEAR_SI, "si"),

    # time
    "time_s":     ("time", 1.0, "si"),
    "time_ms":    ("time", 1e-3, "si"),

    # mass
    "mass_kg":    ("mass", 1.0, "si"),
    "mass_g":     ("mass", 1e-3, "si"),
    "mass_msun":  ("mass", M_SOL_SI, "si"),

    # frequency
    "hz":         ("freq", 1.0, "si"),
    "khz":        ("freq", 1e3, "si"),
}


def _convert_units_quantity(value: float, quantity_attr: str, from_unit: str, to_unit: str) -> float:
    """Convert using Units-derived attributes."""
    fu = SYSTEMS[from_unit]
    tu = SYSTEMS[to_unit]
    factor = getattr(fu / tu, quantity_attr)   # quantity in "to_unit" per 1 "from_unit"
    return value * factor


def _convert_magnetic_field(value: float, from_unit: str, to_unit: str) -> float:
    """
    Convert magnetic field among:
      - CU
      - CGS (Gauss)
      - SI  (Tesla)
    """
    from_unit = from_unit.lower()
    to_unit = to_unit.lower()

    # Convert input to Gauss as an intermediate
    if from_unit == "cgs":
        b_gauss = value
    elif from_unit == "si":
        b_gauss = value * 1e4  # Tesla -> Gauss
    elif from_unit == "cu":
        b_gauss = value * (1.0 / GAUSS_CU)  # CU -> Gauss
    else:
        raise ValueError(f"Unsupported system for magnetic_field: {from_unit}")

    # Convert Gauss -> target
    if to_unit == "cgs":
        return b_gauss
    elif to_unit == "si":
        return b_gauss / 1e4   # Gauss -> Tesla
    elif to_unit == "cu":
        return b_gauss * GAUSS_CU
    else:
        raise ValueError(f"Unsupported system for magnetic_field: {to_unit}")


def convert_value(value: float, quantity: str, from_unit: str, to_unit: str) -> float:
    """
    Generic converter.

    quantity options include:
      - direct Units quantities: length, time, mass, density, pressure, ...
      - magnetic_field
      - convenience aliases: length_km, length_cm, time_ms, hz, khz, mass_msun, ...
    """
    q = quantity.lower()
    f = from_unit.lower()
    t = to_unit.lower()

    if f not in SYSTEMS or t not in SYSTEMS:
        raise ValueError(f"Unknown system. Allowed: {sorted(SYSTEMS.keys())}")

    if q in ("magnetic_field", "bfield", "b", "gauss"):
        return _convert_magnetic_field(value, f, t)

    if q in UNITS_QUANTITIES:
        return _convert_units_quantity(value, UNITS_QUANTITIES[q], f, t)

    if q in SPECIAL_QUANTITIES:
        # SPECIAL_QUANTITIES interpret input as a named physical unit encoded in SI scale.
        # Example: quantity='length_km', from_unit='si' means "value km" -> SI meters -> target system
        base_attr, si_scale, canonical_system = SPECIAL_QUANTITIES[q]
        if f != canonical_system:
            raise ValueError(
                f"Quantity '{quantity}' is defined as a convenience alias in '{canonical_system.upper()}' input form. "
                f"Use --from-unit {canonical_system}"
            )
        # Convert the raw physical unit into canonical system base (SI) quantity first:
        # e.g. 10 km -> 10 * 1000 (meters), then meters -> target
        value_in_canonical_base = value * si_scale
        return _convert_units_quantity(value_in_canonical_base, base_attr, canonical_system, t)

    raise ValueError(
        f"Unknown quantity '{quantity}'. "
        f"Supported direct quantities: {sorted(set(UNITS_QUANTITIES.keys()))}; "
        f"special: {sorted(SPECIAL_QUANTITIES.keys())}; magnetic_field"
    )


# =============================================================================
# PostCactus-style conversion factor report (standalone)
# =============================================================================

def write_conversion_report(filename: str = "conversion_factors_CU_CGS.txt") -> None:
    """
    Write conversion factors from Cactus Units (CU)
    """
    uc   = CU / CGS  # 1 CU in CGS
    ucsi = CU / SI   # 1 CU in SI

    with open(filename, "w", encoding="utf-8") as fn:
        fn.write("Conversion factors from Cactus Units (CU):- \n \n \n")

        fn.write("Length: \n")
        fn.write("1 CU = %10.16e cm = %10.16e km \n" % (uc.length, uc.length / 10**5))
        fn.write("1 km = 1e+05 cm = %10.16e CU \n \n" % (10**5 / uc.length))

        fn.write("Time: \n")
        fn.write("1 CU = %10.16e ms \n" % (uc.time / 10**-3))
        fn.write("1 ms = %10.16e CU \n \n" % (10**-3 / uc.time))

        fn.write("Mass: \n")
        fn.write("1 CU = %10.16e g \n" % (uc.mass))
        fn.write("1 g = %10.16e CU \n \n" % (1 / uc.mass))

        fn.write("Speed: \n")
        fn.write("1 CU = %10.16e cm/s \n" % (uc.velocity))
        fn.write("1 cm/s = %10.16e CU \n \n" % (1 / uc.velocity))

        fn.write("Rest mass density: \n")
        fn.write("1 CU = %10.16e g/cm^3 \n" % (uc.density))
        fn.write("1 g/cm^3 = %10.16e CU \n \n" % (1 / uc.density))

        fn.write("Magnetic field (B): \n")
        fn.write("1 CU = %10.16e Gauss \n" % (1 / GAUSS_CU))
        fn.write("1 Gauss  = %10.16e CU \n \n \n" % (GAUSS_CU))

        fn.write("Other useful conversions:- \n \n \n")

        fn.write("Frequency: \n")
        fn.write("1 CU = %10.16e Hz \n" % (1 / uc.time))
        fn.write("1 Hz = %10.16e CU \n \n" % (uc.time))

        fn.write("Acceleration: \n")
        fn.write("1 CU = %10.16e cm/s^2 \n" % (uc.accel))
        fn.write("1 cm/s^2 = %10.16e CU \n \n" % (1 / uc.accel))

        fn.write("Area: \n")
        fn.write("1 CU = %10.16e cm^2 \n" % (uc.area))
        fn.write("1 cm^2 = %10.16e CU \n \n" % (1 / uc.area))

        fn.write("Volume: \n")
        fn.write("1 CU = %10.16e cm^3 \n" % (uc.volume))
        fn.write("1 cm^3 = %10.16e CU \n \n" % (1 / uc.volume))

        fn.write("Force: \n")
        fn.write("1 CU = %10.16e g cm/s^2 \n" % (uc.force))
        fn.write("1 g cm/s^2 = %10.16e CU \n \n" % (1 / uc.force))

        fn.write("Pressure: \n")
        fn.write("1 CU = %10.16e g/(cm s^2) \n" % (uc.pressure))
        fn.write("1 g/(cm s^2) = %10.16e CU \n \n" % (1 / uc.pressure))

        fn.write("Energy: \n")
        fn.write("1 CU = %10.16e g cm^2/s^2  [erg] \n" % (uc.energy))
        fn.write("1 erg = %10.16e CU \n \n" % (1 / uc.energy))

        fn.write("Luminosity (Power): \n")
        fn.write("1 CU = %10.16e g cm^2/ s^3  [erg/s] \n" % (uc.power))
        fn.write("1 erg/s = %10.16e CU \n \n" % (1 / uc.power))

        fn.write("Energy density: \n")
        fn.write("1 CU = %10.16e erg/cm^3 \n" % (uc.edens))
        fn.write("1 erg/cm^3 = %10.16e CU \n \n" % (1 / uc.edens))

        fn.write("Angular momentum: \n")
        fn.write("1 CU = %10.16e g cm^2/s \n" % (uc.angmom))
        fn.write("1 g cm^2/s = %10.16e CU \n \n" % (1 / uc.angmom))

        fn.write("Moment of inertia: \n")
        fn.write("1 CU = %10.16e g cm^2 \n" % (uc.minertia))
        fn.write("1 g cm^2  = %10.16e CU \n \n" % (1 / uc.minertia))

        fn.write("k, polytropic proportionality constant (P = k rho^2): \n")
        fn.write("1 CU = %10.16e g^-1 cm^5 s^-2 \n" % (uc.pressure / (uc.density ** 2)))
        fn.write("1 g^-1 cm^5 s^-2  = %10.16e CU \n \n" % (uc.density ** 2 / uc.pressure))

        fn.write("For k in CU, poly_rho in SI: \n")
        fn.write("poly_rho = (%10.16e kg/m^3)/k  \n" % (ucsi.density))


# =============================================================================
# CLI
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Standalone SI/CGS/Cactus (CU) unit converter"
    )

    p.add_argument(
        "--write-report",
        action="store_true",
        help="Write conversion factor report to text file"
    )
    p.add_argument(
        "--report-file",
        default="conversion_factors_CU_CGS.txt",
        help="Output filename for --write-report (default: %(default)s)"
    )

    sub = p.add_subparsers(dest="command")

    c = sub.add_parser("convert", help="Convert a value between SI/CGS/CU")
    c.add_argument("--value", type=float, required=True, help="Numeric value to convert")
    c.add_argument(
        "--quantity", required=True,
        help=(
            "Quantity name, e.g. length, time, mass, density, pressure, velocity, energy, "
            "magnetic_field, length_km, length_cm, time_ms, hz, khz, mass_msun"
        )
    )
    c.add_argument(
        "--from-unit", required=True, choices=sorted(SYSTEMS.keys()),
        help="Input unit system: si, cgs, cu"
    )
    c.add_argument(
        "--to-unit", required=True, choices=sorted(SYSTEMS.keys()),
        help="Output unit system: si, cgs, cu"
    )

    sub.add_parser("list", help="List supported quantities")

    return p


def print_supported_quantities():
    print("\nDirect quantities (derived from Units):")
    for k in sorted(set(UNITS_QUANTITIES.keys())):
        print(f"  - {k}")

    print("\nSpecial convenience quantities:")
    for k in sorted(SPECIAL_QUANTITIES.keys()):
        print(f"  - {k}")

    print("\nMagnetic field:")
    print("  - magnetic_field (CGS interpreted as Gauss, SI interpreted as Tesla)")
    print("  - bfield, b, gauss (aliases)")


def main():
    parser = build_parser()
    args = parser.parse_args()

    did_something = False

    if args.write_report:
        write_conversion_report(args.report_file)
        print(f"✅ Wrote conversion report to: {args.report_file}")
        did_something = True

    if args.command == "list":
        print_supported_quantities()
        did_something = True

    elif args.command == "convert":
        result = convert_value(
            value=args.value,
            quantity=args.quantity,
            from_unit=args.from_unit,
            to_unit=args.to_unit
        )
        print(f"{args.value:.16e} [{args.quantity}] in {args.from_unit.upper()} = {result:.16e} in {args.to_unit.upper()}")
        did_something = True

    if not did_something:
        # Friendly default behavior: write the report if no args are given
        write_conversion_report("conversion_factors_CU_CGS.txt")
        print(" No command given; wrote default conversion report: conversion_factors_CU_CGS.txt")
        print("   Use '--help' for CLI conversion options.")


if __name__ == "__main__":
    main()
