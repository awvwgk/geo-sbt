/*
 * This file is part of `geo-sbt'
 *
 * Copyright (C) 2019 S. Ehlert
 * 
 * This program is free software: you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License as 
 * published by the Free Software Foundation, either version 3 of 
 * the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.
 * If not, see <http://www.gnu.org/licenses/>.";
 */

use std::io::BufRead;

extern crate structopt;
use structopt::StructOpt;

extern crate ndarray;
use ndarray::prelude::*;

extern crate ndarray_linalg;
use ndarray_linalg::lapack_traits::*;
use ndarray_linalg::eigh::*;
use ndarray_linalg::Solve;

extern crate statrs;
use statrs::function::erf::*;

#[derive(StructOpt)]
struct Args
{
    #[structopt(long)]
    /// print GPL license and exit
    license: bool,
    #[structopt(long)]
    /// obtain partial charges from classical model density
    ///
    /// calculate partial charges from Electronegativity Equilibrium (EEQ) model
    /// as in
    /// E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher,
    /// C. Bannwarth and S. Grimme, ChemRxiv, 2018, preprint.
    /// https://doi.org/10.26434/chemrxiv.7430216.v2
    charges: bool,
    #[structopt(long)]
    /// find relevant bonds
    stretch: bool,
    #[structopt(long,default_value = "0.350")]
    /// force constant for stretches in au
    kstretch: f64,
    #[structopt(long)]
    /// find relevant bond angles
    bend: bool,
    #[structopt(long,default_value = "0.150")]
    /// force constant for stretches in au
    kbend: f64,
    #[structopt(long)]
    /// find relevant dihedral angles
    torsion: bool,
    #[structopt(long,default_value = "0.005")]
    /// force constant for stretches in au
    ktorsion: f64,
    #[structopt(long)]
    #[structopt(short,long,default_value = "1.2")]
    /// threshold in multiples of the covalent radius for determining bonds
    thr: f64,
    #[structopt(short,long,default_value = "0.0")]
    /// molecular charge of the system in e⁻
    chrg: f64,
    #[structopt(parse(from_os_str))]
    /// geometry input file for calculation
    path: std::path::PathBuf,
}

// Enumerator for PSE and related functions
mod elements;
use crate::elements::*;

// charge model parameters
mod chrgeq;
use crate::chrgeq::*;

// unit conversion factors and important constants
mod econv;
use crate::econv::*;

mod printout;
use crate::printout::*;

// define molecule type and ways to obtain basic properties
mod molecule;
use crate::molecule::*;

const debug : bool = true;

// parser for xmol formatted geometry input
fn read_geometry(file : std::fs::File)
    -> Result<(Molecule), Box<dyn std::error::Error>>
{
    let mut reader = std::io::BufReader::new(file);

    let mut line = String::new();
    reader.read_line(&mut line)?;
    //if debug { println!("First line read {:?}", line); }
    let nat = line.trim().parse::<usize>()?;
    //if debug { println!("First line parsed {:?}", nat); }
    let mut xyz = Vec::new();
    let mut atoms = Vec::new();

    reader.read_line(&mut line)?;
    //if debug { println!("Second line read {:?} (skipped)", line); }

    for line in reader.lines()
    {
        //if debug { println!("{:?}", line); }
        let line = line?;
        let mut parts = line.split_whitespace();

        let sym = match parts.next()
           { Some(value) => value.to_string(),     None => continue };
        let x   = match parts.next()
           { Some(value) => value.parse::<f64>()?, None => continue };
        let y   = match parts.next()
           { Some(value) => value.parse::<f64>()?, None => continue };
        let z   = match parts.next()
           { Some(value) => value.parse::<f64>()?, None => continue };

        let at = element(sym);
        //if debug { println!("{:2} {:?} {:?} {:?}", symbol(&at), x, y, z); }
        atoms.push(at);
        xyz.push(x*aatoau); xyz.push(y*aatoau); xyz.push(z*aatoau);
    }
    let xyz = Array::from_shape_vec((atoms.len(),3),xyz)?;
    //if debug { println!("{:?}", xyz); }

    Ok(Molecule{n : atoms.len(), at : atoms, xyz : xyz})
}

// fractional coordination number
fn coordination_number(mol : &Molecule)
    -> Result<Array1<f64>, Box<dyn std::error::Error>>
{
    let k2 = 4.0/3.0;
    let kn = 7.5;
    let mut cn = Vec::new();
    for iatom in mol.xyz.axis_iter(Axis(0)).zip(mol.at.iter())
    {
        let mut cni = 0.0;

        let (ixyz, iat) = iatom;
        let rcovi = covalent_radius(&iat) * k2;

        for jatom in mol.xyz.axis_iter(Axis(0)).zip(mol.at.iter())
        {
            let (jxyz, jat) = jatom;
            //if debug {println!("{:?} {:?}",iatom,jatom);}
            if ixyz == jxyz { continue; }
            let rcovj = covalent_radius(&jat) * k2;
            let r = get_distance(&ixyz, &jxyz);

            let rcovij = rcovi + rcovj;
            
            let count = 0.5 * (1.0 + erf(-kn * (r/rcovij - 1.0))); 
            //if debug {println!("{:12.7} {:12.7} {:12.7} +{:12.7}",r,rcovij,count,cni);}
            cni += count;
        }
        //if debug {println!("-> {}",cni);}
        cn.push(cni);
    }

    Ok(Array1::<f64>::from_vec(cn))
}

fn eeq_chrgeq(mol : &Molecule, chrg : &f64, cn : &Array1<f64>)
    -> Result<(Array1<f64>), Box<dyn std::error::Error>>
{
    let mut amat = Array::<f64, _>::zeros((mol.n+1,mol.n+1));
    let mut xvec = Array::<f64, _>::zeros((mol.n+1));
    for iter in (0..mol.n).zip(mol.xyz.axis_iter(Axis(0)).zip(mol.at.iter()))
    {
        let (i, (ixyz, iat)) = iter;
        let ipar = get_charge_parameter(&iat);
        let alp2 = f64::powi(ipar.alpha,2);
        for jter in (0..i).zip(mol.xyz.axis_iter(Axis(0)).zip(mol.at.iter()))
        {
            let (j, (jxyz, jat)) = jter;
            let jpar = get_charge_parameter(&jat);
            let r = get_distance(&ixyz, &jxyz);
            let gamij = 1.0/f64::sqrt(alp2+f64::powi(jpar.alpha,2));
            let tmp = erf(gamij*r)/r;
            amat[[i,j]] = tmp;
            amat[[j,i]] = tmp;
        }
        xvec[i] = -ipar.en + ipar.kappa * f64::sqrt(cn[i]);
        amat[[i,i]] = ipar.gam + f64::sqrt(2.0/pi)/ipar.alpha;
        amat[[i,mol.n]] = 1.0;
        amat[[mol.n,i]] = 1.0;
    }
    amat[[mol.n,mol.n]] = 0.0;
    xvec[mol.n] = *chrg;

    let q = amat.solve_into(xvec)?;

    Ok(q)
}

fn get_dipole_from_charges(mol : &Molecule, q : &Array1<f64>)
    -> Result<Array1<f64>, Box<dyn std::error::Error>>
{
    let mut dipole = Array::<f64, _>::zeros((3));
    for iter in mol.xyz.axis_iter(Axis(0)).zip(q.iter())
    {
        let (xyz, qi) = iter;
        dipole += &xyz.mapv(|x|x * qi);
    }

    Ok(dipole)
}

fn main()
    -> Result<(), Box<dyn std::error::Error>>
{
    let args = Args::from_args();
    if args.license { gpl_license(); return Ok(()) }

    header();
    //gpl_license();

    let file = std::fs::File::open(&args.path)?;
    let mol = read_geometry(file)?;

    let cn = coordination_number(&mol)?;

    let molmass = total_molecular_mass(&mol)?;
    let cma     = center_of_mass(&mol)?;
    let moments = moments_of_interia(&mol, &cma)?;

    let q      = eeq_chrgeq(&mol, &args.chrg, &cn)?;
    let dipole = get_dipole_from_charges(&mol, &q)?;

    print_geometry_information(&mol, &cn, &q, &molmass, &cma, &moments, &dipole);

    if args.stretch
    {
        let stretch = find_all_stretches(&mol, &args.kstretch)?;
        print_stretch_information(&mol, &stretch);
    }
    if args.bend || args.torsion
    {
        let bonds = get_bonds(&mol, &args.thr)?;
        if args.bend
        {
            let bend = find_all_bends(&mol, &bonds, &args.kbend)?;
            print_bend_information(&mol, &bend);
        }
        if args.torsion
        {
            let torsion = find_all_torsions(&mol, &bonds)?;
            print_torsion_information(&mol, &torsion);
        }
    }


    Ok(())
}
