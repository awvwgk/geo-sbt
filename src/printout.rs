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

extern crate ndarray;
use ndarray::prelude::*;

pub use crate::molecule::*;
pub use crate::elements::*;
pub use crate::econv::*;

pub fn print_geometry_information(mol : &Molecule, cn : &Array1<f64>,
    q : &Array1<f64>, molmass : &f64, cma : &Array1<f64>, moments : &Array1<f64>,
    dipole : &Array1<f64>)
{
    println!("{:>5} {:>3}    {:>10} {:>10} {:>11} {:>11} {:>11}",
             "#","Z","CN","q","x","y","z");
    for iter in (0..mol.n).zip(mol.xyz.axis_iter(Axis(0)).zip(mol.at.iter()
                        .zip(cn.iter().zip(q.iter()))))
    {
        let (i, (xyz, (at, (cni, qi)))) = iter;
        println!("{:5} {:3} {:2} {:10.5} {:10.5} {:10.6}", i+1,
            ordinal_number(&at), symbol(&at), cni, qi, xyz.mapv(|x|x*autoaa) );
    }
    println!("");
    println!("      molecular mass/u   : {:12.7}", molmass*autoamu);
    println!("   center of mass at/Å   : {:11.7}", cma*autoaa);
    println!("rotational constants/cm⁻¹: {:11.7}", moments.mapv(|x|0.5*autorcm/x));
    println!("       dipole moment/D   : {:11.7}", dipole*autod);
}

pub fn print_stretch_information(mol : &Molecule, stretch : &Array1<Bond>)
{
    println!("\n * {} selected distances\n", stretch.len());
    println!("{:>5} {:>3}    {:>5} {:>3}    {:>14} {:>14}",
             "#","Z","#","Z","distance/Å","tension/kcal·Å⁻²");
    for bond in stretch.iter()
    {
        let i = bond.iatom;
        let iat = &mol.at[i];
        let j = bond.jatom;
        let jat = &mol.at[j];
        println!("{:5} {:3} {:2} {:5} {:3} {:2} {:14.7} {:14.7}",
                 i+1, ordinal_number(&iat), symbol(&iat),
                 j+1, ordinal_number(&jat), symbol(&jat),
                 bond.length*autoaa,
                 bond.force*autokcal/autoaa/autoaa);
    }
}
pub fn print_bend_information(mol : &Molecule, bend : &Array1<Angle>)
{
    println!("\n * {} selected bends\n", bend.len());
    println!("{:>5} {:>3}    {:>5} {:>3}    {:>5} {:>3}    {:>14} {:>14}",
             "#","Z","#","Z","#","Z","angle/°","bend/kcal·Å⁻²");
    for angle in bend.iter()
    {
        let i = angle.iatom;
        let iat = &mol.at[i];
        let j = angle.jatom;
        let jat = &mol.at[j];
        let k = angle.katom;
        let kat = &mol.at[k];
        println!("{:5} {:3} {:2} {:5} {:3} {:2} {:5} {:3} {:2} {:14.7} {:14.7}",
                 i+1, ordinal_number(&iat), symbol(&iat),
                 j+1, ordinal_number(&jat), symbol(&jat),
                 k+1, ordinal_number(&kat), symbol(&kat),
                 angle.angle*180.0/pi,
                 angle.force*autokcal/autoaa/autoaa);
    }

}
pub fn print_torsion_information(mol : &Molecule, torsion : &Array1<Dihedral>)
{
    println!("\n * {} selected dihedral angles\n", torsion.len());
    println!("{:>5} {:>3}    {:>5} {:>3}    {:>5} {:>3}    {:>5} {:>3}    {:>14}",
             "#","Z","#","Z","#","Z","#","Z","angle/°");
    for angle in torsion.iter()
    {
        let i = angle.iatom;
        let iat = &mol.at[i];
        let j = angle.jatom;
        let jat = &mol.at[j];
        let k = angle.katom;
        let kat = &mol.at[k];
        let l = angle.katom;
        let lat = &mol.at[l];
        println!("{:5} {:3} {:2} {:5} {:3} {:2} {:5} {:3} {:2} {:5} {:3} {:2} {:14.7}",
                 i+1, ordinal_number(&iat), symbol(&iat),
                 j+1, ordinal_number(&jat), symbol(&jat),
                 k+1, ordinal_number(&kat), symbol(&kat),
                 l+1, ordinal_number(&lat), symbol(&lat),
                 angle.angle*180.0/pi);
    }
}


pub fn header()
{
    let printout =
        /*< < < < < < < < < < < < < > > > > > > > > > > > > >*/"
           -------------------------------------------------
          |                  G E O - S B T                  |
           -------------------------------------------------";
        /*< < < < < < < < < < < < < > > > > > > > > > > > > >*/
    println!("{}\n", printout);
}

pub fn gpl_license()
{
    let printout = "
   Copyright (C) 2019 S. Ehlert
   
   This program is free software: you can redistribute it and/or 
   modify it under the terms of the GNU General Public License as 
   published by the Free Software Foundation, either version 3 of 
   the License, or (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.
   If not, see <http://www.gnu.org/licenses/>.";

    println!("{}", printout);
}

