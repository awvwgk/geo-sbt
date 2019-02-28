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

extern crate ndarray_linalg;
use ndarray_linalg::lapack_traits::*;
use ndarray_linalg::eigh::*;

pub use crate::elements::*;
pub use crate::econv::*;

// standard wrapper type for molecule information used in this program
pub type Coord = Array2<f64>;
pub struct Molecule
{
    // number of atoms == dimension
    pub n   : usize,
    // identifier for element type
    pub at  : Vec<Elem>,
    // cartesian coordinates in Bohr
    pub xyz : Coord,
}

pub fn total_molecular_mass(mol : &Molecule)
    -> Result<f64, Box<dyn std::error::Error>>
{
    let mut mass = 0.0;
    for atom in mol.at.iter()
    {
        //if debug { println!("{}", symbol(&atom)); }
        mass += atommass(&atom);
    }
    //if debug { println!("{}", mass * autoamu); }

    Ok(mass)
}

pub fn center_of_mass(mol : &Molecule)
    -> Result<(Array1<f64>), Box<dyn std::error::Error>>
{
    let mut mass = 0.0;
    let mut cma  = arr1(&[0.0,0.0,0.0]);
    for atom in mol.xyz.axis_iter(Axis(0)).zip(mol.at.iter())
    {
        let (xyz, at) = atom;
        let atmass = atommass(&at);
        cma += &xyz.mapv(|x| x*atmass);
        //if debug { println!("{} {:?}", symbol(&at), cma); }
        mass += atmass;
    }
    let cma = cma / mass;
    //if debug { println!("{} {:?}", mass * autoamu, &cma * autoaa); }

    Ok(cma)
}

// get the moments of intertia for a given molecular geometry
// point should be the center of mass of the molecule
pub fn moments_of_interia(mol : &Molecule, point : &Array1<f64>)
    -> Result<(Array1<f64>), Box<dyn std::error::Error>>
{
    let mut inertia = Array2::<f64>::zeros((3,3));

    for atom in mol.xyz.axis_iter(Axis(0)).zip(mol.at.iter())
    {
        let (xyz, at) = atom;
        let atmass = atommass(&at);
        let x = xyz[(0)] - point[(0)]; let x2 = x*x;
        let y = xyz[(1)] - point[(1)]; let y2 = y*y;
        let z = xyz[(2)] - point[(2)]; let z2 = z*z;
        inertia[(0,0)] += atmass * (y2 + z2);
        inertia[(1,0)] += atmass * x * y;
        inertia[(2,0)] += atmass * x * z;
        inertia[(1,1)] += atmass * (x2 + z2);
        inertia[(2,1)] += atmass * y * z;
        inertia[(2,2)] += atmass * (x2 + y2);
    }
    inertia[(0,1)] = inertia[(1,0)];
    inertia[(0,2)] = inertia[(2,0)];
    inertia[(1,2)] = inertia[(2,1)];

    //if debug { println!("{}", inertia); }

    let moments = inertia.eigvalsh(UPLO::Upper)?;
    //if debug { println!("{:?}", moments); }

    Ok(moments)
}

pub struct Bond
{
    pub iatom  : usize,
    pub jatom  : usize,
    pub length : f64,
    pub force  : f64,
}

pub fn get_distance(ixyz : &ArrayView1<f64>,
                    jxyz : &ArrayView1<f64>)
    -> f64
{
    let rij : Array1<f64> = jxyz - ixyz;
    f64::sqrt(rij.dot(&rij))
}

fn fk_swart(alpha : &f64,c : &f64,r : &f64)
    -> f64
{
    f64::exp(-alpha*(r/c - 1.0))
}

pub fn find_all_stretches(mol : &Molecule, kr : &f64)
    -> Result<Array1<Bond>, Box<dyn std::error::Error>>
{
    let mut stretch = Vec::new();
    for iter in (0..mol.n).zip(mol.xyz.axis_iter(Axis(0)).zip(mol.at.iter()))
    {
        let (i, (ixyz, iat)) = iter;
        let rcovi = covalent_radius(&iat);
        for jter in (0..i).zip(mol.xyz.axis_iter(Axis(0)).zip(mol.at.iter()))
        {
            let (j, (jxyz, jat)) = jter;
            let rcovj = covalent_radius(&jat);
            let rcovij = 1.2*(rcovi + rcovj);
            let r = get_distance(&ixyz, &jxyz);
            if r < rcovij
            {
                let c = rcovi + rcovj;
                let gij = kr * fk_swart(&1.0,&c,&r);
                stretch.push(Bond{iatom:i,jatom:j,length:r,force:gij});
            }
        }
    }
    Ok(Array1::<Bond>::from_vec(stretch))
}

pub fn get_bonds(mol : &Molecule, scale : &f64)
    -> Result<Array2<usize>, Box<dyn std::error::Error>>
{
    let mut bonds = Array::<usize, _>::zeros((mol.n,mol.n));
    for iter in (0..mol.n).zip(mol.xyz.axis_iter(Axis(0)).zip(mol.at.iter()))
    {
        let (i, (ixyz, iat)) = iter;
        let rcovi = covalent_radius(&iat);
        for jter in (0..i).zip(mol.xyz.axis_iter(Axis(0)).zip(mol.at.iter()))
        {
            let (j, (jxyz, jat)) = jter;
            let rcovj = covalent_radius(&jat);
            let rcovij = scale*(rcovi + rcovj);
            let r = get_distance(&ixyz, &jxyz);
            if r < rcovij
            {
                bonds[[i,i]] += 1;
                bonds[[j,j]] += 1;
                bonds[[j,i]] += 1;
                bonds[[i,j]] += 1;
            }
        }
    }
    Ok(bonds)
}

pub struct Angle
{
    pub iatom : usize,
    pub jatom : usize,
    pub katom : usize,
    pub angle : f64,
    pub force : f64,
}

pub fn get_angle(ixyz : &ArrayView1<f64>,
                 jxyz : &ArrayView1<f64>,
                 kxyz : &ArrayView1<f64>)
    -> f64
{
    /*
    let rij : Array1<f64> = jxyz - ixyz;
    let nij = f64::sqrt(rij.dot(&rij));
    let rjk : Array1<f64> = kxyz - jxyz;
    let njk = f64::sqrt(rjk.dot(&rjk));

    let rij = rij/nij;
    let rjk = rjk/njk;

    let rijk = rij.dot(&rjk);

    if rijk >  1.0 { return f64::acos( 1.0); }
    if rijk < -1.0 { return f64::acos(-1.0); }

    f64::acos(rijk)
    */
    let rij : Array1<f64> = jxyz - ixyz;
    let r2ij = rij.dot(&rij);
    let rjk : Array1<f64> = kxyz - jxyz;
    let r2jk = rjk.dot(&rjk);
    let rik : Array1<f64> = kxyz - ixyz;
    let r2ik = rik.dot(&rik);
    let rijk = 0.5*(r2ij+r2jk-r2ik)/f64::sqrt(r2ij*r2jk);

    if rijk >  1.0 { return f64::acos( 1.0); }
    if rijk < -1.0 { return f64::acos(-1.0); }

    f64::acos(rijk)
}

pub fn find_all_bends(mol : &Molecule, bonds : &Array2::<usize>, kr : &f64)
    -> Result<Array1<Angle>, Box<dyn std::error::Error>>
{
    let mut bend = Vec::new();
    for iter in (0..mol.n).zip(mol.xyz.axis_iter(Axis(0)).zip(mol.at.iter()))
    {
        let (i, (ixyz, iat)) = iter;
        if bonds[[i,i]] < 1 { continue; }
        let rcovi = covalent_radius(&iat);
        for jter in (0..mol.n).zip(mol.xyz.axis_iter(Axis(0)).zip(mol.at.iter()))
        {
            let (j, (jxyz, jat)) = jter;
            if i == j { continue; }
            if bonds[[j,j]] < 2 { continue; }
            if bonds[[i,j]] < 1 { continue; }
            let rcovj = covalent_radius(&jat);
            let cij = rcovi + rcovj;
            let gij = fk_swart(&1.0,&cij,&get_distance(&ixyz,&jxyz));
            for kter in (0..mol.n).zip(mol.xyz.axis_iter(Axis(0)).zip(mol.at.iter()))
            {
                let (k, (kxyz, kat)) = kter;
                if (i == k) || (j == k) { continue; }
                if bonds[[k,k]] < 1 { continue; }
                if bonds[[j,k]] < 1 { continue; }
                let rcovk = covalent_radius(&kat);
                let cjk = rcovj + rcovk;
                let gjk = fk_swart(&1.0,&cjk,&get_distance(&jxyz,&kxyz));
                let gik = kr * gij * gjk;
                let phi = get_angle(&ixyz, &jxyz, &kxyz);
                bend.push(Angle{iatom:i,jatom:j,katom:k,angle:phi,force:gik});
            }
        }
    }
    Ok(Array1::<Angle>::from_vec(bend))
}

fn cross_product(ri : &ArrayView1<f64>,
                 rj : &ArrayView1<f64>)
    -> Array1<f64>
{
    let x = ri[1]*rj[2] - ri[2]*rj[1];
    let y = ri[2]*rj[0] - ri[0]*rj[2];
    let z = ri[0]*rj[1] - ri[1]*rj[0];

    arr1(&[x,y,z])
}

fn vec_determinant(ri : &ArrayView1<f64>,
                   rj : &ArrayView1<f64>,
                   rk : &ArrayView1<f64>)
    -> f64
{
    let det = ri[0]*(rj[1]*rk[2]-rj[2]*rk[1])
             -ri[1]*(rj[0]*rk[2]-rj[2]*rk[0])
             +ri[2]*(rj[0]*rk[1]-rj[1]*rk[0]);
    det
}

pub fn get_dihedral(ixyz : &ArrayView1<f64>,
                    jxyz : &ArrayView1<f64>,
                    kxyz : &ArrayView1<f64>,
                    lxyz : &ArrayView1<f64>)
    -> f64
{
    let rij : Array1<f64> = jxyz - ixyz;
    let nij = f64::sqrt(rij.dot(&rij));
    let rjk : Array1<f64> = kxyz - jxyz;
    let njk = f64::sqrt(rjk.dot(&rjk));
    let rkl : Array1<f64> = lxyz - kxyz;
    let nkl = f64::sqrt(rkl.dot(&rkl));
    // normalize all distances
    let rij = rij/nij;
    let rjk = rjk/njk;
    let rkl = rkl/nkl;

    let det = vec_determinant(&rij.view(),&rjk.view(),&rkl.view());
    let rji = -&rij;
    let tst = cross_product(&rji.view(),&rkl.view());
    let tst = tst.dot(&rjk);

    let rijk = cross_product(&rij.view(),&rjk.view());
    let rjkl = cross_product(&rjk.view(),&rkl.view());
    let rijkl = rijk.dot(&rjkl);

    if rijkl >  1.0 { return f64::acos( 1.0); }
    if rijkl < -1.0 { return f64::acos(-1.0); }

    let phi = f64::acos(rijkl);
    if tst < 0.0
    { 
        -phi
    }
    else
    {
        phi
    }
}

pub struct Dihedral
{
    pub iatom : usize,
    pub jatom : usize,
    pub katom : usize,
    pub latom : usize,
    pub angle : f64,
}

pub fn find_all_torsions(mol : &Molecule, bonds : &Array2::<usize>)
    -> Result<Array1<Dihedral>, Box<dyn std::error::Error>>
{
    let mut torsion = Vec::new();
    for iter in (0..mol.n).zip(mol.xyz.axis_iter(Axis(0)))
    {
        let (i, ixyz) = iter;
        if bonds[[i,i]] < 1 { continue; }
        for jter in (0..mol.n).zip(mol.xyz.axis_iter(Axis(0)))
        {
            let (j, jxyz) = jter;
            if i == j { continue; }
            if bonds[[j,j]] < 2 { continue; }
            if bonds[[i,j]] < 1 { continue; }
            for kter in (0..mol.n).zip(mol.xyz.axis_iter(Axis(0)))
            {
                let (k, kxyz) = kter;
                if (i == k) || (j == k) { continue; }
                if bonds[[k,k]] < 2 { continue; }
                if bonds[[j,k]] < 1 { continue; }
                for lter in (0..mol.n).zip(mol.xyz.axis_iter(Axis(0)))
                {
                    let (l, lxyz) = lter;
                    if (i == l) || (j == l) || (k == l) { continue; }
                    if bonds[[l,l]] < 1 { continue; }
                    if bonds[[k,l]] < 1 { continue; }
                    let phi = get_dihedral(&ixyz,&jxyz,&kxyz,&lxyz);
                    torsion.push(Dihedral{iatom:i,jatom:j,katom:k,latom:l,angle:phi});
                }
            }
        }
    }

    Ok(Array1::<Dihedral>::from_vec(torsion))
}



