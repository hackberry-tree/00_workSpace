#! /usr/bin/env python


# vasp_unfold
#
# Author: Tomic Milan
#
# vasp_unfold is a code whose purpose is to perform the unfolding
# of band structures calculated with VASP. It is based on the method
# described in (reference in preparation). If this code have been
# used to obtain the data in your publication please cite the given
# reference.
#
# This code is provided as is. There is no guarantee that it will work.
# However, if you encounter errors or unexpected behavor, please report
# it to tomic@itp.uni-frankfurt.de and I will try to improve it.


import numpy as np
import argparse
import sys


def post_error(error_info):
    sys.stderr.write('\nError: '+error_info+'\n\n')

    exit()


def translation(tstring):
    try:
        tstring = tstring.replace('0', ' 1/1')

        return np.array([int(s.split('/')[1]) for s in tstring.split(',')])
    except:
        post_error('Unable to parse string: "{0}". Check help for valid '
                   'translation generator specification'.format(tstring))



def gcd(a, b):
    '''Return greatest common divisor using Euclid's Algorithm.'''
    while b:
        a, b = b, a % b
    return a


def lcm(a, b):
    '''Return lowest common multiple.'''
    return a * b // gcd(a, b)


def lcmm(*args):
    '''Return lcm of args.'''
    return reduce(lcm, args)




def parse_poscar(poscar):
    '''Parse the VASP's POSCAR file and return the unit cell
    and the list of fractional coordinates.
    '''
    try:
        plines = open(poscar).readlines()
    except:
        post_error('Unable to open "{0}" for reading.'.format(poscar))

    # Read the scaling factor
    f = float(plines[1].strip())

    # Read the unit cell vectors
    cell = f*np.array([l.strip().split() for l in plines[2:5]], float)

    # Read the chemical symbols
    syms = plines[5].strip().split()

    # Read the atom counts
    counts = np.array(plines[6].strip().split(), int)

    # Rearrange into the long list of chemical symbols
    symbols = []

    for s, c in zip(syms, counts):
        symbols += [s]*c

    # Cartesian or fractional coordinates?
    if l[7].lower() == 'c':
        mult = np.linalg.inv(cell)
    else:
        mult = np.eye(3)

    spos = np.zeros((len(symbols), 3))

    i = 0

    for l in plines[8:]:
        l = l.strip().split()

        if len(l) != 3:
            continue

        # Read the position
        spos[i] = np.dot(np.array(l[:3], float), mult)

        i += 1

    return cell, spos


def build_operators(spos, trans, eps=1e-6):
    '''Given a list for fractional positions, and a list of
    fractional translations, produce a set of matrices,
    describing how fractional translations permute atomic
    positions within the unit cell. Two fractional positions
    si and sj are considered to be identical when |si-sj|<eps.
    '''
    ntrans = len(trans)
    natoms = len(spos)

    ops = np.zeros((ntrans, natoms, natoms), int)

    for i, ti in enumerate(trans):
        for j, sj in enumerate(spos):
            for k, sk in enumerate(spos):
                # If displacement between two atomic positions
                # differs from the fractional translations by
                # by a lattice translation+-eps we consider the
                # that two atomic positions to be map onto each
                # other by the fractional translation
                disp = sk-sj-ti
                ops[i, j, k] = np.linalg.norm(disp-np.rint(disp)) < eps

    # Every row and every column of every operator must
    # contain exactly one unity, otherwise translations
    # are not maping atoms one-to-one.
    if np.any(np.sum(ops, axis=1) != 1) or np.any(np.sum(ops, axis=2) != 1):
        post_error('Translations are not one-to-one. '
            'Try changing the matching tolerance, or try using '
            'the POSCAR file with more regular positions.')

    return ops


def build_translations(tgens):
    '''Build a list of translations and irreps from at most
    three linearly independent generators specified in the
    form [nx,ny,nz], representing translation [1/nx,1/ny,1/nz].
    '''
    if len(tgens) > 3:
        post_error('There can be at most three generators '
            'of fractional translations.')

    # Check if generators are linearly independent
    if len(tgens) == 2 and np.all(np.cross(tgens[0], tgens[1]) == 0):
        post_error('Generators are not linearly independant.')
    elif len(tgens) == 3 and np.linalg.det(tgens) == 0:
        post_error('Generators are not linearly independant.')

    # Expand the generator list to have a 3x3 matrix
    tgens = np.append(tgens, np.ones((3-len(tgens), 3)), axis=0)

    # Get the order of every generator
    order = np.array([lcmm(*g) for g in tgens], int)

    # Get the corresponding roots of unity which will be
    # used to construct irreps
    deltag = np.exp(-2*np.pi*1j/order)

    # Calculate translation vectors
    tgens = (1.0/tgens)%1

    # Total number of translations
    ntrans = np.prod(order)

    # Storage for translations
    trans = np.zeros((ntrans, 3), float)

    # Calculate irreps of every generator
    irrepg = np.zeros((ntrans, 3), complex)

    # Loop over all possible products of powers of generators
    # and store the irreps
    for i in xrange(order[0]):
        for j in xrange(order[1]):
            for k in xrange(order[2]):
                ii = i*order[1]*order[2]+j*order[2]+k

                irrepg[ii] = deltag**[i, j, k]

    # Irrep table for all translations
    irreps = np.zeros((ntrans, ntrans), complex)

    # Loop over all posible products of powers of generators
    for i in xrange(order[0]):
        for j in xrange(order[1]):
            for k in xrange(order[2]):
                ti = i*order[1]*order[2]+j*order[2]+k

                # Get the translation vector
                trans[ti] = np.dot([i, j, k], tgens)%1

                # Get all irreps of that translations
                for l in xrange(ntrans):
                    irreps[l, ti] = np.prod(irrepg[l]**[i, j, k])

    return trans, irreps


def build_projectors(irreps, ops, intdim=1):
    '''Builds irrep projection operators, given irreps
    and corresponding translation operator matrices.
    Intdim specifies how many orbitals per atomic site
    there are.
    '''
    projs = np.zeros((len(irreps), ops.shape[1], ops.shape[2]), complex)

    # Loop over irreps
    for i, ii in enumerate(irreps):
        # Sum operators multiplied by their irrep
        for j, oj in enumerate(ops):
            projs[i] += irreps[i, j]*oj

    # Expand onto the orbital space and normalize
    return np.kron(projs, np.eye(intdim))/len(ops)


def parse_procar(filename):
    try:
        procar = open(filename)
    except:
        post_error('Unable to open "{0}" for reading.'.format(filename))

    # Read the first line
    line = procar.readline()

    # Check if phase information is included
    if 'phase' not in line:
        post_error('Phase information is needed in the input PROCAR file. '
            'Please specify LORBIT=12 in the INCAR file and recalculate '
            'the bandstructure.')

    # Read the second line and split it
    line = procar.readline().strip().split()

    npoints = int(line[3])
    nbands = int(line[7])
    nions = int(line[11])

    # Skip one line
    procar.readline()

    # We are now at the first k-point data block
    # Remeber the position
    last_pos = procar.tell()

    # Skip four lines
    procar.readline()
    procar.readline()
    procar.readline()
    procar.readline()

    # Figure out the number of the orbitals from the
    # column titles of the first k-point block
    line = procar.readline()

    orbitals = np.array(line.strip().split()[1:-1])

    norbs = len(orbitals)

    # Go back to the beginning of the k-point block
    procar.seek(last_pos)

    # Allocate the storage
    kpoints = np.zeros((npoints, 3), float)
    kweights = np.zeros(npoints, float)
    bands = np.zeros((npoints, nbands), float)
    occupations = np.zeros((npoints, nbands), float)
    weights = np.zeros((npoints, nions*norbs, nbands), complex)

    for i in xrange(npoints):
        # Read the k-point and its weight
        line = procar.readline()
        lspl = line.strip().split()

        kpoints[i] = np.array(lspl[3:6], float)
        kweights[i] = float(lspl[-1])

        # Skip the empty line
        procar.readline()

        for j in xrange(nbands):
            # Read the band energy and the occupations
            line = procar.readline()
            lspl = line.strip().split()

            bands[i, j] = float(lspl[4])
            occupations[i, j] = float(lspl[-1])

            # Skip the empty line
            procar.readline()

            # Skip the orbit labels
            procar.readline()

            # Skip the weights
            for k in xrange(nions):
                procar.readline()

            # Skip the totals
            procar.readline()

            # Skip the orbit labels
            procar.readline()

            for k in xrange(nions):
                # Get the real part
                line = procar.readline()
                weights[i, k*norbs:(k+1)*norbs, j] = \
                    np.array(line.split()[1:1+norbs], float)

                # Get the imaginary parth
                line = procar.readline()
                weights[i, k*norbs:(k+1)*norbs, j] += \
                np.array(line.split()[1:1+norbs], float)*1j

            # Skip the empty line
            procar.readline()

        #Skip the empty line
        procar.readline()

    return orbitals, kpoints, kweights, bands, occupations, weights


def write_procar(fname, orbitals, kpoints, kweights, bands,
                 occupations, weights):
    # Labels for orbitals
    orblabels = ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2',
                 'f-3', 'f-2', 'f-1', 'f0', 'f1', 'f2', 'f3']

    norb = len(orbitals)
    npoints = len(kpoints)
    nbands = bands.shape[1]
    nions = weights.shape[1]/norb

    try:
        out = open(fname, 'w')
    except:
        post_error('Unable to open "{0}" for writing'.format(fname))

    # Write the first line of the PROCAR file
    out.write('PROCAR lm decomposed + phase\n')

    # Write the second line containing the sizes
    out.write('# of k-points:  {0}         # of bands:  {1}'
              '         # of ions:   {2}\n\n'.format(npoints, nbands, nions))

    # Format out the column title line for orbital weights
    orb_ttl_1 = 'ion '
    orb_ttl_2 = 'ion '

    for i in xrange(norb):
        orb_ttl_1 += '{0: >6} '.format(orblabels[i])
        orb_ttl_2 += '{0: >6} '.format(orblabels[i])

    orb_ttl_1 += '{0: >6}\n'.format('tot')
    orb_ttl_2 += '\n'

    for i, k in enumerate(kpoints):
        # Write the k-point info
        out.write(' k-point {0: >4} :    '.format(i+1))
        out.write('{0:.8f} {1:.8f} {2:.8f}     '.format(*k))
        out.write('weight = {0:.8f}\n\n'.format(kweights[i]))

        for j, b in enumerate(bands[i]):
            # Write the band info
            out.write('band {0: >4} # '.format(j+1))
            out.write('energy {0: >13.8f} # '.format(b))
            out.write('occ. {0: >11.8f}\n\n'.format(occupations[i, j]))

            out.write(orb_ttl_1)

            tot_orb = np.zeros(norb, float)

            for k in xrange(nions):
                out.write('{0: >3} '.format(k+1))

                cw = weights[i, k*norb:(k+1)*norb, j]
                aw = np.abs(cw)**2

                tot_orb += aw

                for l in xrange(norb):
                    out.write('{0: >6.3f} '.format(aw[l]))

                out.write('{0: >6.3f}\n'.format(np.sum(aw)))

            out.write('tot ')

            for l in xrange(norb):
                out.write('{0: >6.3f} '.format(tot_orb[l]))

            out.write('{0: >6.3f}\n'.format(np.sum(tot_orb)))

            out.write(orb_ttl_2)

            for k in xrange(nions):
                out.write('{0: >3} '.format(k+1))

                cw = weights[i, k*norb:(k+1)*norb, j]

                for l in xrange(norb):
                    out.write('{0: >6.3f} '.format(cw[l].real))

                out.write('\n{0: >3} '.format(k+1))

                for l in xrange(norb):
                    out.write('{0: >6.3f} '.format(cw[l].imag))

                out.write('\n')

            out.write('\n')

        out.write('\n')

    out.close()



def main():
    desc_str = 'Unfold bands calculated by VASP. For this, phase '\
               'information needs to be present in the PROCAR file '\
               'which means that bands need to be calculated with '\
               'the option LORBIT=12 specified in the INCAR file. '\
               'The required input consists of POSCAR file, PROCAR '\
               'file and a set of up to three fractional translations. '\
               'There is no need to specify all translations, just the '\
               'generators. The programs will generate all distinct '\
               'translations from the generators and generate all irreps. '


    parser = argparse.ArgumentParser(prog='vasp_unfold', description = desc_str)

    parser.add_argument('poscar', type=str, help='POSCAR file')
    parser.add_argument('procar', type=str, help='PROCAR file')

    parser.add_argument('--tgen', type=translation, action='append',
                        metavar='SX,SY,SZ', help='Fractional translation '
                        'generator. No whitespaces are allowed between the '
                        'components! SX, SY and SZ can be either 0 or 1/n, '
                        'where n is an integer. Up to three linearly '
                        'independant generators can be specified.')

    parser.add_argument('--out', type=str, help='Output filename. If left '
                        'unspecified  output is writen to PROCAR.irrep.n '
                        'where PROCAR is location of the input PROCAR file. '
                        'If specified, output is written to OUT.irrep.n')

    parser.add_argument('--eps', type=float, default=1e-6, help='Numerical '
                        'precision. When building permutation representation '
                        'of the fractional translations this parameter is used '
                        'to determine whether two fractional positions are '
                        'identical. For irregular structures, it may need to '
                        'be increased. If this does not help, try to tweak '
                        'the atomic positions in POSCAR to make the structure '
                        'more regular. Default is 1e-6')

    args = parser.parse_args()

    tgens = args.tgen

    trans, irreps = build_translations(tgens)

    try:
        cell, spos = parse_poscar(args.poscar)
    except:
        post_error('Unable to parse the input POSCAR file. Please '
            'check if the file exists and is formatted properly')

    ops = build_operators(spos, trans, args.eps)

    try:
        data = parse_procar(args.procar)
    except:
        post_error('Unable to parse the input PROCAR file. Please '
            'check if the file exists and is formatted properly.')

    weights = np.copy(data[-1])

    norbs = weights.shape[1]/len(spos)

    projs = build_projectors(irreps, ops, norbs)

    if args.out is None:
        output = args.procar
    else:
        output = args.out

    for i, p in enumerate(projs):
        for j in xrange(data[-1].shape[0]):
            try:
                data[-1][j] = np.dot(p, weights[j])
            except:
                post_error('Unable to apply projectors. Are you sure '
                    'that specified POSCAR and PROCAR file belong to '
                    'the same crystal structure?')

        write_procar('{0}.irrep.{1}'.format(output, i), *data)


if __name__ == '__main__':
    main()

