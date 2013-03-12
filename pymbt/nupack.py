# nupack 3.0 python wrapper. Loosely based on the nupack 2.0
# wrapper distributed as part of the RBS Calculator:
# (https://github.com/hsalis/Ribosome-Binding-Site-Calculator)

import multiprocessing
import time
from subprocess import call
from tempfile import mkdtemp
from shutil import rmtree
from os.path import isdir
from os import environ
from pymbt.sequence_manipulation import check_alphabet

# TODO: rmdir doesn't do anything?

if 'NUPACKHOME' in environ:
    NUPACKHOME = environ['NUPACKHOME']
else:
    raise EnvironmentError('NUPACKHOME environment variable must be set')


class Nupack:
    '''NUPACK class'''
    def __init__(self, sequence_list, material):
        # Set up nupack environment variable
        self.nupack_home = NUPACKHOME
        # Store reused information
        if type(sequence_list) != list:
            if type(sequence_list) == str:
                # if user input a string, make it a list for processing later
                sequence_list = [sequence_list]
            else:
                raise ValueError('Did not supply a sequence string or list')
        self.sequence_list = sequence_list
        self.outdir = mkdtemp()

        self.material = material
        if self.material == 'rna' or self.material == 'rna1999':
            mat = 'rna'
        elif material == 'dna':
            mat = 'dna'
        else:
            raise ValueError("material must be 'dna', 'rna', or 'rna1999'.")

        for sequence in self.sequence_list:
            check_alphabet(sequence, material=mat)

        self.complexes_run = False

    def complexes(self, max_complexes, mfe=True, rmdir=True, T=50):
        """
        Runs 'complexes' command on the inputted sequence list.

        'max_complexes': maximum complex size (integer).
        'mfe': include mfe calculations (boolean, defaults to True).
        'T': sets the temperature.
        """

        # This version can only handle max_complex 2 or below - I think
        self._open()  # Recreate temp dir if it was removed

        # Prepare input file
        complexes_input = '\n'.join([str(len(self.sequence_list)),
                                     '\n'.join(self.sequence_list),
                                     str(max_complexes)])

        f = open(self.outdir + '/nupack.in', 'w')
        f.write(complexes_input)
        f.close()

        # Prepare arguments
        if mfe:
            mfe = ' -mfe '
        else:
            mfe = ''

        # Calculate complexes
        complexes_args = ' -T %f -material %s %s' % (T, self.material, mfe)
        self._run_nupack('complexes', complexes_args)

        # Parse the output
        f_out = open(self.outdir + '/nupack.cx', 'r+')
        eq_results = f_out.readlines()

        eq_results = [x for x in eq_results if '%' not in x]
        eq_results = [v.split() for v in eq_results]

        # Separate results into complex type and energy
        # Remove rank entry
        for i in eq_results:
            i.pop(0)

        eq_results_cx = []
        for x in eq_results:
            cx_i = []
            for i, v in enumerate(self.sequence_list):
                cx_i.append(int(float(x.pop(0))))
            eq_results_cx.append(cx_i)

        eq_results_en = [float(x.pop(0)) for x in eq_results]

        self.complexes_run = max_complexes
        return {'complexes': eq_results_cx,
                'complex_energy': eq_results_en}

    def concentrations(self, max_complexes, conc=0.5e-6, mfe=True, rmdir=True):
        """
        Runs 'concentrations' command on the inputted sequence list.

        'conc': the concentration of each species. Defaults to 0.5e-6.
                Can be either a single value or list.
        'max_complexes': maximum complex size (integer).
        'mfe': include mfe calculations (boolean, defaults to True).
        """
        if self.complexes_run and self.complexes_run == max_complexes:
            pass  # complexes has already been run
        else:
            self._open()  # Recreate temp dir if it was removed
            # Run complexes in order to get input file
            self.complexes(max_complexes=max_complexes,
                           mfe=True,
                           rmdir=False)

        # Prepare input file
        if type(conc) == list:
            if len(conc) != len(self.sequence_list):
                raise ValueError("len(conc) must be len(sequences)")
            concstring = '\n'.join(conc)
        else:
            concstring = '\n'.join([str(conc) for x in self.sequence_list])

        f = open(self.outdir + '/nupack.con', 'w')
        f.write(concstring)
        f.close()

        # Run 'complexes'
        concentrations_args = '-sort 3 '
        self._run_nupack('concentrations', concentrations_args)

        # Parse the output of 'complexes'
        f_out = open(self.outdir + '/nupack.eq', 'r+')
        eq_results = f_out.readlines()
        eq_results = [x for x in eq_results if '%' not in x]
        eq_results = [v.split() for v in eq_results]

        # Format results: complex type, concentration, and energy
        # Remove rank information
        for i in eq_results:
            i.pop(0)

        eq_results_cx = []
        for x in eq_results:
            eq_cx_i = []
            for i, v in enumerate(self.sequence_list):
                eq_cx_i.append(int(float(x.pop(0))))
            eq_results_cx.append(eq_cx_i)

        eq_results_en = [x.pop(0) for x in eq_results]
        eq_results_conc = [float(x[0]) for x in eq_results]

        return {'types': eq_results_cx,
                'concentration': eq_results_conc,
                'energy': eq_results_en}

    def mfe(self, strand, T=50):
        """
        Runs 'mfe' command on the inputted sequence list.

        'T' sets the temperature.
        """

        self._open()  # Recreate temp dir if it was removed

        if type(strand) != int:
            raise TypeError('strand value must be integer')
        elif strand > len(self.sequence_list):
            raise ValueError('strand value exceeds the number of sequences')

        sequence = self.sequence_list[strand]
        # Prepare input file
        f = open(self.outdir + '/nupack.in', 'w')
        f.write(sequence)
        f.close()

        # Run 'mfe'
        mfe_args = ' -T ' + str(T) + ' -material ' + self.material
        self._run_nupack('mfe', mfe_args)

        # Parse the output of 'mfe'
        mfe_raw = open(self.outdir + '/nupack.mfe', 'r+').readlines()[14]
        mfe = float(mfe_raw.strip())

        # Return the mfe of the designated strand
        return mfe

    def pairs(self, strand, T=50):
        """
        Runs 'pairs' command on the inputted sequence list.

        'T' sets the temperature.
        """

        self._open()  # Recreate temp dir if it was removed

        sequence = self.sequence_list[strand]
        f = open(self.outdir + '/nupack.in', 'w')
        f.write(sequence)
        f.close()

        # Calculate pair probabilities with 'pairs'
        pairs_args = '-T ' + str(T) + ' -material ' + self.material
        self._run_nupack('pairs', pairs_args)

        # Parse the 'pairs' output
        # Only look at the last n rows - unbound probabilities
        pairs = open(self.outdir + '/nupack.ppairs', 'r+').readlines()
        pairs = pairs[-len(sequence):]
        pairs = [x for x in pairs if '%' not in x]
        #########################
        # Why was this code here?
        #pairs.pop(0) # Remove gap
        #pair_n = pairs.pop(0)
        #########################
        types = [(int(v.split()[0]), int(v.split()[1])) for v in pairs]
        pp = [float(v.split()[2]) for v in pairs]

        # Return the pair probabilities
        return {'type': types, 'probabilities': pp}

    def _open(self):
        # Create the input data file in a temp dir if it doesn't exist
        if not isdir(self.outdir):
            self.outdir = mkdtemp()

    def close(self):
        rmtree(self.outdir)

    def _run_nupack(self, cmd, arguments):
        known = ['complexes', 'concentrations', 'mfe', 'pairs']
        if cmd not in known:
            raise ValueError('Command must be one of the following:' + known)

        env_line = "export NUPACKHOME=" + self.nupack_home + ' && '
        command_prefix = '%s/bin/%s ' % (self.nupack_home, cmd)
        command_line = 'cd %s && ' % (self.outdir) + command_prefix
        arguments_line = arguments + ' ' + self.outdir + '/nupack'

        run_command = env_line + command_line + arguments_line
        call(run_command + ' > /dev/null', shell=True)


def nupack_multiprocessing(inputs, material, cmd, arguments, report=True):
    '''Provides access to NUPACK commands with multiprocessing support.
       Inputs: list of sequences
       material: 'dna' or 'rna'
       cmd: string of command - 'mfe', 'pairs', 'complexes', 'concentrations'
       arguments: keyword list of arguments
    '''
    nupack_pool = multiprocessing.Pool()
    try:
        args = [{'seq': x,
                 'cmd': cmd,
                 'material': material,
                 'arguments': arguments} for x in inputs]
        nupack_iterator = nupack_pool.imap(run_nupack, args)
        total = len(inputs)
        msg = ' calculations complete.'
        t = 4
        while report:
            completed = nupack_iterator._index
            if (completed == total):
                break
            else:
                if t >= 4:
                    print '(%s/%s) ' % (completed, total) + msg
                    t = 0
                t += 1
                time.sleep(1)
        multi_output = [x for x in nupack_iterator]
        nupack_pool.close()
        nupack_pool.join()
    except KeyboardInterrupt:
        print "Caught KeyboardInterrupt, terminating workers"
        nupack_pool.terminate()
        nupack_pool.close()
    return multi_output


def run_nupack(kwargs):
    run = Nupack(kwargs['seq'], material=kwargs['material'])
    output = getattr(run, kwargs['cmd'])(**kwargs['arguments'])
    run.close()
    return output
