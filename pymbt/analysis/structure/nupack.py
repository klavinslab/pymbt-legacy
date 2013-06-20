'''Wrapper for NUPACK 3.0.'''

import multiprocessing
import time
from subprocess import call
from tempfile import mkdtemp
from shutil import rmtree
from os.path import isdir
from os import environ
from pymbt.sequence.utils import check_instance


if 'NUPACKHOME' in environ:
    NUPACKHOME = environ['NUPACKHOME']
else:
    raise EnvironmentError('NUPACKHOME environment variable must be set')


class Nupack(object):
    '''Contain sequence inputs and use NUPACK computation methods.'''

    def __init__(self, dna_list):
        '''
        :param dna_list: Input sequences. Can also be list of sequences.
        :type dna_list: str

        '''

        # TODO: automatically figure out material based on input
        # Hard-coded to 'dna' for now.
        self.material = 'dna'

        # Set up nupack environment variable
        self.nupack_home = NUPACKHOME

        # Hard-coded to 'dna' for now.
        self.material = 'dna'
        if type(dna_list) != list:
            dna_list = [dna_list]
        else:
            dna_list = dna_list

        for seq in dna_list:
            check_instance(seq)

        # convert dna object to string
        dna_list = [str(dna) for dna in dna_list]
        self.dna_list = dna_list
        self.outdir = mkdtemp()

        # Track whether complexes has been run
        self.complexes_run = False

    def complexes(self, max_complexes, mfe=True, temp=50):
        '''
        Run \'complexes\'.

        :param max_complexes: Maximum complex size (integer).
        :type max_complexes: int
        :param mfe: Include mfe calculations (boolean, defaults to True).
        :type mfe: bool
        :param temp: Sets the temperature.
        :type temp: float

        '''

        self._open()  # Recreate temp dir if it was removed

        # Prepare input file
        complexes_input = '\n'.join([str(len(self.dna_list)),
                                     '\n'.join(self.dna_list),
                                     str(max_complexes)])

        handle = open(self.outdir + '/nupack.in', 'w')
        handle.write(complexes_input)
        handle.close()

        # Prepare arguments
        if mfe:
            mfe = ' -mfe '
        else:
            mfe = ''

        # Run 'complexes'
        args = ' -T {} -material {} {}'.format(temp, self.material, mfe)
        self._run_nupack('complexes', args)

        # Parse the output
        eq_results = open(self.outdir + '/nupack.cx', 'r+').readlines()
        eq_results = [x.split() for x in eq_results if '%' not in x]

        # Remove rank entry
        for i in eq_results:
            i.pop(0)

        # Extract complexes
        eq_results_cx = []
        for result in eq_results:
            cx_i = []
            for i in range(len(self.dna_list)):
                cx_i.append(int(float(result.pop(0))))
            eq_results_cx.append(cx_i)

        # Extract energies
        eq_results_en = [float(x.pop(0)) for x in eq_results]

        self.complexes_run = max_complexes

        return {'complexes': eq_results_cx,
                'complex_energy': eq_results_en}

    def concentrations(self, max_complexes, conc=0.5e-6, mfe=True):
        '''
        Run \'concentrations\'.

        :param max_complexes: Maximum complex size.
        :type max_complexes: int
        :param conc: Oligo concentrations.
        :type conc: float
        :param mfe: Include mfe calculations.
        :type mfe: bool

        '''

        if self.complexes_run and self.complexes_run == max_complexes:
            pass  # complexes has already been run
        else:
            self._open()  # Recreate temp dir if it was removed
            # Run complexes in order to get input file
            self.complexes(max_complexes=max_complexes,
                           mfe=mfe)

        # Prepare input file
        if type(conc) == list:
            if len(conc) != len(self.dna_list):
                raise ValueError("len(conc) must be len(sequences)")
            concstring = '\n'.join(conc)
        else:
            concstring = '\n'.join([str(conc) for x in self.dna_list])

        handle = open(self.outdir + '/nupack.con', 'w')
        handle.write(concstring)
        handle.close()

        # Run 'complexes'
        self._run_nupack('concentrations', '-sort 3 ')

        # Parse the output of 'complexes'
        eq_results = open(self.outdir + '/nupack.eq', 'r+').readlines()
        eq_results = [x.split() for x in eq_results if '%' not in x]

        # Format results: complex type, concentration, and energy
        # Remove rank information
        for i in eq_results:
            i.pop(0)

        eq_results_cx = []
        for result in eq_results:
            eq_cx_i = []
            for i in range(len(self.dna_list)):
                eq_cx_i.append(int(float(result.pop(0))))
            eq_results_cx.append(eq_cx_i)

        eq_results_en = [x.pop(0) for x in eq_results]
        eq_results_conc = [float(x[0]) for x in eq_results]

        return {'types': eq_results_cx,
                'concentration': eq_results_conc,
                'energy': eq_results_en}

    def mfe(self, index=0, temp=50):
        '''
        Run 'mfe'.

        :param index: Index of strand to analyze.
        :tyep index: int
        :param temp: Temperature in degrees C.
        :type temp: float

        '''
        # TODO: should return multiple mfe results by unless index is specified

        self._open()  # Recreate temp dir if it was removed

        sequence = self.dna_list[index]
        # Prepare input file
        handle = open(self.outdir + '/nupack.in', 'w')
        handle.write(sequence)
        handle.close()

        # Run 'mfe'
        args = ' -T {} -material {}'.format(temp, self.material)
        self._run_nupack('mfe', args)

        # Parse the output of 'mfe'
        mfe_raw = open(self.outdir + '/nupack.mfe', 'r+').readlines()[14]
        mfe = float(mfe_raw.strip())

        # Return the mfe
        return mfe

    def pairs(self, strand, temp=50):
        '''
        Run 'pairs'.

        :param temp: Sets the temperature in degrees C.
        :type temp: float

        '''

        self._open()  # Recreate temp dir if it was removed

        sequence = self.dna_list[strand]
        handle = open(self.outdir + '/nupack.in', 'w')
        handle.write(sequence)
        handle.close()

        # Calculate pair probabilities with 'pairs'
        args = '-T {} -material {}'.format(temp, self.material)
        self._run_nupack('pairs', args)

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
        pair_probabilities = [float(v.split()[2]) for v in pairs]

        # Return the pair probabilities
        return {'type': types, 'probabilities': pair_probabilities}

    def _open(self):
        '''Check for temp dir. If it doesn't exist, create it.'''

        if not isdir(self.outdir):
            self.outdir = mkdtemp()

    def close(self):
        '''Delete the temp dir. This prevents filling up /tmp.'''

        rmtree(self.outdir)

    def _run_nupack(self, cmd, arguments):
        '''Wrapper for running NUPACK commands.'''

        known = ['complexes', 'concentrations', 'mfe', 'pairs']
        if cmd not in known:
            raise ValueError('Command must be one of the following:' + known)

        env_line = "export NUPACKHOME=" + self.nupack_home + ' && '
        command_prefix = '{0}/bin/{1} '.format(self.nupack_home, cmd)
        command_line = 'cd {} && '.format(self.outdir) + command_prefix
        arguments_line = arguments + ' ' + self.outdir + '/nupack'

        run_command = env_line + command_line + arguments_line
        call(run_command + ' > /dev/null', shell=True)


def nupack_multiprocessing(inputs, material, cmd, arguments, report=True):
    '''
    Provides access to NUPACK commands with multiprocessing support.

    :param inputs: List of sequences.
    :type inpus: list
    :param material: Input material: 'dna' or 'rna'.
    :type material: str
    :param cmd: Command: 'mfe', 'pairs', 'complexes', or 'concentrations'.
    :type cmd: str
    :param arguments: Arguments for the command.
    :type arguments: str

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
        passed = 4
        while report:
            completed = nupack_iterator._index
            if (completed == total):
                break
            else:
                if passed >= 4:
                    print '({0}/{1}) '.format(completed, total) + msg
                    passed = 0
                passed += 1
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
    '''
    Create Nupack instance, run command with arguments.

    :param kwargs: keyword arguments to pass to Nupack

    '''

    run = Nupack(kwargs['seq'])
    output = getattr(run, kwargs['cmd'])(**kwargs['arguments'])
    run.close()
    return output
