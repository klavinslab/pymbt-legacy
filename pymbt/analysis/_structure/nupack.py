'''Wrapper for NUPACK 3.0.'''
import multiprocessing
import time
from subprocess import Popen, PIPE
from tempfile import mkdtemp
from shutil import rmtree
from os.path import isdir
from os import environ
from pymbt.analysis.utils import sequence_type


class Nupack(object):
    '''Use several NUPACK functions on a set of input sequences.'''
    def __init__(self, seq_list, rna1999=False, temp=50, nupack_home=None):
        '''
        :param seq_list: Input sequence(s).
        :type seq_list: pymbt.sequence.DNA, pymbt.sequence.RNA, or list of
                        either.
        :param rna1999: Use RNA 1999 settings.
        :type rna1999: bool
        :param temp: Temperature in C.
        :type temp: float
        :param nupack_home: NUPACK home dir. The script attempts to find the
                            NUPACKHOME environment variable if this isn't set.
        :type nupack_home: str

        '''
        # Set up nupack environment variable
        if nupack_home is None:
            if 'NUPACKHOME' in environ:
                nupack_home = environ['NUPACKHOME']
            else:
                msg1 = 'Must have NUPACKHOME environment variable set or '
                msg2 = 'set nupack_home argument.'
                raise ValueError(msg1 + msg2)
        self._nupack_home = nupack_home

        # If input isn't list, make it one
        if type(seq_list) != list:
            self._seq_list = [seq_list]
        else:
            self._seq_list = seq_list

        # Figure out material based on input and ensure it's consistent
        self._material = sequence_type(self._seq_list[0])
        if not all([sequence_type(seq) == self._material for seq in
                   self._seq_list]):
            raise ValueError('Sequence inputs were of mixed types.')

        # Convert seq object(s) to string(s)
        self._seq_list = [str(seq) for seq in self._seq_list]

        # Shared temperature input
        self._temp = temp

        # Handle rna1999 setting
        if self._material == 'rna' and rna1999:
            self._material = 'rna1999'
            if self._temp != 37:
                raise Exception('To use rna1999 settings, temp must be 37.')

        # Create temp dir
        self._tempdir = mkdtemp()

        # Track whether complexes has been run to avoid redundant computation
        self._complexes_run = False
        self._complexes_file = None

    def complexes(self, max_complexes, mfe=True):
        '''Find properties of polymer complexes.

        :param max_complexes: Maximum complex size (integer).
        :type max_complexes: int
        :param mfe: Include mfe calculations (boolean, defaults to True).
        :type mfe: bool

        '''
        self._temp_dir()

        # Prepare input file
        n_seqs = str(len(self._seq_list))
        seqs = '\n'.join(self._seq_list)
        max_complexes = str(max_complexes)
        complexes_input = '\n'.join([n_seqs, seqs, max_complexes])
        with open(self._tempdir + '/nupack.in', 'w') as input_handle:
            input_handle.write(complexes_input)

        # Run 'complexes'
        args = ['-T', str(self._temp), '-material', self._material]
        if mfe:
            args += ['-mfe']
        self._run_cmd('complexes', args)

        # Parse the output
        with open(self._tempdir + '/nupack.cx', 'r+') as output:
            self._complexes_file = output.readlines()
        complexes_results = [x.split() for x in self._complexes_file if '%'
                             not in x]

        # Remove rank entry
        for complexes_result in complexes_results:
            complexes_result.pop(0)

        # Extract complexes
        complexes = []
        for result in complexes_results:
            complex_i = []
            for seq in self._seq_list:
                complex_i.append(int(result.pop(0)))
            complexes.append(complex_i)

        # Extract energies
        energies = [float(x.pop(0)) for x in complexes_results]

        self._complexes_run = int(max_complexes), mfe

        self._close()

        return {'complexes': complexes, 'complex_energy': energies}

    def concentrations(self, max_complexes, conc=[0.5e-6], mfe=True):
        '''Find the predicted concentrations of polymer complexes.

        :param max_complexes: Maximum complex size.
        :type max_complexes: int
        :param conc: Oligo concentrations.
        :type conc: list
        :param mfe: Include mfe calculations.
        :type mfe: bool

        '''
        # If complexes hasn't been run, run it
        if self._complexes_run != (max_complexes, mfe):
            self.complexes(max_complexes=max_complexes, mfe=mfe)

        self._temp_dir()
        with open(self._tempdir + '/nupack.cx', 'w') as cx_file:
            cx_file.writelines(self._complexes_file)

        # Prepare input file
        if len(conc) > 1:
            input_concs = '\n'.join([str(x) for x in conc])
        else:
            input_concs = '\n'.join([str(conc[0]) for x in self._seq_list])
        with open(self._tempdir + '/nupack.con', 'w') as input_handle:
            input_handle.write(input_concs)

        args = ['-sort', str(3)]

        # Run 'concentrations'
        self._run_cmd('concentrations', args)

        # Parse the output of 'complexes'
        with open(self._tempdir + '/nupack.eq', 'r+') as output:
            con_results = output.readlines()
            con_results = [x.split() for x in con_results if '%' not in x]

        # Format results: complex type, concentration, and energy
        # Remove rank information
        for concentration in con_results:
            concentration.pop(0)
        con_types = []
        for result in con_results:
            eq_cx_i = []
            for i in range(len(self._seq_list)):
                eq_cx_i.append(int(float(result.pop(0))))
            con_types.append(eq_cx_i)

        # Extract energies and concentrations
        energies = [x.pop(0) for x in con_results]
        concentrations = [float(x[0]) for x in con_results]

        self._close()

        return {'types': con_types,
                'concentrations': concentrations,
                'energy': energies}

    def mfe(self, index=0):
        '''Calculate the minimum free energy of a single polymer.

        :param index: Index of strand to analyze.
        :type index: int

        '''
        self._temp_dir()

        # Prepare input file
        with open(self._tempdir + '/nupack.in', 'w') as input_handle:
            input_handle.write(self._seq_list[index])

        args = ['-T', str(self._temp), '-material', self._material]
        self._run_cmd('mfe', args)
        # Parse the output of 'mfe'

        with open('{}/nupack.mfe'.format(self._tempdir), 'r+') as output:
            mfe = float(output.readlines()[14].strip())

        self._close()

        return mfe

    def pairs(self, index=0):
        '''Calculate unbound pair probabilities.

        :param index: Index of strand to analyze.
        :type index: int

        '''
        self._temp_dir()

        # Sequence at the specified index
        sequence = self._seq_list[index]

        # Prepare input file
        with open(self._tempdir + '/nupack.in', 'w') as input_handle:
            input_handle.write(sequence)

        # Calculate pair probabilities with 'pairs'
        args = ['-T', str(self._temp), '-material', self._material]
        self._run_cmd('pairs', args)

        # Parse the output of 'pairs'
        # Only look at the last n rows - unbound probabilities
        with open(self._tempdir + '/nupack.ppairs', 'r+') as handle:
            pairs = handle.readlines()[-len(sequence):]
            pairs = [pair for pair in pairs if '%' not in pair]

        # Extract pair probability types and pair_probabilities. These are in
        # Nupack's raw text format
        #types = [(int(x.split()[0]), int(x.split()[1])) for x in pairs]
        pair_probabilities = [float(x.split()[2]) for x in pairs]

        self._close()

        return pair_probabilities

    def _close(self):
        '''Delete the temporary dir to prevent filling up /tmp.'''
        rmtree(self._tempdir)

    def _temp_dir(self):
        '''If temporary dir doesn't exist, create it.'''
        if not isdir(self._tempdir):
            self._tempdir = mkdtemp()

    def _run_cmd(self, cmd, cmd_args):
        '''Run NUPACK command line programs.

        :param cmd: NUPACK command line tool to run ('mfe', 'complexes',
                    'concentrations', 'pairs').
        :type cmd: str
        :param cmd_args: Arguments to pass to the command line.
        :type cmd_args: str

        '''
        known_cmds = ['complexes', 'concentrations', 'mfe', 'pairs']
        if cmd not in known_cmds:
            msg = 'Command must be one of: {}'.format(known_cmds)
            raise ValueError(msg)

        cmd_path = '{0}/bin/{1}'.format(self._nupack_home, cmd)
        cmd_args += ['nupack']
        cmd_list = [cmd_path] + cmd_args
        process = Popen(cmd_list,
                        env={'NUPACKHOME': self._nupack_home},
                        cwd=str(self._tempdir), stdout=PIPE)
        process.wait()


def nupack_multiprocessing(seqs, material, cmd, arguments, report=True):
    '''Run NUPACK commands with the benefits of multiprocessing.

    :param inputs: List of sequences, same format as for pymbt.analysis.Nupack.
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
        args = [{'seq': seq,
                 'cmd': cmd,
                 'material': material,
                 'arguments': arguments} for seq in seqs]
        nupack_iterator = nupack_pool.imap(run_nupack, args)
        total = len(seqs)
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
        nupack_pool.terminate()
        nupack_pool.close()
        raise KeyboardInterrupt

    return multi_output


def run_nupack(kwargs):
    '''Run Nupack command in a picklable way.

    :param kwargs: keyword arguments to pass to Nupack as well as 'cmd'.

    '''
    run = Nupack(kwargs['seq'])
    output = getattr(run, kwargs['cmd'])(**kwargs['arguments'])
    return output
