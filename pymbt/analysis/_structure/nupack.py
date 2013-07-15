'''
Wrapper for NUPACK 3.0.

'''

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
        if not nupack_home:
            if 'NUPACKHOME' in environ:
                nupack_home = environ['NUPACKHOME']
            else:
                msg1 = 'Must have NUPACKHOME environment variable set or '
                msg2 = 'set nupack_home argument.'
                raise ValueError(msg1 + msg2)
        self.nupack_home = nupack_home

        # If input isn't list, make it one
        if type(seq_list) != list:
            self.seq_list = [seq_list]
        else:
            self.seq_list = seq_list

        # Figure out material based on input and ensure it's consistent
        self.material = sequence_type(self.seq_list[0])
        if not all([sequence_type(seq) == self.material for seq in
                   self.seq_list]):
            raise ValueError('Sequence inputs were of mixed types.')

        # Convert seq object(s) to string(s)
        self.seq_list = [str(seq) for seq in self.seq_list]

        # Shared temperature input
        self.temp = temp

        # Handle rna1999 setting
        if self.material == 'rna' and rna1999:
            self.material = 'rna1999'
            if self.temp != 37:
                raise Exception('To use rna1999 settings, temp must be 37.')

        # Create temp dir
        self.tmpdir = mkdtemp()

        # Track whether complexes has been run to avoid redundant computation
        self._complexes_run = False

    def complexes(self, max_complexes, mfe=True):
        '''Find properties of polymer complexes.

        :param max_complexes: Maximum complex size (integer).
        :type max_complexes: int
        :param mfe: Include mfe calculations (boolean, defaults to True).
        :type mfe: bool

        '''
        self._temp_dir()

        # Prepare input file
        n_seqs = str(len(self.seq_list))
        seqs = '\n'.join(self.seq_list)
        max_complexes = str(max_complexes)
        complexes_input = '\n'.join([n_seqs, seqs, max_complexes])
        with open(self.tmpdir + '/nupack.in', 'w') as input_handle:
            input_handle.write(complexes_input)

        # Run 'complexes'
        args = ['-T', str(self.temp), '-material', self.material]
        if mfe:
            args += ['-mfe']
        self._run_cmd('complexes', args)

        # Parse the output
        with open(self.tmpdir + '/nupack.cx', 'r+') as output_handle:
            complexes_results = output_handle.readlines()
            complexes_results = [x.split() for x in complexes_results if '%'
                                 not in x]

        # Remove rank entry
        for complexes_result in complexes_results:
            complexes_result.pop(0)

        # Extract complexes
        complexes = []
        for result in complexes_results:
            complex_i = []
            for seq in self.seq_list:
                complex_i.append(int(float(result.pop(0))))
            complexes.append(complex_i)

        # Extract energies
        energies = [float(x.pop(0)) for x in complexes_results]

        self._complexes_run = int(max_complexes), mfe

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

        self._temp_dir()

        # If complexes has already been run with the same settings, keep the
        # result (more efficient). Otherwise, run complexes.
        if self._complexes_run == (max_complexes, mfe):
            pass
        else:
            self.complexes(max_complexes=max_complexes, mfe=mfe)

        # Prepare input file
        if len(conc) > 1:
            input_concs = '\n'.join([str(x) for x in conc])
        else:
            input_concs = '\n'.join([str(conc[0]) for x in self.seq_list])
        with open(self.tmpdir + '/nupack.con', 'w') as input_handle:
            input_handle.write(input_concs)

        args = ['-sort', str(3)]

        # Run 'concentrations'
        self._run_cmd('concentrations', args)

        # Parse the output of 'complexes'
        with open(self.tmpdir + '/nupack.eq', 'r+') as output_handle:
            con_results = output_handle.readlines()
            con_results = [x.split() for x in con_results if '%' not in x]

        # Format results: complex type, concentration, and energy
        # Remove rank information
        for concentration in con_results:
            concentration.pop(0)
        con_types = []
        for result in con_results:
            eq_cx_i = []
            for i in range(len(self.seq_list)):
                eq_cx_i.append(int(float(result.pop(0))))
            con_types.append(eq_cx_i)

        # Extract energies and concentrations
        energies = [x.pop(0) for x in con_results]
        concentrations = [float(x[0]) for x in con_results]

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
        with open(self.tmpdir + '/nupack.in', 'w') as input_handle:
            input_handle.write(self.seq_list[index])

        args = ['-T', str(self.temp), '-material', self.material]
        self._run_cmd('mfe', args)
        # Parse the output of 'mfe'

        with open('{}/nupack.mfe'.format(self.tmpdir), 'r+') as output_handle:
            mfe = float(output_handle.readlines()[14].strip())
        # Return the mfe

        return mfe

    def pairs(self, index=0):
        '''Calculate unbound pair probabilities.

        :param index: Index of strand to analyze.
        :type index: int

        '''

        self._temp_dir()

        # Sequence at the specified index
        sequence = self.seq_list[index]

        # Prepare input file
        with open(self.tmpdir + '/nupack.in', 'w') as input_handle:
            input_handle.write(sequence)

        # Calculate pair probabilities with 'pairs'
        args = ['-T', str(self.temp), '-material', self.material]
        self._run_cmd('pairs', args)

        # Parse the output of 'pairs'
        # Only look at the last n rows - unbound probabilities
        # TODO: return both unbound and bound under separate keys
        with open(self.tmpdir + '/nupack.ppairs', 'r+') as handle:
            pairs = handle.readlines()[-len(sequence):]
            pairs = [pair for pair in pairs if '%' not in pair]

        # Extract pair probability types and pair_probabilities. These are in
        # Nupack's raw text format
        #types = [(int(x.split()[0]), int(x.split()[1])) for x in pairs]
        pair_probabilities = [float(x.split()[2]) for x in pairs]

        # Return the pair probabilities
        return pair_probabilities
#        return {'type': types, 'probabilities': pair_probabilities}

    def close(self):
        '''Delete the temporary dir to prevent filling up /tmp.'''
        rmtree(self.tmpdir)

    def _temp_dir(self):
        '''If temporary dir doesn't exist, create it.'''
        if not isdir(self.tmpdir):
            self.tmpdir = mkdtemp()

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

        cmd_path = '{0}/bin/{1}'.format(self.nupack_home, cmd)
        cmd_args += ['nupack']
        cmd_list = [cmd_path] + cmd_args
        process = Popen(cmd_list,
                        env={'NUPACKHOME': self.nupack_home},
                        cwd=str(self.tmpdir), stdout=PIPE)
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
    run.close()
    return output
