'''Vienna RNA module.'''
from subprocess import Popen, PIPE, STDOUT
from tempfile import mkdtemp
from os.path import isdir
from shutil import rmtree


class Vienna(object):
    '''Vienna functions.'''
    def __init__(self, seq, temp=50.0):
        '''
        :param seq: DNA or RNA sequence to evaluate.
        :type seq: pymbt.sequence.DNA or pymbt.sequence.RNA
        :param temp: Temperature at which to run calculations.
        :type temp: float
        '''
        self._seq = str(seq)
        self._temp = temp
        self._tempdir = mkdtemp()

    def mfe(self):
        '''Calculate minimum free energy of sequence.'''
        self._check_tempdir()
        process = Popen(['RNAfold', '-T', str(self._temp)], stdin=PIPE,
                        stdout=PIPE, stderr=STDOUT, cwd=self._tempdir)
        output = process.communicate(input=self._seq)[0]
        lines = output.splitlines()
        lines = lines[-1].split('(')[-1].split(')')[0].strip()
        mfe = float(lines)
        self._close()
        return mfe

    def pairs(self):
        '''Calculate pair probabilities.'''
        self._check_tempdir()
        process = Popen(['RNAfold', '-p', '-T', str(self._temp)], stdin=PIPE,
                        stdout=PIPE, stderr=STDOUT, cwd=self._tempdir)
        process.communicate(input=self._seq)[0]
        with open('{}/dot.ps'.format(self._tempdir), 'r') as dot:
            text = dot.read()
            text = text[text.find('%data'):]
        split = text.split('\n')
        data = [x for x in split if x.endswith('ubox')]
        data = [x.rstrip(' ubox') for x in data]
        data = [x.split() for x in data]
        data = [(int(a), int(b), float(c)) for a, b, c in data]
        unbound = [1.0] * len(self._seq)
        # TODO: this is not the right way to calculate the unbound
        # probability
        for base1, base2, prob_sqr in data:
            probability = prob_sqr**2
            unbound[base1 - 1] -= probability
            unbound[base2 - 1] -= probability
        self._close()
        return unbound

    def _close(self):
        '''Close the temporary dir (keeps /tmp clean).'''
        rmtree(self._tempdir)

    def _check_tempdir(self):
        '''If temp dir has been removed, create a new one.'''
        if not isdir(self._tempdir):
            self._tempdir = mkdtemp()
