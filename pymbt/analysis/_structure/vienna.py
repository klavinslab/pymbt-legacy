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
        self.seq = seq
        self.temp = temp
        self.tmpdir = mkdtemp()

    def mfe(self):
        '''Calculate minimum free energy of sequence.'''
        process = Popen(['RNAfold', '-T', str(self.temp)], stdin=PIPE,
                        stdout=PIPE, stderr=STDOUT, cwd=self.tmpdir)
        output = process.communicate(input=self.seq)[0]
        lines = output.splitlines()
        lines = lines[-1].split('(')[-1].split(')')[0].strip()
        mfe = float(lines)
        return mfe

    def pairs(self):
        '''Calculate pair probabilities.'''
        process = Popen(['RNAfold', '-p', '-T', str(self.temp)], stdin=PIPE,
                        stdout=PIPE, stderr=STDOUT, cwd=self.tmpdir)
        process.communicate(input=self.seq)[0]
        with open('{}/dot.ps'.format(self.tmpdir), 'r') as dot:
            text = dot.read()
            text = text[text.find('%data'):]
        split = text.split('\n')
        data = [x for x in split if x.endswith('ubox')]
        data = [x.rstrip(' ubox') for x in data]
        data = [x.split() for x in data]
        data = [(int(a), int(b), float(c)) for a, b, c in data]
        data = sorted(data, key=lambda triple: triple[2], reverse=True)
        return data

    def close(self):
        '''Close the temporary dir (keeps /tmp clean).'''
        rmtree(self.tmpdir)

    def _check_tmpdir(self):
        '''If temp dir has been removed, create a new one.'''
        if not isdir(self.tmpdir):
            self.tmpdir = mkdtemp()
