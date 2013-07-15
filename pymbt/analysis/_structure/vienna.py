'''Vienna RNA module.'''
from subprocess import Popen, PIPE, STDOUT


class Vienna(object):
    '''Vienna functions.'''
    def __init__(self, seq, temp=50):
        self.seq = seq
        self.temp = temp

    def mfe(self):
        process = Popen(['RNAfold', '-T', str(self.temp)], stdin=PIPE,
                        stdout=PIPE,
                        stderr=STDOUT)
        output = process.communicate(input=self.seq)[0]
        lines = output.splitlines()
        lines = lines[-1].split('(')[-1].split(')')[0].strip()
        mfe = float(lines)
        return mfe
